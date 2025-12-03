//! PCISPH - Predictive-Corrective Incompressible SPH
//!
//! Implementation of the PCISPH algorithm from Solenthaler & Pajarola (2009).
//! This solver maintains incompressibility through iterative pressure correction.
//!
//! ## Algorithm Overview
//!
//! 1. Compute non-pressure forces (gravity, viscosity, external)
//! 2. Predict positions based on non-pressure forces
//! 3. Iteratively correct pressure until density error is below threshold
//! 4. Update velocities and positions
//!
//! ## References
//!
//! - Solenthaler & Pajarola (2009): Predictive-Corrective Incompressible SPH

use crate::config::PhysicsConstants;
use crate::kernels::{grad_wspiky, lap_wviscosity, wpoly6};
use crate::particle::{ParticleBuffer, ParticleType};
use crate::spatial_hash::SpatialHash;
use rayon::prelude::*;

/// PCISPH Solver
pub struct PcisphSolver {
    /// Physics constants
    physics: PhysicsConstants,
    /// Spatial hash for neighbor search
    spatial_hash: SpatialHash,
    /// Predicted positions buffer
    predicted_positions: Vec<f32>,
    /// Pressure forces buffer
    pressure_forces: Vec<f32>,
    /// Predicted densities buffer
    predicted_densities: Vec<f32>,
}

impl PcisphSolver {
    /// Create a new PCISPH solver
    pub fn new(physics: PhysicsConstants, spatial_hash: SpatialHash, max_particles: usize) -> Self {
        Self {
            physics,
            spatial_hash,
            predicted_positions: vec![0.0; max_particles * 4],
            pressure_forces: vec![0.0; max_particles * 4],
            predicted_densities: vec![0.0; max_particles],
        }
    }

    /// Perform one simulation step
    pub fn step(&mut self, particles: &mut ParticleBuffer, dt: f32) {
        let n = particles.len();
        if n == 0 {
            return;
        }

        // Ensure buffers are large enough
        self.ensure_capacity(n);

        // 1. Build spatial hash and find neighbors
        self.spatial_hash.hash_particles(particles);
        self.spatial_hash.sort_particles(particles);
        self.spatial_hash.build_cell_index(particles);

        // 2. Compute density
        self.compute_density(particles);

        // 3. Compute non-pressure forces (gravity, viscosity)
        self.compute_non_pressure_forces(particles, dt);

        // 4. PCISPH pressure correction loop
        self.pressure_correction_loop(particles, dt);

        // 5. Integrate (update velocities and positions)
        self.integrate(particles, dt);
    }

    /// Ensure internal buffers have enough capacity
    fn ensure_capacity(&mut self, n: usize) {
        if self.predicted_positions.len() < n * 4 {
            self.predicted_positions.resize(n * 4, 0.0);
            self.pressure_forces.resize(n * 4, 0.0);
            self.predicted_densities.resize(n, 0.0);
        }
    }

    /// Compute density for all particles
    pub fn compute_density(&self, particles: &mut ParticleBuffer) {
        let n = particles.len();
        let h = self.physics.h_scaled();
        let h_sq = h * h;
        let mass = self.physics.mass;
        let wpoly6_coeff = self.physics.wpoly6_coeff as f32;

        // Parallel density computation
        let densities: Vec<f32> = (0..n)
            .into_par_iter()
            .map(|i| {
                let pos_i = [
                    particles.sorted_positions[i * 4],
                    particles.sorted_positions[i * 4 + 1],
                    particles.sorted_positions[i * 4 + 2],
                ];

                let mut density = 0.0;

                // Self contribution
                density += mass * wpoly6(0.0, h, h_sq, wpoly6_coeff);

                // Get neighbors from spatial hash
                for cell_id in self.spatial_hash.neighbor_cells(pos_i) {
                    let cell_idx = cell_id as usize;
                    if cell_idx >= self.spatial_hash.num_cells() {
                        continue;
                    }

                    let start = particles.grid_cell_index[cell_idx];
                    let end = particles.grid_cell_index_end[cell_idx];

                    if start == u32::MAX {
                        continue;
                    }

                    for j in start..end {
                        let j = j as usize;
                        if j == i {
                            continue;
                        }

                        let pos_j = [
                            particles.sorted_positions[j * 4],
                            particles.sorted_positions[j * 4 + 1],
                            particles.sorted_positions[j * 4 + 2],
                        ];

                        let dx = pos_i[0] - pos_j[0];
                        let dy = pos_i[1] - pos_j[1];
                        let dz = pos_i[2] - pos_j[2];
                        let r_sq = dx * dx + dy * dy + dz * dz;

                        if r_sq < h_sq {
                            density += mass * wpoly6(r_sq, h, h_sq, wpoly6_coeff);
                        }
                    }
                }

                density
            })
            .collect();

        // Copy back to particle buffer
        for (i, d) in densities.iter().enumerate() {
            particles.densities[i] = *d;
        }
    }

    /// Compute non-pressure forces (gravity, viscosity)
    fn compute_non_pressure_forces(&self, particles: &mut ParticleBuffer, _dt: f32) {
        let n = particles.len();
        let h = self.physics.h_scaled();
        let viscosity = self.physics.viscosity;
        let gravity = self.physics.gravity;
        let mass = self.physics.mass;
        let lap_coeff = self.physics.lap_wviscosity_coeff as f32;

        // Reset accelerations
        for i in 0..n {
            let base = i * 4;
            particles.accelerations[base] = 0.0;
            particles.accelerations[base + 1] = 0.0;
            particles.accelerations[base + 2] = 0.0;
        }

        // Parallel force computation
        let forces: Vec<[f32; 3]> = (0..n)
            .into_par_iter()
            .map(|i| {
                let ptype_i = particles.types[i];
                if !ptype_i.is_dynamic() {
                    return [0.0, 0.0, 0.0];
                }

                let pos_i = [
                    particles.sorted_positions[i * 4],
                    particles.sorted_positions[i * 4 + 1],
                    particles.sorted_positions[i * 4 + 2],
                ];
                let vel_i = [
                    particles.sorted_velocities[i * 4],
                    particles.sorted_velocities[i * 4 + 1],
                    particles.sorted_velocities[i * 4 + 2],
                ];

                let density_i = particles.densities[i];
                if density_i < 1e-10 {
                    return gravity;
                }

                let mut force = [gravity[0], gravity[1], gravity[2]];

                // Viscosity force
                for cell_id in self.spatial_hash.neighbor_cells(pos_i) {
                    let cell_idx = cell_id as usize;
                    if cell_idx >= self.spatial_hash.num_cells() {
                        continue;
                    }

                    let start = particles.grid_cell_index[cell_idx];
                    let end = particles.grid_cell_index_end[cell_idx];

                    if start == u32::MAX {
                        continue;
                    }

                    for j in start..end {
                        let j = j as usize;
                        if j == i {
                            continue;
                        }

                        let pos_j = [
                            particles.sorted_positions[j * 4],
                            particles.sorted_positions[j * 4 + 1],
                            particles.sorted_positions[j * 4 + 2],
                        ];

                        let dx = pos_i[0] - pos_j[0];
                        let dy = pos_i[1] - pos_j[1];
                        let dz = pos_i[2] - pos_j[2];
                        let r = (dx * dx + dy * dy + dz * dz).sqrt();

                        if r < h && r > 1e-10 {
                            let vel_j = [
                                particles.sorted_velocities[j * 4],
                                particles.sorted_velocities[j * 4 + 1],
                                particles.sorted_velocities[j * 4 + 2],
                            ];
                            let density_j = particles.densities[j];

                            if density_j > 1e-10 {
                                let lap = lap_wviscosity(r, h, lap_coeff);

                                // Viscosity: μ * m * (v_j - v_i) / ρ_j * ∇²W
                                let factor = viscosity * mass / density_j * lap;
                                force[0] += factor * (vel_j[0] - vel_i[0]);
                                force[1] += factor * (vel_j[1] - vel_i[1]);
                                force[2] += factor * (vel_j[2] - vel_i[2]);
                            }
                        }
                    }
                }

                force
            })
            .collect();

        // Copy forces to accelerations
        for (i, f) in forces.iter().enumerate() {
            let base = i * 4;
            particles.accelerations[base] = f[0];
            particles.accelerations[base + 1] = f[1];
            particles.accelerations[base + 2] = f[2];
        }
    }

    /// PCISPH pressure correction loop
    fn pressure_correction_loop(&mut self, particles: &mut ParticleBuffer, dt: f32) {
        let n = particles.len();
        let h = self.physics.h_scaled();
        let h_sq = h * h;
        let rho0 = self.physics.rho0;
        let mass = self.physics.mass;
        let grad_coeff = self.physics.grad_wspiky_coeff as f32;

        // Initialize pressures to zero
        for i in 0..n {
            particles.pressures[i] = 0.0;
        }

        // Predict initial positions
        for i in 0..n {
            let base = i * 4;
            self.predicted_positions[base] = particles.sorted_positions[base]
                + dt * particles.sorted_velocities[base]
                + dt * dt * particles.accelerations[base];
            self.predicted_positions[base + 1] = particles.sorted_positions[base + 1]
                + dt * particles.sorted_velocities[base + 1]
                + dt * dt * particles.accelerations[base + 1];
            self.predicted_positions[base + 2] = particles.sorted_positions[base + 2]
                + dt * particles.sorted_velocities[base + 2]
                + dt * dt * particles.accelerations[base + 2];
        }

        // Correction iterations
        for _iter in 0..self.physics.max_iterations {
            // Predict density at predicted positions
            self.predict_density(particles, h, h_sq, mass);

            // Compute density error and update pressure
            let mut max_error: f32 = 0.0;
            for i in 0..n {
                if particles.types[i].has_pressure() {
                    let density_error = self.predicted_densities[i] - rho0;
                    max_error = max_error.max(density_error.abs() / rho0);

                    // Update pressure: p += δ * error
                    // δ is derived from the discretized pressure-density relation
                    let delta = self.physics.beta as f32 / (dt * dt);
                    particles.pressures[i] += delta * density_error.max(0.0);
                }
            }

            // Check convergence
            if max_error < self.physics.density_error_threshold {
                break;
            }

            // Compute pressure force
            self.compute_pressure_force(particles, h, grad_coeff);

            // Update predicted positions
            for i in 0..n {
                if particles.types[i].is_dynamic() {
                    let base = i * 4;
                    let density_i = particles.densities[i];
                    if density_i > 1e-10 {
                        let factor = dt * dt / density_i;
                        self.predicted_positions[base] += factor * self.pressure_forces[base];
                        self.predicted_positions[base + 1] += factor * self.pressure_forces[base + 1];
                        self.predicted_positions[base + 2] += factor * self.pressure_forces[base + 2];
                    }
                }
            }
        }
    }

    /// Predict density at predicted positions
    fn predict_density(&mut self, particles: &ParticleBuffer, h: f32, h_sq: f32, mass: f32) {
        let n = particles.len();
        let wpoly6_coeff = self.physics.wpoly6_coeff as f32;

        for i in 0..n {
            let pos_i = [
                self.predicted_positions[i * 4],
                self.predicted_positions[i * 4 + 1],
                self.predicted_positions[i * 4 + 2],
            ];

            let mut density = mass * wpoly6(0.0, h, h_sq, wpoly6_coeff);

            for cell_id in self.spatial_hash.neighbor_cells(pos_i) {
                let cell_idx = cell_id as usize;
                if cell_idx >= self.spatial_hash.num_cells() {
                    continue;
                }

                let start = particles.grid_cell_index[cell_idx];
                let end = particles.grid_cell_index_end[cell_idx];

                if start == u32::MAX {
                    continue;
                }

                for j in start..end {
                    let j = j as usize;
                    if j == i {
                        continue;
                    }

                    let pos_j = [
                        self.predicted_positions[j * 4],
                        self.predicted_positions[j * 4 + 1],
                        self.predicted_positions[j * 4 + 2],
                    ];

                    let dx = pos_i[0] - pos_j[0];
                    let dy = pos_i[1] - pos_j[1];
                    let dz = pos_i[2] - pos_j[2];
                    let r_sq = dx * dx + dy * dy + dz * dz;

                    if r_sq < h_sq {
                        density += mass * wpoly6(r_sq, h, h_sq, wpoly6_coeff);
                    }
                }
            }

            self.predicted_densities[i] = density;
        }
    }

    /// Compute pressure force from pressure values
    fn compute_pressure_force(&mut self, particles: &ParticleBuffer, h: f32, grad_coeff: f32) {
        let n = particles.len();
        let mass = self.physics.mass;

        // Reset pressure forces
        for i in 0..n * 4 {
            self.pressure_forces[i] = 0.0;
        }

        for i in 0..n {
            if !particles.types[i].is_dynamic() {
                continue;
            }

            let pos_i = [
                self.predicted_positions[i * 4],
                self.predicted_positions[i * 4 + 1],
                self.predicted_positions[i * 4 + 2],
            ];

            let pressure_i = particles.pressures[i];
            let density_i = particles.densities[i];

            if density_i < 1e-10 {
                continue;
            }

            for cell_id in self.spatial_hash.neighbor_cells(pos_i) {
                let cell_idx = cell_id as usize;
                if cell_idx >= self.spatial_hash.num_cells() {
                    continue;
                }

                let start = particles.grid_cell_index[cell_idx];
                let end = particles.grid_cell_index_end[cell_idx];

                if start == u32::MAX {
                    continue;
                }

                for j in start..end {
                    let j = j as usize;
                    if j == i {
                        continue;
                    }

                    let pos_j = [
                        self.predicted_positions[j * 4],
                        self.predicted_positions[j * 4 + 1],
                        self.predicted_positions[j * 4 + 2],
                    ];

                    let r = [
                        pos_i[0] - pos_j[0],
                        pos_i[1] - pos_j[1],
                        pos_i[2] - pos_j[2],
                    ];
                    let r_len = (r[0] * r[0] + r[1] * r[1] + r[2] * r[2]).sqrt();

                    if r_len < h && r_len > 1e-10 {
                        let pressure_j = particles.pressures[j];
                        let density_j = particles.densities[j];

                        if density_j > 1e-10 {
                            // Pressure force: -m * (p_i/ρ_i² + p_j/ρ_j²) * ∇W
                            let factor = -mass
                                * (pressure_i / (density_i * density_i)
                                    + pressure_j / (density_j * density_j));

                            let grad = grad_wspiky(r, r_len, h, grad_coeff);

                            let base = i * 4;
                            self.pressure_forces[base] += factor * grad[0];
                            self.pressure_forces[base + 1] += factor * grad[1];
                            self.pressure_forces[base + 2] += factor * grad[2];
                        }
                    }
                }
            }
        }
    }

    /// Integrate velocities and positions
    fn integrate(&self, particles: &mut ParticleBuffer, dt: f32) {
        let n = particles.len();

        for i in 0..n {
            if !particles.types[i].is_dynamic() {
                continue;
            }

            let base = i * 4;
            let density_i = particles.densities[i];

            // Total acceleration = non-pressure + pressure
            let mut acc = [
                particles.accelerations[base],
                particles.accelerations[base + 1],
                particles.accelerations[base + 2],
            ];

            if density_i > 1e-10 {
                acc[0] += self.pressure_forces[base] / density_i;
                acc[1] += self.pressure_forces[base + 1] / density_i;
                acc[2] += self.pressure_forces[base + 2] / density_i;
            }

            // Update velocity: v = v + dt * a
            particles.velocities[base] = particles.sorted_velocities[base] + dt * acc[0];
            particles.velocities[base + 1] = particles.sorted_velocities[base + 1] + dt * acc[1];
            particles.velocities[base + 2] = particles.sorted_velocities[base + 2] + dt * acc[2];

            // Update position: x = x + dt * v
            particles.positions[base] = particles.sorted_positions[base]
                + dt * particles.velocities[base];
            particles.positions[base + 1] = particles.sorted_positions[base + 1]
                + dt * particles.velocities[base + 1];
            particles.positions[base + 2] = particles.sorted_positions[base + 2]
                + dt * particles.velocities[base + 2];
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::particle::Particle;

    #[test]
    fn test_density_computation() {
        let physics = PhysicsConstants::water();
        let spatial_hash = SpatialHash::new(
            physics.hash_grid_cell_size,
            [-10.0, -10.0, -10.0],
            [10.0, 10.0, 10.0],
        );
        let solver = PcisphSolver::new(physics.clone(), spatial_hash, 1000);

        let mut particles = ParticleBuffer::with_capacity(100, 32);

        // Add particles in a small cube
        for x in 0..3 {
            for y in 0..3 {
                for z in 0..3 {
                    let p = Particle::new(
                        [x as f32 * 0.5, y as f32 * 0.5, z as f32 * 0.5],
                        [0.0, 0.0, 0.0],
                        ParticleType::Liquid,
                    );
                    particles.push(&p);
                }
            }
        }

        particles.init_sorted_buffers();
        solver.compute_density(&mut particles);

        // All particles should have positive density
        for i in 0..particles.len() {
            assert!(particles.densities[i] > 0.0, "Density should be positive");
        }
    }
}
