//! SPH World - Main Simulation Container
//!
//! The SphWorld struct is the main entry point for SPH simulations.
//! It manages particles, physics, and the simulation loop.

use crate::config::{SphConfig, SimulationBounds};
use crate::elastic::ElasticNetwork;
use crate::membrane::Membrane;
use crate::muscle::MuscleActivation;
use crate::particle::{Particle, ParticleBuffer, ParticleType};
use crate::pcisph::PcisphSolver;
use crate::spatial_hash::SpatialHash;
use crate::{Result, SphError};

/// Main SPH simulation world
pub struct SphWorld {
    /// Configuration
    config: SphConfig,
    /// Particle buffer
    particles: ParticleBuffer,
    /// Elastic connections
    elastic: ElasticNetwork,
    /// Membrane surfaces
    membrane: Membrane,
    /// Current muscle activation
    muscle_activation: MuscleActivation,
    /// PCISPH solver
    solver: PcisphSolver,
    /// Current simulation time
    time: f64,
    /// Current step count
    step: u64,
}

impl SphWorld {
    /// Create a new SPH world with the given configuration
    pub fn new(config: SphConfig) -> Self {
        let bounds = &config.bounds;
        let spatial_hash = SpatialHash::from_config(
            &config.physics,
            [bounds.x_min, bounds.y_min, bounds.z_min],
            [bounds.x_max, bounds.y_max, bounds.z_max],
        );

        let max_particles = config.solver.max_particles;
        let max_neighbors = config.solver.max_neighbors;

        let solver = PcisphSolver::new(
            config.physics.clone(),
            spatial_hash,
            max_particles,
        );

        Self {
            config,
            particles: ParticleBuffer::with_capacity(max_particles, max_neighbors),
            elastic: ElasticNetwork::new(),
            membrane: Membrane::new(),
            muscle_activation: MuscleActivation::new(),
            solver,
            time: 0.0,
            step: 0,
        }
    }

    /// Create a world configured for C. elegans simulation
    pub fn worm() -> Self {
        Self::new(SphConfig::celegans())
    }

    /// Create a world for fluid simulation
    pub fn fluid() -> Self {
        Self::new(SphConfig::fluid())
    }

    /// Add a single particle
    pub fn add_particle(&mut self, position: [f32; 3], velocity: [f32; 3], ptype: ParticleType) -> Result<usize> {
        if self.particles.len() >= self.config.solver.max_particles {
            return Err(SphError::ParticleLimit {
                count: self.particles.len() + 1,
                max: self.config.solver.max_particles,
            });
        }

        let particle = Particle::new(position, velocity, ptype);
        self.particles.push(&particle);
        Ok(self.particles.len() - 1)
    }

    /// Add multiple particles in a grid pattern
    pub fn add_particle_block(
        &mut self,
        min: [f32; 3],
        max: [f32; 3],
        spacing: f32,
        ptype: ParticleType,
    ) -> Result<Vec<usize>> {
        let mut indices = Vec::new();

        let mut x = min[0];
        while x <= max[0] {
            let mut y = min[1];
            while y <= max[1] {
                let mut z = min[2];
                while z <= max[2] {
                    let idx = self.add_particle([x, y, z], [0.0, 0.0, 0.0], ptype)?;
                    indices.push(idx);
                    z += spacing;
                }
                y += spacing;
            }
            x += spacing;
        }

        Ok(indices)
    }

    /// Add an elastic connection between two particles
    pub fn connect(&mut self, p1: usize, p2: usize, stiffness: f32) {
        let pos1 = self.particles.get_position(p1);
        let pos2 = self.particles.get_position(p2);

        let dx = pos2[0] - pos1[0];
        let dy = pos2[1] - pos1[1];
        let dz = pos2[2] - pos1[2];
        let rest_length = (dx * dx + dy * dy + dz * dz).sqrt();

        self.elastic.connect(p1 as u32, p2 as u32, rest_length, stiffness);
    }

    /// Add a muscle connection
    pub fn connect_muscle(&mut self, p1: usize, p2: usize, stiffness: f32, muscle_id: i32) {
        let pos1 = self.particles.get_position(p1);
        let pos2 = self.particles.get_position(p2);

        let dx = pos2[0] - pos1[0];
        let dy = pos2[1] - pos1[1];
        let dz = pos2[2] - pos1[2];
        let rest_length = (dx * dx + dy * dy + dz * dz).sqrt();

        self.elastic.connect_muscle(p1 as u32, p2 as u32, rest_length, stiffness, muscle_id);
    }

    /// Add a membrane triangle
    pub fn add_membrane(&mut self, p1: usize, p2: usize, p3: usize) {
        self.membrane.add(p1 as u32, p2 as u32, p3 as u32);
    }

    /// Set muscle activation pattern
    pub fn set_muscle_activations(&mut self, activation: &MuscleActivation) {
        self.muscle_activation = activation.clone();

        // Update elastic network with new activations
        let flat = activation.to_flat();
        self.elastic.set_all_activations(&flat);
    }

    /// Set single muscle activation by index
    pub fn set_muscle_activation(&mut self, muscle_id: usize, activation: f32) {
        self.muscle_activation.set_by_index(muscle_id, activation);
        self.elastic.set_muscle_activation(muscle_id, activation);
    }

    /// Perform one simulation step
    pub fn step(&mut self) {
        let dt = self.config.physics.time_step;

        // Initialize sorted buffers if needed
        if self.particles.sorted_positions.len() != self.particles.positions.len() {
            self.particles.init_sorted_buffers();
        }

        // Apply elastic forces (including muscles)
        self.elastic.compute_forces(
            &mut self.particles,
            self.config.physics.elasticity,
            self.config.physics.max_muscle_force,
        );

        // Apply membrane forces
        if self.config.solver.membranes_enabled {
            self.membrane.compute_forces(
                &mut self.particles,
                self.config.physics.surface_tension,
            );
        }

        // Run PCISPH solver
        self.solver.step(&mut self.particles, dt);

        // Apply boundary conditions
        self.apply_boundaries();

        // Update time
        self.time += dt as f64;
        self.step += 1;
    }

    /// Run simulation for a given duration
    pub fn run(&mut self, duration: f64) {
        let dt = self.config.physics.time_step as f64;
        let steps = (duration / dt).ceil() as u64;

        for _ in 0..steps {
            self.step();
        }
    }

    /// Apply boundary conditions
    fn apply_boundaries(&mut self) {
        let bounds = &self.config.bounds;
        let damping = 0.5; // Velocity damping on boundary collision

        for i in 0..self.particles.len() {
            if !self.particles.types[i].is_dynamic() {
                continue;
            }

            let base = i * 4;
            let mut pos = [
                self.particles.positions[base],
                self.particles.positions[base + 1],
                self.particles.positions[base + 2],
            ];
            let mut vel = [
                self.particles.velocities[base],
                self.particles.velocities[base + 1],
                self.particles.velocities[base + 2],
            ];

            // X bounds
            if pos[0] < bounds.x_min {
                pos[0] = bounds.x_min;
                vel[0] = -vel[0] * damping;
            } else if pos[0] > bounds.x_max {
                pos[0] = bounds.x_max;
                vel[0] = -vel[0] * damping;
            }

            // Y bounds
            if pos[1] < bounds.y_min {
                pos[1] = bounds.y_min;
                vel[1] = -vel[1] * damping;
            } else if pos[1] > bounds.y_max {
                pos[1] = bounds.y_max;
                vel[1] = -vel[1] * damping;
            }

            // Z bounds
            if pos[2] < bounds.z_min {
                pos[2] = bounds.z_min;
                vel[2] = -vel[2] * damping;
            } else if pos[2] > bounds.z_max {
                pos[2] = bounds.z_max;
                vel[2] = -vel[2] * damping;
            }

            self.particles.set_position(i, pos);
            self.particles.set_velocity(i, vel);
        }
    }

    // ========== Getters ==========

    /// Get current simulation time
    pub fn time(&self) -> f64 {
        self.time
    }

    /// Get current step count
    pub fn step_count(&self) -> u64 {
        self.step
    }

    /// Get number of particles
    pub fn num_particles(&self) -> usize {
        self.particles.len()
    }

    /// Get particle positions as flat array
    pub fn positions(&self) -> &[f32] {
        self.particles.positions_raw()
    }

    /// Get particle velocities as flat array
    pub fn velocities(&self) -> &[f32] {
        self.particles.velocities_raw()
    }

    /// Get particle densities
    pub fn densities(&self) -> &[f32] {
        &self.particles.densities
    }

    /// Get single particle position
    pub fn get_position(&self, idx: usize) -> [f32; 3] {
        self.particles.get_position(idx)
    }

    /// Get single particle velocity
    pub fn get_velocity(&self, idx: usize) -> [f32; 3] {
        self.particles.get_velocity(idx)
    }

    /// Get configuration
    pub fn config(&self) -> &SphConfig {
        &self.config
    }

    /// Get particle buffer reference
    pub fn particles(&self) -> &ParticleBuffer {
        &self.particles
    }

    /// Get mutable particle buffer reference
    pub fn particles_mut(&mut self) -> &mut ParticleBuffer {
        &mut self.particles
    }

    /// Get elastic network reference
    pub fn elastic(&self) -> &ElasticNetwork {
        &self.elastic
    }

    /// Get mutable elastic network reference
    pub fn elastic_mut(&mut self) -> &mut ElasticNetwork {
        &mut self.elastic
    }

    /// Get membrane reference
    pub fn membrane(&self) -> &Membrane {
        &self.membrane
    }

    /// Get mutable membrane reference
    pub fn membrane_mut(&mut self) -> &mut Membrane {
        &mut self.membrane
    }

    /// Reset simulation to initial state
    pub fn reset(&mut self) {
        self.particles.clear();
        self.elastic.clear();
        self.membrane.clear();
        self.muscle_activation = MuscleActivation::new();
        self.time = 0.0;
        self.step = 0;
    }

    /// Calculate total kinetic energy
    pub fn kinetic_energy(&self) -> f64 {
        let mass = self.config.physics.mass as f64;
        let mut energy = 0.0;

        for i in 0..self.particles.len() {
            let vel = self.particles.get_velocity(i);
            let speed_sq = (vel[0] * vel[0] + vel[1] * vel[1] + vel[2] * vel[2]) as f64;
            energy += 0.5 * mass * speed_sq;
        }

        energy
    }

    /// Calculate center of mass
    pub fn center_of_mass(&self) -> [f64; 3] {
        if self.particles.is_empty() {
            return [0.0, 0.0, 0.0];
        }

        let mut com = [0.0_f64; 3];
        let n = self.particles.len();

        for i in 0..n {
            let pos = self.particles.get_position(i);
            com[0] += pos[0] as f64;
            com[1] += pos[1] as f64;
            com[2] += pos[2] as f64;
        }

        let inv_n = 1.0 / n as f64;
        [com[0] * inv_n, com[1] * inv_n, com[2] * inv_n]
    }
}

impl Default for SphWorld {
    fn default() -> Self {
        Self::new(SphConfig::default())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_world_creation() {
        let world = SphWorld::new(SphConfig::default());
        assert_eq!(world.num_particles(), 0);
        assert_eq!(world.time(), 0.0);
    }

    #[test]
    fn test_add_particles() {
        let mut world = SphWorld::fluid();

        let idx = world.add_particle([0.0, 0.0, 0.0], [0.0, 0.0, 0.0], ParticleType::Liquid).unwrap();
        assert_eq!(idx, 0);
        assert_eq!(world.num_particles(), 1);
    }

    #[test]
    fn test_particle_block() {
        let mut world = SphWorld::fluid();

        let indices = world.add_particle_block(
            [0.0, 0.0, 0.0],
            [1.0, 1.0, 1.0],
            0.5,
            ParticleType::Liquid,
        ).unwrap();

        // 3x3x3 = 27 particles
        assert_eq!(indices.len(), 27);
        assert_eq!(world.num_particles(), 27);
    }

    #[test]
    fn test_simulation_step() {
        let mut world = SphWorld::fluid();

        world.add_particle_block(
            [0.0, 10.0, 0.0],
            [5.0, 15.0, 5.0],
            1.0,
            ParticleType::Liquid,
        ).unwrap();

        let initial_com = world.center_of_mass();

        // Run a few steps
        for _ in 0..10 {
            world.step();
        }

        let final_com = world.center_of_mass();

        // Particles should fall due to gravity
        assert!(final_com[1] < initial_com[1], "Particles should fall");
    }

    #[test]
    fn test_elastic_connections() {
        let mut world = SphWorld::fluid();

        let p1 = world.add_particle([0.0, 0.0, 0.0], [0.0, 0.0, 0.0], ParticleType::Elastic).unwrap();
        let p2 = world.add_particle([1.0, 0.0, 0.0], [0.0, 0.0, 0.0], ParticleType::Elastic).unwrap();

        world.connect(p1, p2, 100.0);

        assert_eq!(world.elastic().num_connections(), 1);
    }
}
