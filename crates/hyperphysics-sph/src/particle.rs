//! Particle types and buffers for SPH simulation
//!
//! Particles are the fundamental simulation units in SPH. This module defines
//! particle types, state, and efficient buffer layouts for both CPU and GPU.

use serde::{Deserialize, Serialize};
use smallvec::SmallVec;

/// Particle type classification
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
#[repr(u8)]
pub enum ParticleType {
    /// Liquid particle (water, blood, etc.)
    Liquid = 1,
    /// Elastic particle (deformable solid)
    Elastic = 2,
    /// Boundary particle (static obstacle)
    Boundary = 3,
}

impl ParticleType {
    /// Check if this particle type can move
    #[inline]
    pub fn is_dynamic(&self) -> bool {
        !matches!(self, Self::Boundary)
    }

    /// Check if this particle type participates in pressure calculation
    #[inline]
    pub fn has_pressure(&self) -> bool {
        matches!(self, Self::Liquid | Self::Elastic)
    }
}

/// A single SPH particle
#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
#[repr(C)]
pub struct Particle {
    /// Position (x, y, z) + cell_id (w) - aligned for GPU
    pub position: [f32; 4],
    /// Velocity (x, y, z) + particle_type (w)
    pub velocity: [f32; 4],
    /// Density and inverse density
    pub density: f32,
    pub density_inv: f32,
    /// Pressure
    pub pressure: f32,
    /// Particle index (for tracking)
    pub index: u32,
}

impl Particle {
    /// Create a new particle
    pub fn new(position: [f32; 3], velocity: [f32; 3], particle_type: ParticleType) -> Self {
        Self {
            position: [position[0], position[1], position[2], 0.0],
            velocity: [velocity[0], velocity[1], velocity[2], particle_type as u8 as f32],
            density: 0.0,
            density_inv: 0.0,
            pressure: 0.0,
            index: 0,
        }
    }

    /// Get position as 3D vector
    #[inline]
    pub fn pos(&self) -> [f32; 3] {
        [self.position[0], self.position[1], self.position[2]]
    }

    /// Get velocity as 3D vector
    #[inline]
    pub fn vel(&self) -> [f32; 3] {
        [self.velocity[0], self.velocity[1], self.velocity[2]]
    }

    /// Get particle type
    #[inline]
    pub fn particle_type(&self) -> ParticleType {
        match self.velocity[3] as u8 {
            1 => ParticleType::Liquid,
            2 => ParticleType::Elastic,
            3 => ParticleType::Boundary,
            _ => ParticleType::Liquid,
        }
    }

    /// Set position
    #[inline]
    pub fn set_pos(&mut self, pos: [f32; 3]) {
        self.position[0] = pos[0];
        self.position[1] = pos[1];
        self.position[2] = pos[2];
    }

    /// Set velocity
    #[inline]
    pub fn set_vel(&mut self, vel: [f32; 3]) {
        self.velocity[0] = vel[0];
        self.velocity[1] = vel[1];
        self.velocity[2] = vel[2];
    }

    /// Get cell ID from position.w
    #[inline]
    pub fn cell_id(&self) -> i32 {
        self.position[3] as i32
    }

    /// Set cell ID
    #[inline]
    pub fn set_cell_id(&mut self, cell_id: i32) {
        self.position[3] = cell_id as f32;
    }
}

impl Default for Particle {
    fn default() -> Self {
        Self::new([0.0, 0.0, 0.0], [0.0, 0.0, 0.0], ParticleType::Liquid)
    }
}

/// Structure-of-Arrays particle buffer for cache-efficient access
///
/// This layout is optimal for SIMD operations and GPU memory coalescing.
#[derive(Debug, Clone, Default)]
pub struct ParticleBuffer {
    /// Positions: x, y, z, cell_id (interleaved float4)
    pub positions: Vec<f32>,
    /// Velocities: vx, vy, vz, particle_type (interleaved float4)
    pub velocities: Vec<f32>,
    /// Sorted positions (for neighbor search)
    pub sorted_positions: Vec<f32>,
    /// Sorted velocities
    pub sorted_velocities: Vec<f32>,
    /// Density values
    pub densities: Vec<f32>,
    /// Pressure values
    pub pressures: Vec<f32>,
    /// Acceleration values (x, y, z, _)
    pub accelerations: Vec<f32>,
    /// Neighbor map: particle_id, distance pairs
    pub neighbor_map: Vec<f32>,
    /// Particle index for sorting: cell_id, original_id
    pub particle_index: Vec<u32>,
    /// Back-reference to original particle order
    pub particle_index_back: Vec<u32>,
    /// Grid cell start indices
    pub grid_cell_index: Vec<u32>,
    /// Grid cell end indices
    pub grid_cell_index_end: Vec<u32>,
    /// Particle types
    pub types: Vec<ParticleType>,
    /// Maximum neighbors per particle
    pub max_neighbors: usize,
}

impl ParticleBuffer {
    /// Create a new particle buffer with given capacity
    pub fn with_capacity(capacity: usize, max_neighbors: usize) -> Self {
        Self {
            positions: Vec::with_capacity(capacity * 4),
            velocities: Vec::with_capacity(capacity * 4),
            sorted_positions: Vec::with_capacity(capacity * 4),
            sorted_velocities: Vec::with_capacity(capacity * 4),
            densities: Vec::with_capacity(capacity),
            pressures: Vec::with_capacity(capacity),
            accelerations: Vec::with_capacity(capacity * 4),
            neighbor_map: Vec::with_capacity(capacity * max_neighbors * 2),
            particle_index: Vec::with_capacity(capacity * 2),
            particle_index_back: Vec::with_capacity(capacity),
            grid_cell_index: Vec::new(), // Sized based on grid
            grid_cell_index_end: Vec::new(),
            types: Vec::with_capacity(capacity),
            max_neighbors,
        }
    }

    /// Get number of particles
    #[inline]
    pub fn len(&self) -> usize {
        self.positions.len() / 4
    }

    /// Check if buffer is empty
    #[inline]
    pub fn is_empty(&self) -> bool {
        self.positions.is_empty()
    }

    /// Clear all particles
    pub fn clear(&mut self) {
        self.positions.clear();
        self.velocities.clear();
        self.sorted_positions.clear();
        self.sorted_velocities.clear();
        self.densities.clear();
        self.pressures.clear();
        self.accelerations.clear();
        self.neighbor_map.clear();
        self.particle_index.clear();
        self.particle_index_back.clear();
        self.types.clear();
    }

    /// Add a particle to the buffer
    pub fn push(&mut self, particle: &Particle) {
        self.positions.extend_from_slice(&particle.position);
        self.velocities.extend_from_slice(&particle.velocity);
        self.densities.push(particle.density);
        self.pressures.push(particle.pressure);
        self.accelerations.extend_from_slice(&[0.0, 0.0, 0.0, 0.0]);
        self.types.push(particle.particle_type());

        // Extend neighbor map for this particle
        for _ in 0..self.max_neighbors * 2 {
            self.neighbor_map.push(-1.0);
        }

        // Update indices
        let idx = self.len() - 1;
        self.particle_index.extend_from_slice(&[0, idx as u32]);
        self.particle_index_back.push(idx as u32);
    }

    /// Get position of particle at index
    #[inline]
    pub fn get_position(&self, idx: usize) -> [f32; 3] {
        let base = idx * 4;
        [
            self.positions[base],
            self.positions[base + 1],
            self.positions[base + 2],
        ]
    }

    /// Set position of particle at index
    #[inline]
    pub fn set_position(&mut self, idx: usize, pos: [f32; 3]) {
        let base = idx * 4;
        self.positions[base] = pos[0];
        self.positions[base + 1] = pos[1];
        self.positions[base + 2] = pos[2];
    }

    /// Get velocity of particle at index
    #[inline]
    pub fn get_velocity(&self, idx: usize) -> [f32; 3] {
        let base = idx * 4;
        [
            self.velocities[base],
            self.velocities[base + 1],
            self.velocities[base + 2],
        ]
    }

    /// Set velocity of particle at index
    #[inline]
    pub fn set_velocity(&mut self, idx: usize, vel: [f32; 3]) {
        let base = idx * 4;
        self.velocities[base] = vel[0];
        self.velocities[base + 1] = vel[1];
        self.velocities[base + 2] = vel[2];
    }

    /// Get acceleration of particle at index
    #[inline]
    pub fn get_acceleration(&self, idx: usize) -> [f32; 3] {
        let base = idx * 4;
        [
            self.accelerations[base],
            self.accelerations[base + 1],
            self.accelerations[base + 2],
        ]
    }

    /// Set acceleration of particle at index
    #[inline]
    pub fn set_acceleration(&mut self, idx: usize, acc: [f32; 3]) {
        let base = idx * 4;
        self.accelerations[base] = acc[0];
        self.accelerations[base + 1] = acc[1];
        self.accelerations[base + 2] = acc[2];
    }

    /// Add to acceleration of particle at index
    #[inline]
    pub fn add_acceleration(&mut self, idx: usize, acc: [f32; 3]) {
        let base = idx * 4;
        self.accelerations[base] += acc[0];
        self.accelerations[base + 1] += acc[1];
        self.accelerations[base + 2] += acc[2];
    }

    /// Get neighbors for particle at index
    pub fn get_neighbors(&self, idx: usize) -> SmallVec<[(u32, f32); 32]> {
        let mut neighbors = SmallVec::new();
        let base = idx * self.max_neighbors * 2;

        for i in 0..self.max_neighbors {
            let neighbor_idx = self.neighbor_map[base + i * 2] as i32;
            if neighbor_idx >= 0 {
                let distance = self.neighbor_map[base + i * 2 + 1];
                neighbors.push((neighbor_idx as u32, distance));
            }
        }

        neighbors
    }

    /// Resize grid cell indices
    pub fn resize_grid(&mut self, grid_size: usize) {
        self.grid_cell_index.resize(grid_size, u32::MAX);
        self.grid_cell_index_end.resize(grid_size, 0);
    }

    /// Get raw position slice for GPU upload
    pub fn positions_raw(&self) -> &[f32] {
        &self.positions
    }

    /// Get raw velocity slice for GPU upload
    pub fn velocities_raw(&self) -> &[f32] {
        &self.velocities
    }

    /// Initialize sorted buffers (copy from unsorted)
    pub fn init_sorted_buffers(&mut self) {
        self.sorted_positions = self.positions.clone();
        self.sorted_velocities = self.velocities.clone();
    }

    /// Convert to a vector of Particle structs
    pub fn to_particles(&self) -> Vec<Particle> {
        let n = self.len();
        let mut particles = Vec::with_capacity(n);

        for i in 0..n {
            let base = i * 4;
            particles.push(Particle {
                position: [
                    self.positions[base],
                    self.positions[base + 1],
                    self.positions[base + 2],
                    self.positions[base + 3],
                ],
                velocity: [
                    self.velocities[base],
                    self.velocities[base + 1],
                    self.velocities[base + 2],
                    self.velocities[base + 3],
                ],
                density: self.densities[i],
                density_inv: if self.densities[i] > 0.0 { 1.0 / self.densities[i] } else { 0.0 },
                pressure: self.pressures[i],
                index: i as u32,
            });
        }

        particles
    }
}

/// Neighbor data for a single particle
#[derive(Debug, Clone, Default)]
pub struct NeighborList {
    /// Neighbor particle indices
    pub indices: SmallVec<[u32; 32]>,
    /// Distances to neighbors
    pub distances: SmallVec<[f32; 32]>,
}

impl NeighborList {
    /// Create empty neighbor list
    pub fn new() -> Self {
        Self::default()
    }

    /// Clear the neighbor list
    pub fn clear(&mut self) {
        self.indices.clear();
        self.distances.clear();
    }

    /// Add a neighbor
    pub fn push(&mut self, index: u32, distance: f32) {
        self.indices.push(index);
        self.distances.push(distance);
    }

    /// Get number of neighbors
    pub fn len(&self) -> usize {
        self.indices.len()
    }

    /// Check if empty
    pub fn is_empty(&self) -> bool {
        self.indices.is_empty()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_particle_creation() {
        let p = Particle::new([1.0, 2.0, 3.0], [0.1, 0.2, 0.3], ParticleType::Liquid);
        assert_eq!(p.pos(), [1.0, 2.0, 3.0]);
        assert_eq!(p.vel(), [0.1, 0.2, 0.3]);
        assert_eq!(p.particle_type(), ParticleType::Liquid);
    }

    #[test]
    fn test_particle_buffer() {
        let mut buffer = ParticleBuffer::with_capacity(100, 32);

        let p1 = Particle::new([0.0, 0.0, 0.0], [0.0, 0.0, 0.0], ParticleType::Liquid);
        let p2 = Particle::new([1.0, 1.0, 1.0], [0.0, 0.0, 0.0], ParticleType::Elastic);

        buffer.push(&p1);
        buffer.push(&p2);

        assert_eq!(buffer.len(), 2);
        assert_eq!(buffer.get_position(0), [0.0, 0.0, 0.0]);
        assert_eq!(buffer.get_position(1), [1.0, 1.0, 1.0]);
    }
}
