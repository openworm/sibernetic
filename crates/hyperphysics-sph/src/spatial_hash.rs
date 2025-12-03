//! Spatial hashing for efficient neighbor search
//!
//! Implements a uniform grid spatial hash for O(n) neighbor finding.
//! Based on the algorithm from Sibernetic's sphFluid.cl.

use crate::config::PhysicsConstants;
use crate::particle::ParticleBuffer;
use hashbrown::HashMap;

/// Spatial hash grid for efficient neighbor queries
#[derive(Debug, Clone)]
pub struct SpatialHash {
    /// Cell size (should be >= smoothing radius h)
    cell_size: f32,
    /// Inverse cell size for fast division
    cell_size_inv: f32,
    /// Grid dimensions
    grid_dims: [u32; 3],
    /// Grid origin (minimum bounds)
    origin: [f32; 3],
    /// Total number of cells
    num_cells: usize,
}

impl SpatialHash {
    /// Create a new spatial hash grid
    pub fn new(
        cell_size: f32,
        bounds_min: [f32; 3],
        bounds_max: [f32; 3],
    ) -> Self {
        let cell_size_inv = 1.0 / cell_size;

        let grid_dims = [
            ((bounds_max[0] - bounds_min[0]) * cell_size_inv).ceil() as u32 + 1,
            ((bounds_max[1] - bounds_min[1]) * cell_size_inv).ceil() as u32 + 1,
            ((bounds_max[2] - bounds_min[2]) * cell_size_inv).ceil() as u32 + 1,
        ];

        let num_cells = grid_dims[0] as usize * grid_dims[1] as usize * grid_dims[2] as usize;

        Self {
            cell_size,
            cell_size_inv,
            grid_dims,
            origin: bounds_min,
            num_cells,
        }
    }

    /// Create from physics constants and bounds
    pub fn from_config(
        physics: &PhysicsConstants,
        bounds_min: [f32; 3],
        bounds_max: [f32; 3],
    ) -> Self {
        Self::new(physics.hash_grid_cell_size, bounds_min, bounds_max)
    }

    /// Get cell factors (integer grid coordinates) for a position
    #[inline]
    pub fn cell_factors(&self, position: [f32; 3]) -> [i32; 3] {
        [
            ((position[0] - self.origin[0]) * self.cell_size_inv) as i32,
            ((position[1] - self.origin[1]) * self.cell_size_inv) as i32,
            ((position[2] - self.origin[2]) * self.cell_size_inv) as i32,
        ]
    }

    /// Convert cell factors to cell ID
    #[inline]
    pub fn cell_id(&self, factors: [i32; 3]) -> i32 {
        if factors[0] < 0 || factors[1] < 0 || factors[2] < 0
            || factors[0] >= self.grid_dims[0] as i32
            || factors[1] >= self.grid_dims[1] as i32
            || factors[2] >= self.grid_dims[2] as i32
        {
            return -1;
        }

        factors[0]
            + factors[1] * self.grid_dims[0] as i32
            + factors[2] * self.grid_dims[0] as i32 * self.grid_dims[1] as i32
    }

    /// Get cell ID directly from position
    #[inline]
    pub fn hash(&self, position: [f32; 3]) -> i32 {
        let factors = self.cell_factors(position);
        self.cell_id(factors)
    }

    /// Get neighboring cell IDs for a position (3x3x3 = 27 cells)
    pub fn neighbor_cells(&self, position: [f32; 3]) -> impl Iterator<Item = i32> + '_ {
        let factors = self.cell_factors(position);

        // Generate all 27 neighboring cell offsets
        let offsets: [(i32, i32, i32); 27] = [
            (-1, -1, -1), (-1, -1, 0), (-1, -1, 1),
            (-1, 0, -1), (-1, 0, 0), (-1, 0, 1),
            (-1, 1, -1), (-1, 1, 0), (-1, 1, 1),
            (0, -1, -1), (0, -1, 0), (0, -1, 1),
            (0, 0, -1), (0, 0, 0), (0, 0, 1),
            (0, 1, -1), (0, 1, 0), (0, 1, 1),
            (1, -1, -1), (1, -1, 0), (1, -1, 1),
            (1, 0, -1), (1, 0, 0), (1, 0, 1),
            (1, 1, -1), (1, 1, 0), (1, 1, 1),
        ];

        offsets.into_iter().filter_map(move |(dx, dy, dz)| {
            let neighbor = [
                factors[0] + dx,
                factors[1] + dy,
                factors[2] + dz,
            ];
            let id = self.cell_id(neighbor);
            if id >= 0 { Some(id) } else { None }
        })
    }

    /// Get total number of cells
    pub fn num_cells(&self) -> usize {
        self.num_cells
    }

    /// Get grid dimensions
    pub fn grid_dims(&self) -> [u32; 3] {
        self.grid_dims
    }

    /// Hash all particles and update their cell IDs
    pub fn hash_particles(&self, particles: &mut ParticleBuffer) {
        let n = particles.len();

        for i in 0..n {
            let pos = particles.get_position(i);
            let cell_id = self.hash(pos);

            // Store cell ID in position.w
            let base = i * 4;
            particles.positions[base + 3] = cell_id as f32;

            // Update particle index: (cell_id, original_index)
            particles.particle_index[i * 2] = (cell_id & 0xFFFFFF) as u32;
            particles.particle_index[i * 2 + 1] = i as u32;
        }
    }

    /// Sort particles by cell ID (counting sort)
    pub fn sort_particles(&self, particles: &mut ParticleBuffer) {
        let n = particles.len();
        if n == 0 {
            return;
        }

        // Build (cell_id, original_index) pairs and sort
        let mut indices: Vec<(u32, u32)> = (0..n as u32)
            .map(|i| {
                let cell_id = particles.particle_index[i as usize * 2];
                (cell_id, i)
            })
            .collect();

        indices.sort_unstable_by_key(|&(cell_id, _)| cell_id);

        // Rearrange particles into sorted order
        for (sorted_idx, &(_, original_idx)) in indices.iter().enumerate() {
            let orig = original_idx as usize;
            let base_orig = orig * 4;
            let base_sorted = sorted_idx * 4;

            // Copy position
            particles.sorted_positions[base_sorted] = particles.positions[base_orig];
            particles.sorted_positions[base_sorted + 1] = particles.positions[base_orig + 1];
            particles.sorted_positions[base_sorted + 2] = particles.positions[base_orig + 2];
            particles.sorted_positions[base_sorted + 3] = particles.positions[base_orig + 3];

            // Copy velocity
            particles.sorted_velocities[base_sorted] = particles.velocities[base_orig];
            particles.sorted_velocities[base_sorted + 1] = particles.velocities[base_orig + 1];
            particles.sorted_velocities[base_sorted + 2] = particles.velocities[base_orig + 2];
            particles.sorted_velocities[base_sorted + 3] = particles.velocities[base_orig + 3];

            // Update back reference
            particles.particle_index_back[orig] = sorted_idx as u32;
        }

        // Update particle_index with sorted order
        for (sorted_idx, &(cell_id, original_idx)) in indices.iter().enumerate() {
            particles.particle_index[sorted_idx * 2] = cell_id;
            particles.particle_index[sorted_idx * 2 + 1] = original_idx;
        }
    }

    /// Build grid cell index (start/end indices for each cell)
    pub fn build_cell_index(&self, particles: &mut ParticleBuffer) {
        let n = particles.len();

        // Reset grid indices
        particles.resize_grid(self.num_cells);
        for i in 0..self.num_cells {
            particles.grid_cell_index[i] = u32::MAX;
            particles.grid_cell_index_end[i] = 0;
        }

        if n == 0 {
            return;
        }

        // Find start and end indices for each cell
        let mut current_cell = particles.particle_index[0];
        particles.grid_cell_index[current_cell as usize] = 0;

        for i in 1..n {
            let cell_id = particles.particle_index[i * 2];
            if cell_id != current_cell {
                particles.grid_cell_index_end[current_cell as usize] = i as u32;
                if (cell_id as usize) < self.num_cells {
                    particles.grid_cell_index[cell_id as usize] = i as u32;
                }
                current_cell = cell_id;
            }
        }
        particles.grid_cell_index_end[current_cell as usize] = n as u32;
    }

    /// Find all neighbors within smoothing radius h
    pub fn find_neighbors(
        &self,
        particles: &ParticleBuffer,
        particle_idx: usize,
        h_sq: f32,
        max_neighbors: usize,
    ) -> Vec<(usize, f32)> {
        let mut neighbors = Vec::with_capacity(max_neighbors);
        let pos = particles.get_position(particle_idx);

        // Check all neighboring cells
        for cell_id in self.neighbor_cells(pos) {
            let cell_idx = cell_id as usize;
            if cell_idx >= self.num_cells {
                continue;
            }

            let start = particles.grid_cell_index[cell_idx];
            let end = particles.grid_cell_index_end[cell_idx];

            if start == u32::MAX {
                continue;
            }

            // Check all particles in this cell
            for j in start..end {
                let j = j as usize;
                if j == particle_idx {
                    continue;
                }

                let other_pos = [
                    particles.sorted_positions[j * 4],
                    particles.sorted_positions[j * 4 + 1],
                    particles.sorted_positions[j * 4 + 2],
                ];

                let dx = pos[0] - other_pos[0];
                let dy = pos[1] - other_pos[1];
                let dz = pos[2] - other_pos[2];
                let dist_sq = dx * dx + dy * dy + dz * dz;

                if dist_sq < h_sq && neighbors.len() < max_neighbors {
                    neighbors.push((j, dist_sq.sqrt()));
                }
            }
        }

        neighbors
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_cell_id() {
        let hash = SpatialHash::new(1.0, [0.0, 0.0, 0.0], [10.0, 10.0, 10.0]);

        let factors = hash.cell_factors([0.5, 0.5, 0.5]);
        assert_eq!(factors, [0, 0, 0]);

        let cell_id = hash.cell_id(factors);
        assert_eq!(cell_id, 0);

        let factors2 = hash.cell_factors([1.5, 0.5, 0.5]);
        assert_eq!(factors2, [1, 0, 0]);
    }

    #[test]
    fn test_neighbor_cells() {
        let hash = SpatialHash::new(1.0, [0.0, 0.0, 0.0], [10.0, 10.0, 10.0]);

        let neighbors: Vec<_> = hash.neighbor_cells([5.0, 5.0, 5.0]).collect();
        assert_eq!(neighbors.len(), 27); // 3x3x3 grid
    }
}
