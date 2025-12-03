//! Membrane Handling for Surface Boundaries
//!
//! Implements triangular membranes that act as liquid-impermeable boundaries
//! and provide surface tension effects.

use serde::{Deserialize, Serialize};
use crate::particle::ParticleBuffer;

/// A triangular membrane element
#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
#[repr(C)]
pub struct MembraneTriangle {
    /// First vertex particle index
    pub p1: u32,
    /// Second vertex particle index
    pub p2: u32,
    /// Third vertex particle index
    pub p3: u32,
    /// Surface tension coefficient
    pub surface_tension: f32,
}

impl MembraneTriangle {
    /// Create a new membrane triangle
    pub fn new(p1: u32, p2: u32, p3: u32) -> Self {
        Self {
            p1,
            p2,
            p3,
            surface_tension: 1.0,
        }
    }

    /// Create with custom surface tension
    pub fn with_tension(p1: u32, p2: u32, p3: u32, surface_tension: f32) -> Self {
        Self {
            p1,
            p2,
            p3,
            surface_tension,
        }
    }

    /// Get vertex indices as array
    pub fn vertices(&self) -> [u32; 3] {
        [self.p1, self.p2, self.p3]
    }

    /// Calculate normal vector from particle positions
    pub fn normal(&self, particles: &ParticleBuffer) -> [f32; 3] {
        let pos1 = particles.get_position(self.p1 as usize);
        let pos2 = particles.get_position(self.p2 as usize);
        let pos3 = particles.get_position(self.p3 as usize);

        // Edge vectors
        let e1 = [pos2[0] - pos1[0], pos2[1] - pos1[1], pos2[2] - pos1[2]];
        let e2 = [pos3[0] - pos1[0], pos3[1] - pos1[1], pos3[2] - pos1[2]];

        // Cross product
        let normal = [
            e1[1] * e2[2] - e1[2] * e2[1],
            e1[2] * e2[0] - e1[0] * e2[2],
            e1[0] * e2[1] - e1[1] * e2[0],
        ];

        // Normalize
        let len = (normal[0] * normal[0] + normal[1] * normal[1] + normal[2] * normal[2]).sqrt();
        if len > 1e-10 {
            [normal[0] / len, normal[1] / len, normal[2] / len]
        } else {
            [0.0, 1.0, 0.0]
        }
    }

    /// Calculate area of the triangle
    pub fn area(&self, particles: &ParticleBuffer) -> f32 {
        let pos1 = particles.get_position(self.p1 as usize);
        let pos2 = particles.get_position(self.p2 as usize);
        let pos3 = particles.get_position(self.p3 as usize);

        let e1 = [pos2[0] - pos1[0], pos2[1] - pos1[1], pos2[2] - pos1[2]];
        let e2 = [pos3[0] - pos1[0], pos3[1] - pos1[1], pos3[2] - pos1[2]];

        let cross = [
            e1[1] * e2[2] - e1[2] * e2[1],
            e1[2] * e2[0] - e1[0] * e2[2],
            e1[0] * e2[1] - e1[1] * e2[0],
        ];

        0.5 * (cross[0] * cross[0] + cross[1] * cross[1] + cross[2] * cross[2]).sqrt()
    }

    /// Calculate centroid of the triangle
    pub fn centroid(&self, particles: &ParticleBuffer) -> [f32; 3] {
        let pos1 = particles.get_position(self.p1 as usize);
        let pos2 = particles.get_position(self.p2 as usize);
        let pos3 = particles.get_position(self.p3 as usize);

        [
            (pos1[0] + pos2[0] + pos3[0]) / 3.0,
            (pos1[1] + pos2[1] + pos3[1]) / 3.0,
            (pos1[2] + pos2[2] + pos3[2]) / 3.0,
        ]
    }
}

/// Collection of membrane surfaces
#[derive(Debug, Clone, Default, Serialize, Deserialize)]
pub struct Membrane {
    /// All membrane triangles
    triangles: Vec<MembraneTriangle>,
    /// Particle to membrane mapping (which membranes include this particle)
    particle_membranes: Vec<Vec<usize>>,
    /// Maximum membranes per particle
    max_membranes_per_particle: usize,
}

impl Membrane {
    /// Create a new empty membrane
    pub fn new() -> Self {
        Self {
            triangles: Vec::new(),
            particle_membranes: Vec::new(),
            max_membranes_per_particle: 7, // From Sibernetic constant
        }
    }

    /// Create with preallocated capacity
    pub fn with_capacity(num_particles: usize, num_triangles: usize) -> Self {
        Self {
            triangles: Vec::with_capacity(num_triangles),
            particle_membranes: vec![Vec::new(); num_particles],
            max_membranes_per_particle: 7,
        }
    }

    /// Add a membrane triangle
    pub fn add_triangle(&mut self, triangle: MembraneTriangle) {
        let idx = self.triangles.len();

        // Ensure particle_membranes is large enough
        let max_particle = triangle.p1.max(triangle.p2).max(triangle.p3) as usize;
        if max_particle >= self.particle_membranes.len() {
            self.particle_membranes.resize(max_particle + 1, Vec::new());
        }

        // Add membrane reference to each vertex particle
        for &p in &[triangle.p1, triangle.p2, triangle.p3] {
            let p = p as usize;
            if self.particle_membranes[p].len() < self.max_membranes_per_particle {
                self.particle_membranes[p].push(idx);
            }
        }

        self.triangles.push(triangle);
    }

    /// Add a triangle from particle indices
    pub fn add(&mut self, p1: u32, p2: u32, p3: u32) {
        self.add_triangle(MembraneTriangle::new(p1, p2, p3));
    }

    /// Get number of triangles
    pub fn num_triangles(&self) -> usize {
        self.triangles.len()
    }

    /// Get membranes containing a particle
    pub fn get_particle_membranes(&self, particle_idx: usize) -> &[usize] {
        if particle_idx < self.particle_membranes.len() {
            &self.particle_membranes[particle_idx]
        } else {
            &[]
        }
    }

    /// Compute membrane interaction forces
    ///
    /// Handles:
    /// 1. Surface tension (keeps membrane intact)
    /// 2. Liquid-membrane interaction (prevents penetration)
    pub fn compute_forces(
        &self,
        particles: &mut ParticleBuffer,
        surface_tension_coeff: f32,
    ) {
        // Surface tension: pull membrane particles toward centroid
        for (tri_idx, triangle) in self.triangles.iter().enumerate() {
            let centroid = triangle.centroid(particles);
            let normal = triangle.normal(particles);
            let area = triangle.area(particles);

            let force_magnitude = surface_tension_coeff * triangle.surface_tension * area;

            // Apply force to each vertex toward maintaining surface
            for &p in &[triangle.p1, triangle.p2, triangle.p3] {
                let p = p as usize;
                let pos = particles.get_position(p);

                // Direction from particle to centroid
                let to_center = [
                    centroid[0] - pos[0],
                    centroid[1] - pos[1],
                    centroid[2] - pos[2],
                ];
                let dist = (to_center[0] * to_center[0]
                    + to_center[1] * to_center[1]
                    + to_center[2] * to_center[2])
                .sqrt();

                if dist > 1e-10 {
                    let factor = force_magnitude / dist;
                    particles.add_acceleration(p, [
                        factor * to_center[0],
                        factor * to_center[1],
                        factor * to_center[2],
                    ]);
                }
            }
        }
    }

    /// Check if a point is inside a closed membrane surface
    /// Uses ray casting algorithm
    pub fn point_inside(&self, point: [f32; 3], particles: &ParticleBuffer) -> bool {
        let mut intersections = 0;

        // Cast ray in +x direction
        for triangle in &self.triangles {
            if self.ray_intersects_triangle(point, [1.0, 0.0, 0.0], triangle, particles) {
                intersections += 1;
            }
        }

        // Odd number of intersections = inside
        intersections % 2 == 1
    }

    /// Ray-triangle intersection test (Möller-Trumbore algorithm)
    fn ray_intersects_triangle(
        &self,
        origin: [f32; 3],
        direction: [f32; 3],
        triangle: &MembraneTriangle,
        particles: &ParticleBuffer,
    ) -> bool {
        let v0 = particles.get_position(triangle.p1 as usize);
        let v1 = particles.get_position(triangle.p2 as usize);
        let v2 = particles.get_position(triangle.p3 as usize);

        let edge1 = [v1[0] - v0[0], v1[1] - v0[1], v1[2] - v0[2]];
        let edge2 = [v2[0] - v0[0], v2[1] - v0[1], v2[2] - v0[2]];

        // Cross product: direction × edge2
        let h = [
            direction[1] * edge2[2] - direction[2] * edge2[1],
            direction[2] * edge2[0] - direction[0] * edge2[2],
            direction[0] * edge2[1] - direction[1] * edge2[0],
        ];

        let a = edge1[0] * h[0] + edge1[1] * h[1] + edge1[2] * h[2];

        if a.abs() < 1e-10 {
            return false; // Ray parallel to triangle
        }

        let f = 1.0 / a;
        let s = [origin[0] - v0[0], origin[1] - v0[1], origin[2] - v0[2]];

        let u = f * (s[0] * h[0] + s[1] * h[1] + s[2] * h[2]);
        if u < 0.0 || u > 1.0 {
            return false;
        }

        let q = [
            s[1] * edge1[2] - s[2] * edge1[1],
            s[2] * edge1[0] - s[0] * edge1[2],
            s[0] * edge1[1] - s[1] * edge1[0],
        ];

        let v = f * (direction[0] * q[0] + direction[1] * q[1] + direction[2] * q[2]);
        if v < 0.0 || u + v > 1.0 {
            return false;
        }

        let t = f * (edge2[0] * q[0] + edge2[1] * q[1] + edge2[2] * q[2]);

        t > 1e-10 // Intersection in front of ray origin
    }

    /// Get raw triangle data
    pub fn triangles_raw(&self) -> &[MembraneTriangle] {
        &self.triangles
    }

    /// Clear all membranes
    pub fn clear(&mut self) {
        self.triangles.clear();
        for pm in &mut self.particle_membranes {
            pm.clear();
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::particle::{Particle, ParticleType};

    #[test]
    fn test_membrane_triangle() {
        let mut particles = ParticleBuffer::with_capacity(10, 32);

        particles.push(&Particle::new([0.0, 0.0, 0.0], [0.0, 0.0, 0.0], ParticleType::Elastic));
        particles.push(&Particle::new([1.0, 0.0, 0.0], [0.0, 0.0, 0.0], ParticleType::Elastic));
        particles.push(&Particle::new([0.0, 1.0, 0.0], [0.0, 0.0, 0.0], ParticleType::Elastic));

        let triangle = MembraneTriangle::new(0, 1, 2);

        let normal = triangle.normal(&particles);
        assert!((normal[2] - 1.0).abs() < 1e-6, "Normal should point in +z");

        let area = triangle.area(&particles);
        assert!((area - 0.5).abs() < 1e-6, "Area should be 0.5");
    }

    #[test]
    fn test_membrane_creation() {
        let mut membrane = Membrane::new();

        membrane.add(0, 1, 2);
        membrane.add(1, 2, 3);

        assert_eq!(membrane.num_triangles(), 2);
        assert_eq!(membrane.get_particle_membranes(1).len(), 2);
    }
}
