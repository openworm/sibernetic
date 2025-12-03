//! Elastic Connections for Deformable Bodies
//!
//! Implements spring-based connections between particles for simulating
//! elastic and deformable materials like the C. elegans body.

use serde::{Deserialize, Serialize};
use crate::particle::ParticleBuffer;

/// A single elastic connection between two particles
#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
#[repr(C)]
pub struct ElasticConnection {
    /// First particle index
    pub p1: u32,
    /// Second particle index
    pub p2: u32,
    /// Rest length of the connection
    pub rest_length: f32,
    /// Spring stiffness coefficient
    pub stiffness: f32,
    /// Damping coefficient
    pub damping: f32,
    /// Muscle index (-1 if not a muscle)
    pub muscle_id: i32,
}

impl ElasticConnection {
    /// Create a new elastic connection
    pub fn new(p1: u32, p2: u32, rest_length: f32, stiffness: f32) -> Self {
        Self {
            p1,
            p2,
            rest_length,
            stiffness,
            damping: 0.1,
            muscle_id: -1,
        }
    }

    /// Create a muscle connection
    pub fn muscle(p1: u32, p2: u32, rest_length: f32, stiffness: f32, muscle_id: i32) -> Self {
        Self {
            p1,
            p2,
            rest_length,
            stiffness,
            damping: 0.1,
            muscle_id,
        }
    }

    /// Check if this connection is a muscle
    pub fn is_muscle(&self) -> bool {
        self.muscle_id >= 0
    }
}

/// Network of elastic connections
#[derive(Debug, Clone, Default)]
pub struct ElasticNetwork {
    /// All connections
    connections: Vec<ElasticConnection>,
    /// Connection indices by particle (for fast lookup)
    particle_connections: Vec<Vec<usize>>,
    /// Muscle activation values (0.0 to 1.0)
    muscle_activations: Vec<f32>,
    /// Number of muscles
    num_muscles: usize,
}

impl ElasticNetwork {
    /// Create a new empty elastic network
    pub fn new() -> Self {
        Self::default()
    }

    /// Create with preallocated capacity
    pub fn with_capacity(num_particles: usize, num_connections: usize, num_muscles: usize) -> Self {
        Self {
            connections: Vec::with_capacity(num_connections),
            particle_connections: vec![Vec::new(); num_particles],
            muscle_activations: vec![0.0; num_muscles],
            num_muscles,
        }
    }

    /// Add an elastic connection
    pub fn add_connection(&mut self, connection: ElasticConnection) {
        let idx = self.connections.len();

        // Ensure particle_connections is large enough
        let max_particle = connection.p1.max(connection.p2) as usize;
        if max_particle >= self.particle_connections.len() {
            self.particle_connections.resize(max_particle + 1, Vec::new());
        }

        self.particle_connections[connection.p1 as usize].push(idx);
        self.particle_connections[connection.p2 as usize].push(idx);
        self.connections.push(connection);
    }

    /// Add a connection between two particles
    pub fn connect(&mut self, p1: u32, p2: u32, rest_length: f32, stiffness: f32) {
        self.add_connection(ElasticConnection::new(p1, p2, rest_length, stiffness));
    }

    /// Add a muscle connection
    pub fn connect_muscle(
        &mut self,
        p1: u32,
        p2: u32,
        rest_length: f32,
        stiffness: f32,
        muscle_id: i32,
    ) {
        self.add_connection(ElasticConnection::muscle(p1, p2, rest_length, stiffness, muscle_id));

        // Ensure muscle activation array is large enough
        let muscle_idx = muscle_id as usize;
        if muscle_idx >= self.muscle_activations.len() {
            self.muscle_activations.resize(muscle_idx + 1, 0.0);
            self.num_muscles = self.muscle_activations.len();
        }
    }

    /// Set activation for a specific muscle
    pub fn set_muscle_activation(&mut self, muscle_id: usize, activation: f32) {
        if muscle_id < self.muscle_activations.len() {
            self.muscle_activations[muscle_id] = activation.clamp(0.0, 1.0);
        }
    }

    /// Set all muscle activations at once
    pub fn set_all_activations(&mut self, activations: &[f32]) {
        let n = activations.len().min(self.muscle_activations.len());
        self.muscle_activations[..n].copy_from_slice(&activations[..n]);
    }

    /// Get number of connections
    pub fn num_connections(&self) -> usize {
        self.connections.len()
    }

    /// Get number of muscles
    pub fn num_muscles(&self) -> usize {
        self.num_muscles
    }

    /// Get connections for a particle
    pub fn get_particle_connections(&self, particle_idx: usize) -> &[usize] {
        if particle_idx < self.particle_connections.len() {
            &self.particle_connections[particle_idx]
        } else {
            &[]
        }
    }

    /// Compute elastic forces for all connections and add to particle accelerations
    pub fn compute_forces(
        &self,
        particles: &mut ParticleBuffer,
        elasticity_coefficient: f32,
        max_muscle_force: f32,
    ) {
        for connection in &self.connections {
            let p1 = connection.p1 as usize;
            let p2 = connection.p2 as usize;

            // Get positions
            let pos1 = particles.get_position(p1);
            let pos2 = particles.get_position(p2);

            // Calculate distance
            let dx = pos2[0] - pos1[0];
            let dy = pos2[1] - pos1[1];
            let dz = pos2[2] - pos1[2];
            let dist = (dx * dx + dy * dy + dz * dz).sqrt();

            if dist < 1e-10 {
                continue;
            }

            // Calculate rest length (may be modified by muscle activation)
            let mut effective_rest_length = connection.rest_length;

            // Muscle contraction
            if connection.is_muscle() {
                let muscle_idx = connection.muscle_id as usize;
                if muscle_idx < self.muscle_activations.len() {
                    let activation = self.muscle_activations[muscle_idx];
                    // Muscles contract when activated (reduce rest length)
                    effective_rest_length *= 1.0 - 0.3 * activation; // 30% max contraction
                }
            }

            // Spring force: F = k * (length - rest_length)
            let stretch = dist - effective_rest_length;
            let force_magnitude = elasticity_coefficient * connection.stiffness * stretch;

            // Add muscle force on top
            let muscle_force = if connection.is_muscle() {
                let muscle_idx = connection.muscle_id as usize;
                if muscle_idx < self.muscle_activations.len() {
                    max_muscle_force * self.muscle_activations[muscle_idx]
                } else {
                    0.0
                }
            } else {
                0.0
            };

            let total_force = force_magnitude + muscle_force;

            // Normalize direction
            let inv_dist = 1.0 / dist;
            let dir = [dx * inv_dist, dy * inv_dist, dz * inv_dist];

            // Calculate damping force
            let vel1 = particles.get_velocity(p1);
            let vel2 = particles.get_velocity(p2);
            let rel_vel = [
                vel2[0] - vel1[0],
                vel2[1] - vel1[1],
                vel2[2] - vel1[2],
            ];
            let vel_along_spring = rel_vel[0] * dir[0] + rel_vel[1] * dir[1] + rel_vel[2] * dir[2];
            let damping_force = connection.damping * vel_along_spring;

            // Apply forces (equal and opposite)
            let f = [
                (total_force + damping_force) * dir[0],
                (total_force + damping_force) * dir[1],
                (total_force + damping_force) * dir[2],
            ];

            particles.add_acceleration(p1, f);
            particles.add_acceleration(p2, [-f[0], -f[1], -f[2]]);
        }
    }

    /// Get raw connection data for GPU upload
    pub fn connections_raw(&self) -> &[ElasticConnection] {
        &self.connections
    }

    /// Clear all connections
    pub fn clear(&mut self) {
        self.connections.clear();
        for pc in &mut self.particle_connections {
            pc.clear();
        }
        self.muscle_activations.fill(0.0);
    }

    /// Serialize to bytes
    pub fn to_bytes(&self) -> Vec<u8> {
        bincode::serialize(self).unwrap_or_default()
    }

    /// Deserialize from bytes
    pub fn from_bytes(bytes: &[u8]) -> Option<Self> {
        bincode::deserialize(bytes).ok()
    }
}

impl Serialize for ElasticNetwork {
    fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
    where
        S: serde::Serializer,
    {
        use serde::ser::SerializeStruct;
        let mut state = serializer.serialize_struct("ElasticNetwork", 3)?;
        state.serialize_field("connections", &self.connections)?;
        state.serialize_field("muscle_activations", &self.muscle_activations)?;
        state.serialize_field("num_muscles", &self.num_muscles)?;
        state.end()
    }
}

impl<'de> Deserialize<'de> for ElasticNetwork {
    fn deserialize<D>(deserializer: D) -> Result<Self, D::Error>
    where
        D: serde::Deserializer<'de>,
    {
        #[derive(Deserialize)]
        struct ElasticNetworkData {
            connections: Vec<ElasticConnection>,
            muscle_activations: Vec<f32>,
            num_muscles: usize,
        }

        let data = ElasticNetworkData::deserialize(deserializer)?;

        // Rebuild particle_connections
        let max_particle = data.connections.iter()
            .map(|c| c.p1.max(c.p2) as usize)
            .max()
            .unwrap_or(0);

        let mut particle_connections = vec![Vec::new(); max_particle + 1];
        for (idx, conn) in data.connections.iter().enumerate() {
            particle_connections[conn.p1 as usize].push(idx);
            particle_connections[conn.p2 as usize].push(idx);
        }

        Ok(Self {
            connections: data.connections,
            particle_connections,
            muscle_activations: data.muscle_activations,
            num_muscles: data.num_muscles,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::particle::{Particle, ParticleType};

    #[test]
    fn test_elastic_connection() {
        let conn = ElasticConnection::new(0, 1, 1.0, 100.0);
        assert!(!conn.is_muscle());

        let muscle = ElasticConnection::muscle(0, 1, 1.0, 100.0, 5);
        assert!(muscle.is_muscle());
    }

    #[test]
    fn test_elastic_network() {
        let mut network = ElasticNetwork::new();

        network.connect(0, 1, 1.0, 100.0);
        network.connect(1, 2, 1.0, 100.0);
        network.connect_muscle(2, 3, 1.0, 100.0, 0);

        assert_eq!(network.num_connections(), 3);
        assert!(network.num_muscles() >= 1);

        network.set_muscle_activation(0, 0.5);
        assert_eq!(network.muscle_activations[0], 0.5);
    }

    #[test]
    fn test_force_computation() {
        let mut particles = ParticleBuffer::with_capacity(10, 32);

        // Two particles 2 units apart
        particles.push(&Particle::new([0.0, 0.0, 0.0], [0.0, 0.0, 0.0], ParticleType::Elastic));
        particles.push(&Particle::new([2.0, 0.0, 0.0], [0.0, 0.0, 0.0], ParticleType::Elastic));

        let mut network = ElasticNetwork::new();
        network.connect(0, 1, 1.0, 100.0); // Rest length 1.0, stretched to 2.0

        network.compute_forces(&mut particles, 1.0, 0.0);

        // Particle 0 should be pulled toward particle 1 (positive x)
        let acc0 = particles.get_acceleration(0);
        assert!(acc0[0] > 0.0, "Particle 0 should accelerate in +x direction");
    }
}
