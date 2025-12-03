//! Neural-Body Coupling Configuration
//!
//! Defines how neural activity is translated to muscle forces and
//! how body state feeds back into the neural network.

#[cfg(feature = "serde")]
use serde::{Deserialize, Serialize};

/// Coupling configuration
#[derive(Debug, Clone)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct CouplingConfig {
    /// Coupling mode
    pub mode: CouplingMode,

    /// Neural to physics time ratio
    /// e.g., 10 means 10 neural steps per physics step
    pub time_ratio: u32,

    /// Muscle force scaling factor
    pub force_scale: f32,

    /// Proprioceptive gain
    pub proprioceptive_gain: f32,

    /// Enable proprioceptive feedback
    pub proprioception_enabled: bool,

    /// Muscle activation smoothing time constant (ms)
    pub activation_tau: f32,

    /// Proprioceptive delay (ms)
    pub proprioceptive_delay: f32,
}

impl Default for CouplingConfig {
    fn default() -> Self {
        Self {
            mode: CouplingMode::Bidirectional,
            time_ratio: 10,
            force_scale: 1.0,
            proprioceptive_gain: 1.0,
            proprioception_enabled: true,
            activation_tau: 20.0,
            proprioceptive_delay: 5.0,
        }
    }
}

impl CouplingConfig {
    /// Configuration for fast forward locomotion
    pub fn locomotion() -> Self {
        Self {
            mode: CouplingMode::Bidirectional,
            time_ratio: 20,
            force_scale: 1.5,
            proprioceptive_gain: 2.0,
            proprioception_enabled: true,
            activation_tau: 15.0,
            proprioceptive_delay: 3.0,
        }
    }

    /// Open-loop (no proprioception)
    pub fn open_loop() -> Self {
        Self {
            mode: CouplingMode::NeuralToBody,
            time_ratio: 10,
            force_scale: 1.0,
            proprioceptive_gain: 0.0,
            proprioception_enabled: false,
            activation_tau: 20.0,
            proprioceptive_delay: 0.0,
        }
    }

    /// High-fidelity simulation
    pub fn high_fidelity() -> Self {
        Self {
            mode: CouplingMode::Bidirectional,
            time_ratio: 40,
            force_scale: 1.0,
            proprioceptive_gain: 1.0,
            proprioception_enabled: true,
            activation_tau: 10.0,
            proprioceptive_delay: 2.0,
        }
    }
}

/// Coupling mode between neural and body simulations
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub enum CouplingMode {
    /// One-way: neural → body (open-loop)
    NeuralToBody,

    /// One-way: body → neural (passive observation)
    BodyToNeural,

    /// Two-way: neural ↔ body (closed-loop)
    Bidirectional,

    /// Decoupled: both run independently
    Decoupled,
}

impl CouplingMode {
    /// Check if neural output affects body
    pub fn neural_drives_body(&self) -> bool {
        matches!(self, Self::NeuralToBody | Self::Bidirectional)
    }

    /// Check if body state affects neural
    pub fn body_drives_neural(&self) -> bool {
        matches!(self, Self::BodyToNeural | Self::Bidirectional)
    }
}

/// Muscle segment mapping
/// Maps the 96 muscles to body segments for force application
#[derive(Debug, Clone)]
pub struct SegmentMapping {
    /// Muscle indices for each body segment
    /// Each segment has 4 muscles (dorsal-right, ventral-right, ventral-left, dorsal-left)
    pub segment_muscles: [[usize; 4]; 24],

    /// Particle indices that belong to each segment
    pub segment_particles: Vec<Vec<usize>>,

    /// Segment centers (computed from particles)
    pub segment_centers: Vec<[f32; 3]>,
}

impl Default for SegmentMapping {
    fn default() -> Self {
        // Default mapping: muscle index = segment * 4 + quadrant
        let mut segment_muscles = [[0usize; 4]; 24];
        for seg in 0..24 {
            for quad in 0..4 {
                segment_muscles[seg][quad] = seg * 4 + quad;
            }
        }

        Self {
            segment_muscles,
            segment_particles: vec![Vec::new(); 24],
            segment_centers: vec![[0.0; 3]; 24],
        }
    }
}

impl SegmentMapping {
    /// Create mapping from particle positions
    /// Assumes particles are organized along the body axis (x-axis)
    pub fn from_particles(positions: &[[f32; 3]], num_segments: usize) -> Self {
        let mut mapping = Self::default();

        if positions.is_empty() {
            return mapping;
        }

        // Find body extent
        let x_min = positions.iter().map(|p| p[0]).fold(f32::INFINITY, f32::min);
        let x_max = positions.iter().map(|p| p[0]).fold(f32::NEG_INFINITY, f32::max);
        let segment_length = (x_max - x_min) / num_segments as f32;

        // Assign particles to segments
        mapping.segment_particles = vec![Vec::new(); num_segments];

        for (i, pos) in positions.iter().enumerate() {
            let segment = ((pos[0] - x_min) / segment_length) as usize;
            let segment = segment.min(num_segments - 1);
            mapping.segment_particles[segment].push(i);
        }

        // Compute segment centers
        mapping.segment_centers = vec![[0.0; 3]; num_segments];
        for (seg, particles) in mapping.segment_particles.iter().enumerate() {
            if particles.is_empty() {
                continue;
            }

            let mut center = [0.0_f32; 3];
            for &p in particles {
                center[0] += positions[p][0];
                center[1] += positions[p][1];
                center[2] += positions[p][2];
            }

            let n = particles.len() as f32;
            mapping.segment_centers[seg] = [center[0] / n, center[1] / n, center[2] / n];
        }

        mapping
    }

    /// Get muscle indices for a segment
    pub fn get_segment_muscles(&self, segment: usize) -> &[usize; 4] {
        &self.segment_muscles[segment]
    }

    /// Get particles in a segment
    pub fn get_segment_particles(&self, segment: usize) -> &[usize] {
        &self.segment_particles[segment]
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_coupling_modes() {
        let mode = CouplingMode::Bidirectional;
        assert!(mode.neural_drives_body());
        assert!(mode.body_drives_neural());

        let mode = CouplingMode::NeuralToBody;
        assert!(mode.neural_drives_body());
        assert!(!mode.body_drives_neural());
    }

    #[test]
    fn test_segment_mapping() {
        // Create a simple line of particles
        let positions: Vec<[f32; 3]> = (0..48)
            .map(|i| [i as f32 * 0.1, 0.0, 0.0])
            .collect();

        let mapping = SegmentMapping::from_particles(&positions, 24);

        // Should have 2 particles per segment
        for seg in 0..24 {
            assert_eq!(mapping.segment_particles[seg].len(), 2);
        }
    }
}
