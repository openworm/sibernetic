//! Muscle Simulation for Biomechanical Models
//!
//! Implements muscle activation and force generation for C. elegans
//! and other biomechanical simulations.

use serde::{Deserialize, Serialize};

/// Muscle quadrant in C. elegans body wall
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum MuscleQuadrant {
    /// Dorsal Right
    MDR,
    /// Ventral Right
    MVR,
    /// Ventral Left
    MVL,
    /// Dorsal Left
    MDL,
}

impl MuscleQuadrant {
    /// Get all quadrants in order
    pub fn all() -> [Self; 4] {
        [Self::MDR, Self::MVR, Self::MVL, Self::MDL]
    }

    /// Get quadrant index (0-3)
    pub fn index(&self) -> usize {
        match self {
            Self::MDR => 0,
            Self::MVR => 1,
            Self::MVL => 2,
            Self::MDL => 3,
        }
    }

    /// Create from index
    pub fn from_index(idx: usize) -> Option<Self> {
        match idx % 4 {
            0 => Some(Self::MDR),
            1 => Some(Self::MVR),
            2 => Some(Self::MVL),
            3 => Some(Self::MDL),
            _ => None,
        }
    }
}

/// A single muscle segment
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct MuscleSegment {
    /// Muscle name (e.g., "MDR01", "MVL12")
    pub name: String,
    /// Quadrant
    pub quadrant: MuscleQuadrant,
    /// Row index (0-23 for C. elegans)
    pub row: usize,
    /// Global muscle index
    pub index: usize,
    /// Current activation level (0.0 to 1.0)
    pub activation: f32,
    /// Associated elastic connection indices
    pub connections: Vec<usize>,
}

impl MuscleSegment {
    /// Create a new muscle segment
    pub fn new(quadrant: MuscleQuadrant, row: usize) -> Self {
        let index = row * 4 + quadrant.index();
        let name = format!(
            "{}{:02}",
            match quadrant {
                MuscleQuadrant::MDR => "MDR",
                MuscleQuadrant::MVR => "MVR",
                MuscleQuadrant::MVL => "MVL",
                MuscleQuadrant::MDL => "MDL",
            },
            row + 1
        );

        Self {
            name,
            quadrant,
            row,
            index,
            activation: 0.0,
            connections: Vec::new(),
        }
    }

    /// Add an elastic connection to this muscle
    pub fn add_connection(&mut self, connection_idx: usize) {
        self.connections.push(connection_idx);
    }

    /// Set activation level
    pub fn set_activation(&mut self, activation: f32) {
        self.activation = activation.clamp(0.0, 1.0);
    }
}

/// Muscle activation pattern for the full body
#[derive(Debug, Clone, Default, Serialize, Deserialize)]
pub struct MuscleActivation {
    /// 24 rows x 4 quadrants = 96 muscle segments
    pub activations: [[f32; 4]; 24],
}

impl MuscleActivation {
    /// Create a new muscle activation pattern (all zeros)
    pub fn new() -> Self {
        Self {
            activations: [[0.0; 4]; 24],
        }
    }

    /// Set activation for a specific muscle
    pub fn set(&mut self, quadrant: MuscleQuadrant, row: usize, activation: f32) {
        if row < 24 {
            self.activations[row][quadrant.index()] = activation.clamp(0.0, 1.0);
        }
    }

    /// Get activation for a specific muscle
    pub fn get(&self, quadrant: MuscleQuadrant, row: usize) -> f32 {
        if row < 24 {
            self.activations[row][quadrant.index()]
        } else {
            0.0
        }
    }

    /// Set activation by global index (0-95)
    pub fn set_by_index(&mut self, index: usize, activation: f32) {
        if index < 96 {
            let row = index / 4;
            let quad = index % 4;
            self.activations[row][quad] = activation.clamp(0.0, 1.0);
        }
    }

    /// Get activation by global index
    pub fn get_by_index(&self, index: usize) -> f32 {
        if index < 96 {
            let row = index / 4;
            let quad = index % 4;
            self.activations[row][quad]
        } else {
            0.0
        }
    }

    /// Get flattened array of all activations
    pub fn to_flat(&self) -> [f32; 96] {
        let mut flat = [0.0; 96];
        for (row, row_activations) in self.activations.iter().enumerate() {
            for (quad, &activation) in row_activations.iter().enumerate() {
                flat[row * 4 + quad] = activation;
            }
        }
        flat
    }

    /// Set from flattened array
    pub fn from_flat(flat: &[f32]) -> Self {
        let mut result = Self::new();
        for (i, &activation) in flat.iter().take(96).enumerate() {
            result.set_by_index(i, activation);
        }
        result
    }

    /// Generate swimming wave pattern
    ///
    /// Based on Sibernetic's main_sim.py parallel_waves function
    pub fn swimming_wave(step: u64, velocity: f32, amplitude: f32) -> Self {
        use std::f32::consts::PI;

        let mut activation = Self::new();
        let n = 12; // Half the rows
        let phase = velocity * step as f32;

        // Wave modulation (amplitude decreases toward tail)
        let wave_m = [0.81, 0.90, 0.97, 1.00, 0.99, 0.95, 0.88, 0.78, 0.65, 0.53, 0.40, 0.25];

        for i in 0..n {
            let row_pos = i as f32 * 0.81 * PI / n as f32;

            // Wave 1 (dorsal right + ventral left)
            let w1 = ((row_pos - phase).sin()).max(0.0) * wave_m[i] * amplitude;
            // Wave 2 (dorsal left + ventral right) - phase shifted
            let w2 = ((row_pos - phase + PI).sin()).max(0.0) * wave_m[i] * amplitude;

            // Front half
            activation.set(MuscleQuadrant::MDR, i, w1);
            activation.set(MuscleQuadrant::MVL, i, w1);
            activation.set(MuscleQuadrant::MDL, i, w2);
            activation.set(MuscleQuadrant::MVR, i, w2);

            // Back half (mirror)
            activation.set(MuscleQuadrant::MDR, i + 12, w2);
            activation.set(MuscleQuadrant::MVL, i + 12, w2);
            activation.set(MuscleQuadrant::MDL, i + 12, w1);
            activation.set(MuscleQuadrant::MVR, i + 12, w1);
        }

        activation
    }

    /// Generate crawling wave pattern
    pub fn crawling_wave(step: u64, velocity: f32, amplitude: f32) -> Self {
        use std::f32::consts::PI;

        let mut activation = Self::new();
        let n = 12;
        let phase = velocity * step as f32;

        // Crawling has higher frequency, different modulation
        let wave_m: [f32; 12] = [1.0, 0.95, 0.90, 0.85, 0.80, 0.75, 0.70, 0.65, 0.60, 0.55, 0.50, 0.45];

        for i in 0..n {
            let row_pos = i as f32 * 2.97 * PI / n as f32;

            let w1 = ((row_pos - phase).sin()).max(0.0) * wave_m[i] * amplitude;
            let w2 = ((row_pos - phase + PI).sin()).max(0.0) * wave_m[i] * amplitude;

            // Front half
            activation.set(MuscleQuadrant::MDR, i, w1);
            activation.set(MuscleQuadrant::MVL, i, w1);
            activation.set(MuscleQuadrant::MDL, i, w2);
            activation.set(MuscleQuadrant::MVR, i, w2);

            // Back half
            activation.set(MuscleQuadrant::MDR, i + 12, w2);
            activation.set(MuscleQuadrant::MVL, i + 12, w2);
            activation.set(MuscleQuadrant::MDL, i + 12, w1);
            activation.set(MuscleQuadrant::MVR, i + 12, w1);
        }

        activation
    }

    /// Smooth start (ramp up over time)
    pub fn with_smooth_start(&mut self, step: u64, ramp_steps: u64) {
        if step < ramp_steps {
            let factor = step as f32 / ramp_steps as f32;
            for row in &mut self.activations {
                for activation in row {
                    *activation *= factor;
                }
            }
        }
    }
}

/// Get list of muscle names in standard order
pub fn get_muscle_names() -> Vec<String> {
    let mut names = Vec::with_capacity(96);
    for row in 0..24 {
        for quadrant in MuscleQuadrant::all() {
            let name = format!(
                "{}{:02}",
                match quadrant {
                    MuscleQuadrant::MDR => "MDR",
                    MuscleQuadrant::MVR => "MVR",
                    MuscleQuadrant::MVL => "MVL",
                    MuscleQuadrant::MDL => "MDL",
                },
                row + 1
            );
            names.push(name);
        }
    }
    names
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_muscle_quadrant() {
        assert_eq!(MuscleQuadrant::MDR.index(), 0);
        assert_eq!(MuscleQuadrant::MVR.index(), 1);
        assert_eq!(MuscleQuadrant::MVL.index(), 2);
        assert_eq!(MuscleQuadrant::MDL.index(), 3);
    }

    #[test]
    fn test_muscle_activation() {
        let mut activation = MuscleActivation::new();

        activation.set(MuscleQuadrant::MDR, 0, 0.5);
        assert!((activation.get(MuscleQuadrant::MDR, 0) - 0.5).abs() < 1e-6);

        activation.set_by_index(4, 0.7); // Row 1, quad 0
        assert!((activation.get(MuscleQuadrant::MDR, 1) - 0.7).abs() < 1e-6);
    }

    #[test]
    fn test_swimming_wave() {
        let activation = MuscleActivation::swimming_wave(0, 0.0001, 1.0);

        // Check that some muscles are activated
        let flat = activation.to_flat();
        let total: f32 = flat.iter().sum();
        assert!(total > 0.0, "Swimming wave should activate some muscles");
    }

    #[test]
    fn test_muscle_names() {
        let names = get_muscle_names();
        assert_eq!(names.len(), 96);
        assert_eq!(names[0], "MDR01");
        assert_eq!(names[3], "MDL01");
        assert_eq!(names[4], "MDR02");
    }
}
