//! Muscle Mapping from Neural to Body
//!
//! Maps motor neuron activity to the 96 body wall muscles.

#[cfg(feature = "serde")]
use serde::{Deserialize, Serialize};

use crate::connectome::Connectome;
use crate::neuron::{NeuronId, NeuronState};

/// Neuromuscular junction
#[derive(Debug, Clone, Copy)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct NeuromuscularJunction {
    /// Motor neuron ID
    pub neuron: NeuronId,
    /// Muscle row (0-23)
    pub muscle_row: u32,
    /// Muscle quadrant (0=MDR, 1=MVR, 2=MVL, 3=MDL)
    pub muscle_quadrant: u32,
    /// Synaptic weight
    pub weight: f32,
}

impl NeuromuscularJunction {
    /// Get global muscle index (0-95)
    pub fn muscle_index(&self) -> usize {
        (self.muscle_row * 4 + self.muscle_quadrant) as usize
    }
}

/// Mapping from motor neurons to muscles
#[derive(Debug, Clone, Default)]
pub struct MuscleMap {
    /// Neuromuscular junctions
    nmjs: Vec<NeuromuscularJunction>,
    /// Muscle indices by motor neuron (for fast lookup)
    neuron_to_muscles: Vec<Vec<(usize, f32)>>, // (muscle_idx, weight)
    /// Activation threshold (mV above rest)
    activation_threshold: f32,
    /// Activation gain
    gain: f32,
    /// Time constant for smoothing (ms)
    tau: f32,
    /// Current smoothed activations
    smoothed: [f32; 96],
}

impl MuscleMap {
    /// Create an empty muscle map
    pub fn new() -> Self {
        Self {
            nmjs: Vec::new(),
            neuron_to_muscles: Vec::new(),
            activation_threshold: 10.0, // 10 mV above rest
            gain: 0.1,
            tau: 50.0, // 50 ms smoothing
            smoothed: [0.0; 96],
        }
    }

    /// Create from connectome NMJ data
    pub fn from_connectome(connectome: &Connectome) -> Self {
        let mut map = Self::new();
        let n_neurons = connectome.num_neurons();

        map.neuron_to_muscles = vec![Vec::new(); n_neurons];

        for nmj in connectome.nmjs() {
            let muscle_idx = nmj.muscle_index();
            map.neuron_to_muscles[nmj.neuron as usize].push((muscle_idx, nmj.weight));
            map.nmjs.push(*nmj);
        }

        map
    }

    /// Add a neuromuscular junction
    pub fn add_nmj(&mut self, nmj: NeuromuscularJunction) {
        let neuron = nmj.neuron as usize;

        // Ensure neuron_to_muscles is large enough
        if neuron >= self.neuron_to_muscles.len() {
            self.neuron_to_muscles.resize(neuron + 1, Vec::new());
        }

        let muscle_idx = nmj.muscle_index();
        self.neuron_to_muscles[neuron].push((muscle_idx, nmj.weight));
        self.nmjs.push(nmj);
    }

    /// Compute muscle activation from neuron states
    pub fn compute_activation(&self, neuron_states: &[NeuronState]) -> [f32; 96] {
        let mut activation = [0.0_f32; 96];
        let v_rest = -65.0; // Typical resting potential

        for (neuron_id, muscles) in self.neuron_to_muscles.iter().enumerate() {
            if neuron_id >= neuron_states.len() {
                continue;
            }

            let state = &neuron_states[neuron_id];
            let v_above_rest = state.v - v_rest;

            // Activation is proportional to voltage above threshold
            let neuron_activation = if v_above_rest > self.activation_threshold {
                ((v_above_rest - self.activation_threshold) * self.gain).clamp(0.0, 1.0)
            } else {
                0.0
            };

            // Distribute to muscles
            for &(muscle_idx, weight) in muscles {
                if muscle_idx < 96 {
                    activation[muscle_idx] += neuron_activation * weight;
                }
            }
        }

        // Clamp to valid range
        for a in &mut activation {
            *a = a.clamp(0.0, 1.0);
        }

        activation
    }

    /// Compute activation with temporal smoothing
    pub fn compute_activation_smoothed(
        &mut self,
        neuron_states: &[NeuronState],
        dt: f32,
    ) -> [f32; 96] {
        let instant = self.compute_activation(neuron_states);

        // Exponential moving average
        let alpha = dt / (self.tau + dt);

        for i in 0..96 {
            self.smoothed[i] = alpha * instant[i] + (1.0 - alpha) * self.smoothed[i];
        }

        self.smoothed
    }

    /// Reset smoothed activations
    pub fn reset(&mut self) {
        self.smoothed = [0.0; 96];
    }

    /// Get number of NMJs
    pub fn num_nmjs(&self) -> usize {
        self.nmjs.len()
    }

    /// Get NMJs for a specific motor neuron
    pub fn get_neuron_nmjs(&self, neuron_id: NeuronId) -> &[(usize, f32)] {
        self.neuron_to_muscles.get(neuron_id as usize)
            .map(|v| v.as_slice())
            .unwrap_or(&[])
    }

    /// Set activation parameters
    pub fn set_params(&mut self, threshold: f32, gain: f32, tau: f32) {
        self.activation_threshold = threshold;
        self.gain = gain;
        self.tau = tau;
    }
}

/// Convert 96-element array to 24x4 grid format
pub fn to_grid(flat: &[f32; 96]) -> [[f32; 4]; 24] {
    let mut grid = [[0.0_f32; 4]; 24];
    for row in 0..24 {
        for quad in 0..4 {
            grid[row][quad] = flat[row * 4 + quad];
        }
    }
    grid
}

/// Convert 24x4 grid to 96-element array
pub fn from_grid(grid: &[[f32; 4]; 24]) -> [f32; 96] {
    let mut flat = [0.0_f32; 96];
    for row in 0..24 {
        for quad in 0..4 {
            flat[row * 4 + quad] = grid[row][quad];
        }
    }
    flat
}

/// Muscle quadrant names
pub const QUADRANT_NAMES: [&str; 4] = ["MDR", "MVR", "MVL", "MDL"];

/// Get muscle name from index
pub fn muscle_name(index: usize) -> String {
    if index >= 96 {
        return "INVALID".to_string();
    }
    let row = index / 4;
    let quad = index % 4;
    format!("{}{:02}", QUADRANT_NAMES[quad], row + 1)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_nmj() {
        let nmj = NeuromuscularJunction {
            neuron: 0,
            muscle_row: 5,
            muscle_quadrant: 2,
            weight: 1.0,
        };

        assert_eq!(nmj.muscle_index(), 5 * 4 + 2);
    }

    #[test]
    fn test_muscle_name() {
        assert_eq!(muscle_name(0), "MDR01");
        assert_eq!(muscle_name(1), "MVR01");
        assert_eq!(muscle_name(4), "MDR02");
        assert_eq!(muscle_name(95), "MDL24");
    }

    #[test]
    fn test_grid_conversion() {
        let mut flat = [0.0_f32; 96];
        flat[5] = 0.5;
        flat[20] = 0.8;

        let grid = to_grid(&flat);
        let back = from_grid(&grid);

        assert!((back[5] - 0.5).abs() < 1e-6);
        assert!((back[20] - 0.8).abs() < 1e-6);
    }

    #[test]
    fn test_activation_computation() {
        let mut map = MuscleMap::new();

        // Add a test NMJ
        map.add_nmj(NeuromuscularJunction {
            neuron: 0,
            muscle_row: 0,
            muscle_quadrant: 0,
            weight: 1.0,
        });

        // Create neuron state with high voltage
        let states = vec![NeuronState {
            v: -40.0, // 25 mV above rest
            ..Default::default()
        }];

        let activation = map.compute_activation(&states);

        // Should have positive activation for muscle 0
        assert!(activation[0] > 0.0, "Should activate muscle");
    }
}
