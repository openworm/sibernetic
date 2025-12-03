//! Homeostatic Plasticity
//!
//! Maintains stable neural activity by regulating firing rates and
//! synaptic weight distributions.

#[cfg(feature = "serde")]
use serde::{Deserialize, Serialize};

use crate::weight_bounds::{WeightBounds, WeightUpdate};
use crate::{PlasticityRule, PlasticityStats};

/// Homeostatic plasticity parameters
#[derive(Debug, Clone)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct HomeostaticParams {
    /// Target firing rate (Hz)
    pub target_rate: f32,

    /// Homeostatic learning rate
    pub learning_rate: f32,

    /// Time constant for rate estimation (ms)
    pub tau_rate: f32,

    /// Time constant for homeostatic adjustment (ms)
    pub tau_homeostatic: f32,

    /// Enable intrinsic plasticity (adjust threshold)
    pub intrinsic_plasticity: bool,

    /// Enable synaptic scaling
    pub synaptic_scaling: bool,

    /// Weight bounds
    pub bounds: WeightBounds,
}

impl Default for HomeostaticParams {
    fn default() -> Self {
        Self {
            target_rate: 5.0, // 5 Hz target
            learning_rate: 0.0001,
            tau_rate: 10000.0,      // 10 seconds for rate estimation
            tau_homeostatic: 60000.0, // 1 minute for homeostatic changes
            intrinsic_plasticity: true,
            synaptic_scaling: true,
            bounds: WeightBounds::default(),
        }
    }
}

impl HomeostaticParams {
    /// Fast homeostasis (for quick adaptation)
    pub fn fast() -> Self {
        Self {
            tau_rate: 1000.0,
            tau_homeostatic: 5000.0,
            learning_rate: 0.001,
            ..Default::default()
        }
    }

    /// Slow homeostasis (for stability)
    pub fn slow() -> Self {
        Self {
            tau_rate: 60000.0,
            tau_homeostatic: 300000.0,
            learning_rate: 0.00001,
            ..Default::default()
        }
    }
}

/// Per-neuron homeostatic state
#[derive(Debug, Clone, Default)]
struct NeuronHomeostasis {
    /// Estimated firing rate (Hz)
    rate_estimate: f32,

    /// Spike count for rate estimation
    spike_count: u32,

    /// Time window start (ms)
    window_start: f64,

    /// Intrinsic excitability modifier
    excitability: f32,

    /// Scaling factor for incoming synapses
    scaling_factor: f32,
}

impl NeuronHomeostasis {
    fn new() -> Self {
        Self {
            rate_estimate: 0.0,
            spike_count: 0,
            window_start: 0.0,
            excitability: 1.0,
            scaling_factor: 1.0,
        }
    }
}

/// Homeostatic plasticity implementation
pub struct HomeostaticStdp {
    /// Parameters
    params: HomeostaticParams,

    /// Per-neuron homeostatic state
    neuron_states: Vec<NeuronHomeostasis>,

    /// Post-neuron to incoming synapses mapping
    post_to_synapses: Vec<Vec<usize>>,

    /// Pending weight updates
    pending_updates: Vec<WeightUpdate>,

    /// Last update time
    last_update: f64,

    /// Statistics
    ltp_count: u64,
    ltd_count: u64,
}

impl HomeostaticStdp {
    /// Create new homeostatic plasticity rule
    pub fn new(num_neurons: usize, num_synapses: usize, params: HomeostaticParams) -> Self {
        // Default synapse distribution
        let mut post_to_synapses = vec![Vec::new(); num_neurons];
        for syn_id in 0..num_synapses {
            let post_id = syn_id % num_neurons;
            post_to_synapses[post_id].push(syn_id);
        }

        Self {
            params,
            neuron_states: vec![NeuronHomeostasis::new(); num_neurons],
            post_to_synapses,
            pending_updates: Vec::new(),
            last_update: 0.0,
            ltp_count: 0,
            ltd_count: 0,
        }
    }

    /// Update homeostatic variables at given time
    pub fn update(&mut self, time: f64) {
        if time <= self.last_update {
            return;
        }

        let dt = (time - self.last_update) as f32;
        let rate_decay = (-dt / self.params.tau_rate).exp();
        let homeostatic_rate = dt / self.params.tau_homeostatic;

        for (neuron_id, state) in self.neuron_states.iter_mut().enumerate() {
            // Update rate estimate
            let window_duration = (time - state.window_start) as f32;
            if window_duration > 100.0 {
                // Minimum 100ms window
                state.rate_estimate = state.rate_estimate * rate_decay
                    + (1.0 - rate_decay) * (state.spike_count as f32 * 1000.0 / window_duration);
                state.spike_count = 0;
                state.window_start = time;
            }

            // Compute rate error
            let rate_error = self.params.target_rate - state.rate_estimate;

            // Intrinsic plasticity: adjust excitability
            if self.params.intrinsic_plasticity {
                let delta_exc = rate_error * self.params.learning_rate * homeostatic_rate;
                state.excitability += delta_exc;
                state.excitability = state.excitability.clamp(0.5, 2.0);
            }

            // Synaptic scaling: adjust all incoming weights
            if self.params.synaptic_scaling {
                let delta_scale = rate_error * self.params.learning_rate * homeostatic_rate * 0.1;
                state.scaling_factor += delta_scale;
                state.scaling_factor = state.scaling_factor.clamp(0.5, 2.0);

                // Generate weight updates for incoming synapses
                if neuron_id < self.post_to_synapses.len() {
                    for &syn_id in &self.post_to_synapses[neuron_id] {
                        if delta_scale.abs() > 1e-9 {
                            self.pending_updates.push(WeightUpdate::new(syn_id, delta_scale));
                            if delta_scale > 0.0 {
                                self.ltp_count += 1;
                            } else {
                                self.ltd_count += 1;
                            }
                        }
                    }
                }
            }
        }

        self.last_update = time;
    }

    /// Get excitability modifier for a neuron
    pub fn get_excitability(&self, neuron_id: usize) -> f32 {
        self.neuron_states
            .get(neuron_id)
            .map(|s| s.excitability)
            .unwrap_or(1.0)
    }

    /// Get scaling factor for a neuron
    pub fn get_scaling_factor(&self, neuron_id: usize) -> f32 {
        self.neuron_states
            .get(neuron_id)
            .map(|s| s.scaling_factor)
            .unwrap_or(1.0)
    }

    /// Get estimated firing rate for a neuron
    pub fn get_rate_estimate(&self, neuron_id: usize) -> f32 {
        self.neuron_states
            .get(neuron_id)
            .map(|s| s.rate_estimate)
            .unwrap_or(0.0)
    }

    /// Set synapse to postsynaptic neuron mapping
    pub fn set_synapse_mapping(&mut self, synapse_id: usize, post_neuron: usize) {
        // Remove from old post
        for syns in &mut self.post_to_synapses {
            syns.retain(|&x| x != synapse_id);
        }

        // Add to new post
        if post_neuron < self.post_to_synapses.len() {
            self.post_to_synapses[post_neuron].push(synapse_id);
        }
    }
}

impl PlasticityRule for HomeostaticStdp {
    fn on_pre_spike(&mut self, _synapse_id: usize, _time: f64) {
        // Homeostatic plasticity doesn't directly respond to pre spikes
    }

    fn on_post_spike(&mut self, neuron_id: usize, time: f64) {
        if neuron_id < self.neuron_states.len() {
            self.neuron_states[neuron_id].spike_count += 1;
        }

        // Periodically update homeostatic variables
        if time - self.last_update > 100.0 {
            self.update(time);
        }
    }

    fn get_weight_updates(&self) -> Vec<WeightUpdate> {
        self.pending_updates.clone()
    }

    fn apply_updates(&mut self, weights: &mut [f32]) {
        for update in &self.pending_updates {
            if update.synapse_id < weights.len() {
                let new_weight = weights[update.synapse_id] + update.delta;
                weights[update.synapse_id] = self.params.bounds.clamp(new_weight);
            }
        }
        self.pending_updates.clear();
    }

    fn reset(&mut self) {
        for state in &mut self.neuron_states {
            *state = NeuronHomeostasis::new();
        }
        self.pending_updates.clear();
        self.last_update = 0.0;
        self.ltp_count = 0;
        self.ltd_count = 0;
    }

    fn stats(&self) -> PlasticityStats {
        PlasticityStats {
            ltp_count: self.ltp_count,
            ltd_count: self.ltd_count,
            avg_weight_change: 0.0,
            max_weight_change: 0.0,
            at_upper_bound: 0,
            at_lower_bound: 0,
        }
    }
}

/// Synaptic scaling function
/// Multiplicatively scales all synapses to maintain target activity
pub fn apply_synaptic_scaling(
    weights: &mut [f32],
    neuron_rates: &[f32],
    target_rate: f32,
    learning_rate: f32,
    neuron_to_synapses: &[Vec<usize>],
) {
    for (neuron_id, synapses) in neuron_to_synapses.iter().enumerate() {
        if neuron_id >= neuron_rates.len() {
            continue;
        }

        let rate = neuron_rates[neuron_id];
        let scaling_factor = 1.0 + learning_rate * (target_rate - rate) / target_rate;

        for &syn_id in synapses {
            if syn_id < weights.len() {
                weights[syn_id] *= scaling_factor;
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_homeostatic_creation() {
        let stdp = HomeostaticStdp::new(10, 100, HomeostaticParams::default());
        assert_eq!(stdp.neuron_states.len(), 10);
    }

    #[test]
    fn test_rate_adaptation() {
        let mut stdp = HomeostaticStdp::new(1, 10, HomeostaticParams::fast());

        // Simulate high firing rate
        for i in 0..100 {
            stdp.on_post_spike(0, i as f64 * 10.0); // 100 Hz
        }
        stdp.update(1000.0);

        // Excitability should decrease (rate too high)
        assert!(stdp.get_excitability(0) < 1.0);
    }

    #[test]
    fn test_synaptic_scaling() {
        let mut weights = vec![1.0; 10];
        let rates = vec![10.0]; // Above target
        let target = 5.0;
        let neuron_to_synapses = vec![(0..10).collect()];

        apply_synaptic_scaling(&mut weights, &rates, target, 0.1, &neuron_to_synapses);

        // Weights should decrease (rate too high)
        assert!(weights[0] < 1.0);
    }
}
