//! Classical STDP
//!
//! Standard spike-timing dependent plasticity with symmetric exponential windows.

#[cfg(feature = "serde")]
use serde::{Deserialize, Serialize};

use crate::traces::{PairDetector, TraceParams};
use crate::weight_bounds::{UpdateAccumulator, WeightBounds, WeightUpdate};
use crate::{PlasticityRule, PlasticityStats};

/// Classical STDP parameters
#[derive(Debug, Clone, Copy)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct ClassicalStdpParams {
    /// LTP amplitude (A+)
    pub a_plus: f32,

    /// LTD amplitude (A-)
    pub a_minus: f32,

    /// LTP time constant (ms)
    pub tau_plus: f32,

    /// LTD time constant (ms)
    pub tau_minus: f32,

    /// Weight bounds
    pub bounds: WeightBounds,
}

impl Default for ClassicalStdpParams {
    fn default() -> Self {
        Self {
            a_plus: 0.005,
            a_minus: 0.00525, // Slightly stronger LTD for stability
            tau_plus: 20.0,
            tau_minus: 20.0,
            bounds: WeightBounds::default(),
        }
    }
}

impl ClassicalStdpParams {
    /// Balanced LTP/LTD
    pub fn balanced() -> Self {
        Self {
            a_plus: 0.005,
            a_minus: 0.005,
            tau_plus: 20.0,
            tau_minus: 20.0,
            bounds: WeightBounds::default(),
        }
    }

    /// LTP-dominated (more potentiation)
    pub fn ltp_dominated() -> Self {
        Self {
            a_plus: 0.01,
            a_minus: 0.005,
            tau_plus: 20.0,
            tau_minus: 20.0,
            bounds: WeightBounds::default(),
        }
    }

    /// LTD-dominated (more depression)
    pub fn ltd_dominated() -> Self {
        Self {
            a_plus: 0.005,
            a_minus: 0.01,
            tau_plus: 20.0,
            tau_minus: 20.0,
            bounds: WeightBounds::default(),
        }
    }

    /// Get trace parameters for LTP window
    fn ltp_trace_params(&self) -> TraceParams {
        TraceParams {
            tau: self.tau_plus,
            increment: 1.0,
            max_value: 10.0,
        }
    }

    /// Get trace parameters for LTD window
    fn ltd_trace_params(&self) -> TraceParams {
        TraceParams {
            tau: self.tau_minus,
            increment: 1.0,
            max_value: 10.0,
        }
    }
}

/// Classical STDP implementation
pub struct ClassicalStdp {
    /// Parameters
    params: ClassicalStdpParams,

    /// Pair detector for timing
    detector: PairDetector,

    /// Update accumulator
    accumulator: UpdateAccumulator,

    /// Synapse to post-neuron mapping
    synapse_to_post: Vec<usize>,

    /// Post-neuron to incoming synapses mapping
    post_to_synapses: Vec<Vec<usize>>,

    /// Number of neurons
    num_neurons: usize,
}

impl ClassicalStdp {
    /// Create new classical STDP rule
    pub fn new(num_synapses: usize, params: ClassicalStdpParams) -> Self {
        // Default mapping: synapse i connects to neuron i % num_neurons
        let num_neurons = (num_synapses / 10).max(1);

        Self::with_mapping(num_synapses, num_neurons, params)
    }

    /// Create with explicit synapse-neuron mapping
    pub fn with_mapping(num_synapses: usize, num_neurons: usize, params: ClassicalStdpParams) -> Self {
        let trace_params = params.ltp_trace_params();

        // Default mapping: distribute synapses evenly
        let synapse_to_post: Vec<usize> = (0..num_synapses)
            .map(|i| i % num_neurons)
            .collect();

        // Build reverse mapping
        let mut post_to_synapses = vec![Vec::new(); num_neurons];
        for (syn_id, &post_id) in synapse_to_post.iter().enumerate() {
            post_to_synapses[post_id].push(syn_id);
        }

        Self {
            params,
            detector: PairDetector::new(num_synapses, num_neurons, trace_params),
            accumulator: UpdateAccumulator::new(num_synapses),
            synapse_to_post,
            post_to_synapses,
            num_neurons,
        }
    }

    /// Set synapse to postsynaptic neuron mapping
    pub fn set_synapse_mapping(&mut self, synapse_id: usize, post_neuron: usize) {
        if synapse_id < self.synapse_to_post.len() {
            // Update reverse mapping
            let old_post = self.synapse_to_post[synapse_id];
            self.post_to_synapses[old_post].retain(|&x| x != synapse_id);

            // Set new mapping
            self.synapse_to_post[synapse_id] = post_neuron;
            if post_neuron < self.post_to_synapses.len() {
                self.post_to_synapses[post_neuron].push(synapse_id);
            }
        }
    }

    /// Get parameters
    pub fn params(&self) -> &ClassicalStdpParams {
        &self.params
    }

    /// Set parameters
    pub fn set_params(&mut self, params: ClassicalStdpParams) {
        self.params = params;
    }
}

impl PlasticityRule for ClassicalStdp {
    fn on_pre_spike(&mut self, synapse_id: usize, time: f64) {
        if synapse_id >= self.synapse_to_post.len() {
            return;
        }

        let post_neuron = self.synapse_to_post[synapse_id];

        // Get post trace (for LTD if post fired recently)
        let post_trace = self.detector.on_pre_spike(synapse_id, post_neuron, time);

        // LTD: post before pre (anti-causal)
        if post_trace > 0.01 {
            let delta = -self.params.a_minus * post_trace;
            self.accumulator.add(synapse_id, delta);
        }
    }

    fn on_post_spike(&mut self, neuron_id: usize, time: f64) {
        if neuron_id >= self.post_to_synapses.len() {
            return;
        }

        let synapse_ids = &self.post_to_synapses[neuron_id];

        // Get pre traces (for LTP if pre fired recently)
        let pre_traces = self.detector.on_post_spike(neuron_id, synapse_ids, time);

        // LTP: pre before post (causal)
        for (&syn_id, &pre_trace) in synapse_ids.iter().zip(pre_traces.iter()) {
            if pre_trace > 0.01 {
                let delta = self.params.a_plus * pre_trace;
                self.accumulator.add(syn_id, delta);
            }
        }
    }

    fn get_weight_updates(&self) -> Vec<WeightUpdate> {
        self.accumulator.get_updates()
    }

    fn apply_updates(&mut self, weights: &mut [f32]) {
        self.accumulator.apply(weights, &self.params.bounds);
        self.accumulator.clear();
    }

    fn reset(&mut self) {
        self.detector.reset();
        self.accumulator.reset();
    }

    fn stats(&self) -> PlasticityStats {
        PlasticityStats {
            ltp_count: self.accumulator.ltp_count(),
            ltd_count: self.accumulator.ltd_count(),
            avg_weight_change: self.accumulator.average_magnitude() as f64,
            max_weight_change: 0.0, // Would need to track
            at_upper_bound: 0,
            at_lower_bound: 0,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_classical_stdp() {
        let mut stdp = ClassicalStdp::new(10, ClassicalStdpParams::default());
        stdp.set_synapse_mapping(0, 0);

        // Pre then post → LTP
        stdp.on_pre_spike(0, 0.0);
        stdp.on_post_spike(0, 10.0);

        let updates = stdp.get_weight_updates();
        assert!(!updates.is_empty());
        assert!(updates[0].delta > 0.0); // Should be LTP
    }

    #[test]
    fn test_ltd() {
        let mut stdp = ClassicalStdp::new(10, ClassicalStdpParams::default());
        stdp.set_synapse_mapping(0, 0);

        // Post then pre → LTD
        stdp.on_post_spike(0, 0.0);
        stdp.on_pre_spike(0, 10.0);

        let updates = stdp.get_weight_updates();
        assert!(!updates.is_empty());
        assert!(updates[0].delta < 0.0); // Should be LTD
    }
}
