//! Triplet STDP
//!
//! Three-factor STDP rule that considers triplets of spikes for more
//! realistic frequency-dependence of plasticity.

#[cfg(feature = "serde")]
use serde::{Deserialize, Serialize};

use crate::traces::{SpikeTrace, TraceParams};
use crate::weight_bounds::{UpdateAccumulator, WeightBounds, WeightUpdate};
use crate::{PlasticityRule, PlasticityStats};

/// Triplet STDP parameters
#[derive(Debug, Clone, Copy)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct TripletStdpParams {
    /// Fast LTP amplitude (pair term)
    pub a2_plus: f32,

    /// Slow LTP amplitude (triplet term)
    pub a3_plus: f32,

    /// Fast LTD amplitude (pair term)
    pub a2_minus: f32,

    /// Slow LTD amplitude (triplet term)
    pub a3_minus: f32,

    /// Fast pre-synaptic trace time constant (ms)
    pub tau_plus: f32,

    /// Slow pre-synaptic trace time constant (ms)
    pub tau_x: f32,

    /// Fast post-synaptic trace time constant (ms)
    pub tau_minus: f32,

    /// Slow post-synaptic trace time constant (ms)
    pub tau_y: f32,

    /// Weight bounds
    pub bounds: WeightBounds,
}

impl Default for TripletStdpParams {
    fn default() -> Self {
        // Parameters from Pfister & Gerstner 2006
        Self {
            a2_plus: 0.0046,
            a3_plus: 0.0091,
            a2_minus: 0.003,
            a3_minus: 0.0,
            tau_plus: 16.8,
            tau_x: 101.0,
            tau_minus: 33.7,
            tau_y: 114.0,
            bounds: WeightBounds::default(),
        }
    }
}

impl TripletStdpParams {
    /// Minimal triplet (pair-based only)
    pub fn minimal() -> Self {
        Self {
            a2_plus: 0.005,
            a3_plus: 0.0,
            a2_minus: 0.005,
            a3_minus: 0.0,
            ..Default::default()
        }
    }

    /// Full triplet model
    pub fn full() -> Self {
        Self::default()
    }
}

/// Per-synapse triplet state
#[derive(Debug, Clone, Default)]
struct SynapseState {
    /// Fast pre trace (r1)
    r1: SpikeTrace,
    /// Slow pre trace (r2)
    r2: SpikeTrace,
}

/// Per-neuron triplet state
#[derive(Debug, Clone, Default)]
struct NeuronState {
    /// Fast post trace (o1)
    o1: SpikeTrace,
    /// Slow post trace (o2)
    o2: SpikeTrace,
}

/// Triplet STDP implementation
pub struct TripletStdp {
    /// Parameters
    params: TripletStdpParams,

    /// Per-synapse state
    synapse_states: Vec<SynapseState>,

    /// Per-neuron state
    neuron_states: Vec<NeuronState>,

    /// Update accumulator
    accumulator: UpdateAccumulator,

    /// Synapse to post-neuron mapping
    synapse_to_post: Vec<usize>,

    /// Post-neuron to incoming synapses mapping
    post_to_synapses: Vec<Vec<usize>>,

    /// Trace parameters for different time scales
    trace_params_fast_pre: TraceParams,
    trace_params_slow_pre: TraceParams,
    trace_params_fast_post: TraceParams,
    trace_params_slow_post: TraceParams,
}

impl TripletStdp {
    /// Create new triplet STDP rule
    pub fn new(num_synapses: usize, num_neurons: usize, params: TripletStdpParams) -> Self {
        // Default synapse mapping
        let synapse_to_post: Vec<usize> = (0..num_synapses).map(|i| i % num_neurons).collect();

        let mut post_to_synapses = vec![Vec::new(); num_neurons];
        for (syn_id, &post_id) in synapse_to_post.iter().enumerate() {
            post_to_synapses[post_id].push(syn_id);
        }

        let trace_params_fast_pre = TraceParams {
            tau: params.tau_plus,
            increment: 1.0,
            max_value: 10.0,
        };
        let trace_params_slow_pre = TraceParams {
            tau: params.tau_x,
            increment: 1.0,
            max_value: 10.0,
        };
        let trace_params_fast_post = TraceParams {
            tau: params.tau_minus,
            increment: 1.0,
            max_value: 10.0,
        };
        let trace_params_slow_post = TraceParams {
            tau: params.tau_y,
            increment: 1.0,
            max_value: 10.0,
        };

        Self {
            params,
            synapse_states: vec![SynapseState::default(); num_synapses],
            neuron_states: vec![NeuronState::default(); num_neurons],
            accumulator: UpdateAccumulator::new(num_synapses),
            synapse_to_post,
            post_to_synapses,
            trace_params_fast_pre,
            trace_params_slow_pre,
            trace_params_fast_post,
            trace_params_slow_post,
        }
    }

    /// Set synapse to postsynaptic neuron mapping
    pub fn set_synapse_mapping(&mut self, synapse_id: usize, post_neuron: usize) {
        if synapse_id < self.synapse_to_post.len() {
            let old_post = self.synapse_to_post[synapse_id];
            self.post_to_synapses[old_post].retain(|&x| x != synapse_id);

            self.synapse_to_post[synapse_id] = post_neuron;
            if post_neuron < self.post_to_synapses.len() {
                self.post_to_synapses[post_neuron].push(synapse_id);
            }
        }
    }
}

impl PlasticityRule for TripletStdp {
    fn on_pre_spike(&mut self, synapse_id: usize, time: f64) {
        if synapse_id >= self.synapse_states.len() {
            return;
        }

        let post_neuron = self.synapse_to_post[synapse_id];
        if post_neuron >= self.neuron_states.len() {
            return;
        }

        // Get post traces at spike time (for LTD)
        let o1 = self.neuron_states[post_neuron]
            .o1
            .value_at(time, &self.trace_params_fast_post);
        let o2 = self.neuron_states[post_neuron]
            .o2
            .value_at(time, &self.trace_params_slow_post);

        // LTD: pair term + triplet term
        let ltd = self.params.a2_minus * o1 + self.params.a3_minus * o1 * o2;
        if ltd > 1e-9 {
            self.accumulator.add(synapse_id, -ltd);
        }

        // Update pre traces
        self.synapse_states[synapse_id]
            .r1
            .on_spike(time, &self.trace_params_fast_pre);
        self.synapse_states[synapse_id]
            .r2
            .on_spike(time, &self.trace_params_slow_pre);
    }

    fn on_post_spike(&mut self, neuron_id: usize, time: f64) {
        if neuron_id >= self.post_to_synapses.len() {
            return;
        }

        // Get slow post trace before updating (for triplet LTP)
        let o2_before = self.neuron_states[neuron_id]
            .o2
            .value_at(time, &self.trace_params_slow_post);

        // LTP for all incoming synapses
        for &syn_id in &self.post_to_synapses[neuron_id] {
            let r1 = self.synapse_states[syn_id]
                .r1
                .value_at(time, &self.trace_params_fast_pre);
            let r2 = self.synapse_states[syn_id]
                .r2
                .value_at(time, &self.trace_params_slow_pre);

            // LTP: pair term + triplet term
            let ltp = self.params.a2_plus * r1 + self.params.a3_plus * r1 * o2_before;
            if ltp > 1e-9 {
                self.accumulator.add(syn_id, ltp);
            }
        }

        // Update post traces
        self.neuron_states[neuron_id]
            .o1
            .on_spike(time, &self.trace_params_fast_post);
        self.neuron_states[neuron_id]
            .o2
            .on_spike(time, &self.trace_params_slow_post);
    }

    fn get_weight_updates(&self) -> Vec<WeightUpdate> {
        self.accumulator.get_updates()
    }

    fn apply_updates(&mut self, weights: &mut [f32]) {
        self.accumulator.apply(weights, &self.params.bounds);
        self.accumulator.clear();
    }

    fn reset(&mut self) {
        for state in &mut self.synapse_states {
            state.r1.reset();
            state.r2.reset();
        }
        for state in &mut self.neuron_states {
            state.o1.reset();
            state.o2.reset();
        }
        self.accumulator.reset();
    }

    fn stats(&self) -> PlasticityStats {
        PlasticityStats {
            ltp_count: self.accumulator.ltp_count(),
            ltd_count: self.accumulator.ltd_count(),
            avg_weight_change: self.accumulator.average_magnitude() as f64,
            max_weight_change: 0.0,
            at_upper_bound: 0,
            at_lower_bound: 0,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_triplet_stdp() {
        let mut stdp = TripletStdp::new(10, 5, TripletStdpParams::default());
        stdp.set_synapse_mapping(0, 0);

        // Pre then post → LTP
        stdp.on_pre_spike(0, 0.0);
        stdp.on_post_spike(0, 10.0);

        let updates = stdp.get_weight_updates();
        assert!(!updates.is_empty());
        assert!(updates[0].delta > 0.0);
    }

    #[test]
    fn test_triplet_frequency() {
        let mut stdp = TripletStdp::new(1, 1, TripletStdpParams::default());
        stdp.set_synapse_mapping(0, 0);

        // High frequency pairing should produce more LTP
        for i in 0..10 {
            let t = i as f64 * 20.0;
            stdp.on_pre_spike(0, t);
            stdp.on_post_spike(0, t + 5.0);
        }

        let updates = stdp.get_weight_updates();
        let total_ltp: f32 = updates.iter().map(|u| u.delta.max(0.0)).sum();

        assert!(total_ltp > 0.01);
    }
}
