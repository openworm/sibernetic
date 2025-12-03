//! Reward-Modulated STDP
//!
//! Three-factor learning rule combining spike timing with reward signals.
//! Used for reinforcement learning in spiking networks.

#[cfg(feature = "serde")]
use serde::{Deserialize, Serialize};

use crate::traces::{EligibilityTrace, PairDetector, TraceParams};
use crate::weight_bounds::{UpdateAccumulator, WeightBounds, WeightUpdate};
use crate::{PlasticityRule, PlasticityStats};

/// Reward signal for modulating plasticity
#[derive(Debug, Clone, Copy)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct RewardSignal {
    /// Reward value (positive = reward, negative = punishment)
    pub value: f32,

    /// Time of reward delivery (ms)
    pub time: f64,

    /// Whether reward is dopamine-like (phasic) or sustained
    pub phasic: bool,
}

impl Default for RewardSignal {
    fn default() -> Self {
        Self {
            value: 0.0,
            time: 0.0,
            phasic: true,
        }
    }
}

impl RewardSignal {
    /// Create reward signal
    pub fn reward(value: f32, time: f64) -> Self {
        Self {
            value,
            time,
            phasic: true,
        }
    }

    /// Create punishment signal
    pub fn punishment(value: f32, time: f64) -> Self {
        Self {
            value: -value.abs(),
            time,
            phasic: true,
        }
    }

    /// Check if this is a reward (positive)
    pub fn is_reward(&self) -> bool {
        self.value > 0.0
    }

    /// Check if this is a punishment (negative)
    pub fn is_punishment(&self) -> bool {
        self.value < 0.0
    }
}

/// Reward-modulated STDP parameters
#[derive(Debug, Clone)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct RewardModulatedParams {
    /// Learning rate for eligibility
    pub learning_rate: f32,

    /// LTP amplitude
    pub a_plus: f32,

    /// LTD amplitude
    pub a_minus: f32,

    /// Eligibility trace time constant (ms)
    pub tau_eligibility: f32,

    /// Spike timing trace time constant (ms)
    pub tau_timing: f32,

    /// Dopamine time constant (ms) - for phasic reward decay
    pub tau_dopamine: f32,

    /// Baseline dopamine level
    pub dopamine_baseline: f32,

    /// Weight bounds
    pub bounds: WeightBounds,
}

impl Default for RewardModulatedParams {
    fn default() -> Self {
        Self {
            learning_rate: 0.01,
            a_plus: 1.0,
            a_minus: 1.0,
            tau_eligibility: 1000.0, // 1 second
            tau_timing: 20.0,
            tau_dopamine: 200.0,
            dopamine_baseline: 0.0,
            bounds: WeightBounds::default(),
        }
    }
}

impl RewardModulatedParams {
    /// Fast learning (short eligibility window)
    pub fn fast() -> Self {
        Self {
            tau_eligibility: 200.0,
            learning_rate: 0.05,
            ..Default::default()
        }
    }

    /// Slow learning (long eligibility window)
    pub fn slow() -> Self {
        Self {
            tau_eligibility: 5000.0,
            learning_rate: 0.001,
            ..Default::default()
        }
    }
}

/// Reward-modulated STDP implementation
pub struct RewardModulatedStdp {
    /// Parameters
    params: RewardModulatedParams,

    /// Pair detector for spike timing
    detector: PairDetector,

    /// Eligibility traces
    eligibility: EligibilityTrace,

    /// Update accumulator
    accumulator: UpdateAccumulator,

    /// Current dopamine level
    dopamine: f32,

    /// Last dopamine update time
    last_dopamine_time: f64,

    /// Synapse to post-neuron mapping
    synapse_to_post: Vec<usize>,

    /// Post-neuron to incoming synapses mapping
    post_to_synapses: Vec<Vec<usize>>,

    /// Total reward received
    total_reward: f64,

    /// Number of reward events
    reward_count: u64,
}

impl RewardModulatedStdp {
    /// Create new reward-modulated STDP rule
    pub fn new(num_synapses: usize, params: RewardModulatedParams) -> Self {
        let num_neurons = (num_synapses / 10).max(1);
        Self::with_mapping(num_synapses, num_neurons, params)
    }

    /// Create with explicit mapping
    pub fn with_mapping(
        num_synapses: usize,
        num_neurons: usize,
        params: RewardModulatedParams,
    ) -> Self {
        let timing_params = TraceParams {
            tau: params.tau_timing,
            increment: 1.0,
            max_value: 10.0,
        };

        let eligibility_params = TraceParams {
            tau: params.tau_eligibility,
            increment: 1.0,
            max_value: 1.0,
        };

        let synapse_to_post: Vec<usize> = (0..num_synapses).map(|i| i % num_neurons).collect();

        let mut post_to_synapses = vec![Vec::new(); num_neurons];
        for (syn_id, &post_id) in synapse_to_post.iter().enumerate() {
            post_to_synapses[post_id].push(syn_id);
        }

        Self {
            params: params.clone(),
            detector: PairDetector::new(num_synapses, num_neurons, timing_params),
            eligibility: EligibilityTrace::new(num_synapses, eligibility_params),
            accumulator: UpdateAccumulator::new(num_synapses),
            dopamine: params.dopamine_baseline,
            last_dopamine_time: 0.0,
            synapse_to_post,
            post_to_synapses,
            total_reward: 0.0,
            reward_count: 0,
        }
    }

    /// Deliver reward signal
    pub fn deliver_reward(&mut self, reward: RewardSignal) {
        // Decay dopamine to current time
        self.update_dopamine(reward.time);

        if reward.phasic {
            // Phasic dopamine burst
            self.dopamine += reward.value;
        } else {
            // Sustained change in dopamine baseline
            self.dopamine = reward.value;
        }

        // Apply reward to eligible synapses
        let updates = self.eligibility.apply_reward(
            self.dopamine - self.params.dopamine_baseline,
            self.params.learning_rate,
        );

        for (syn_id, delta) in updates {
            self.accumulator.add(syn_id, delta);
        }

        self.total_reward += reward.value as f64;
        self.reward_count += 1;
    }

    /// Update dopamine level (decay toward baseline)
    fn update_dopamine(&mut self, time: f64) {
        if time > self.last_dopamine_time {
            let dt = (time - self.last_dopamine_time) as f32;
            let decay = (-dt / self.params.tau_dopamine).exp();
            self.dopamine = self.params.dopamine_baseline
                + (self.dopamine - self.params.dopamine_baseline) * decay;
            self.last_dopamine_time = time;
        }
    }

    /// Get current dopamine level
    pub fn dopamine_level(&self) -> f32 {
        self.dopamine
    }

    /// Get eligibility for a synapse
    pub fn get_eligibility(&self, synapse_id: usize) -> f32 {
        self.eligibility.get(synapse_id)
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

impl PlasticityRule for RewardModulatedStdp {
    fn on_pre_spike(&mut self, synapse_id: usize, time: f64) {
        if synapse_id >= self.synapse_to_post.len() {
            return;
        }

        let post_neuron = self.synapse_to_post[synapse_id];

        // Get post trace (for eligibility if anti-causal)
        let post_trace = self.detector.on_pre_spike(synapse_id, post_neuron, time);

        // Anti-causal pairing creates negative eligibility (LTD candidate)
        if post_trace > 0.01 {
            self.eligibility.add(synapse_id, -self.params.a_minus * post_trace, time);
        }
    }

    fn on_post_spike(&mut self, neuron_id: usize, time: f64) {
        if neuron_id >= self.post_to_synapses.len() {
            return;
        }

        let synapse_ids = &self.post_to_synapses[neuron_id];

        // Get pre traces (for eligibility if causal)
        let pre_traces = self.detector.on_post_spike(neuron_id, synapse_ids, time);

        // Causal pairing creates positive eligibility (LTP candidate)
        for (&syn_id, &pre_trace) in synapse_ids.iter().zip(pre_traces.iter()) {
            if pre_trace > 0.01 {
                self.eligibility.add(syn_id, self.params.a_plus * pre_trace, time);
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
        self.eligibility.reset();
        self.accumulator.reset();
        self.dopamine = self.params.dopamine_baseline;
        self.last_dopamine_time = 0.0;
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
    fn test_reward_stdp() {
        let mut stdp = RewardModulatedStdp::new(10, RewardModulatedParams::default());
        stdp.set_synapse_mapping(0, 0);

        // Create eligibility through spike pairing
        stdp.on_pre_spike(0, 0.0);
        stdp.on_post_spike(0, 10.0);

        assert!(stdp.get_eligibility(0) > 0.0);

        // Deliver reward
        stdp.deliver_reward(RewardSignal::reward(1.0, 50.0));

        // Should have weight updates
        let updates = stdp.get_weight_updates();
        assert!(!updates.is_empty());
        assert!(updates[0].delta > 0.0);
    }

    #[test]
    fn test_punishment() {
        let mut stdp = RewardModulatedStdp::new(10, RewardModulatedParams::default());
        stdp.set_synapse_mapping(0, 0);

        // Create positive eligibility
        stdp.on_pre_spike(0, 0.0);
        stdp.on_post_spike(0, 10.0);

        // Deliver punishment
        stdp.deliver_reward(RewardSignal::punishment(1.0, 50.0));

        // Should have negative weight update
        let updates = stdp.get_weight_updates();
        assert!(!updates.is_empty());
        assert!(updates[0].delta < 0.0);
    }
}
