//! Spike Traces and Eligibility Traces
//!
//! Traces are low-pass filtered versions of spike trains used in various
//! plasticity rules.

#[cfg(feature = "serde")]
use serde::{Deserialize, Serialize};

/// Parameters for spike traces
#[derive(Debug, Clone, Copy)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct TraceParams {
    /// Time constant for trace decay (ms)
    pub tau: f32,

    /// Increment on each spike
    pub increment: f32,

    /// Maximum trace value
    pub max_value: f32,
}

impl Default for TraceParams {
    fn default() -> Self {
        Self {
            tau: 20.0,
            increment: 1.0,
            max_value: 10.0,
        }
    }
}

impl TraceParams {
    /// Fast trace (for precise timing)
    pub fn fast() -> Self {
        Self {
            tau: 10.0,
            increment: 1.0,
            max_value: 10.0,
        }
    }

    /// Slow trace (for rate-based learning)
    pub fn slow() -> Self {
        Self {
            tau: 100.0,
            increment: 1.0,
            max_value: 10.0,
        }
    }

    /// Eligibility trace parameters
    pub fn eligibility() -> Self {
        Self {
            tau: 1000.0, // 1 second eligibility window
            increment: 0.1,
            max_value: 1.0,
        }
    }
}

/// A single exponentially-decaying spike trace
#[derive(Debug, Clone, Copy, Default)]
pub struct SpikeTrace {
    /// Current trace value
    pub value: f32,

    /// Last update time (ms)
    last_time: f64,
}

impl SpikeTrace {
    /// Create new trace
    pub fn new() -> Self {
        Self::default()
    }

    /// Update trace to current time (apply decay)
    pub fn update(&mut self, time: f64, params: &TraceParams) {
        if time > self.last_time {
            let dt = (time - self.last_time) as f32;
            self.value *= (-dt / params.tau).exp();
            self.last_time = time;
        }
    }

    /// Increment trace on spike
    pub fn on_spike(&mut self, time: f64, params: &TraceParams) {
        self.update(time, params);
        self.value = (self.value + params.increment).min(params.max_value);
    }

    /// Get decayed value at given time
    pub fn value_at(&self, time: f64, params: &TraceParams) -> f32 {
        if time <= self.last_time {
            return self.value;
        }
        let dt = (time - self.last_time) as f32;
        self.value * (-dt / params.tau).exp()
    }

    /// Reset trace
    pub fn reset(&mut self) {
        self.value = 0.0;
        self.last_time = 0.0;
    }
}

/// Eligibility trace for reward-modulated learning
#[derive(Debug, Clone)]
pub struct EligibilityTrace {
    /// Per-synapse eligibility
    traces: Vec<f32>,

    /// Parameters
    params: TraceParams,

    /// Last update time
    last_time: f64,
}

impl EligibilityTrace {
    /// Create new eligibility trace buffer
    pub fn new(num_synapses: usize, params: TraceParams) -> Self {
        Self {
            traces: vec![0.0; num_synapses],
            params,
            last_time: 0.0,
        }
    }

    /// Decay all traces to current time
    pub fn decay(&mut self, time: f64) {
        if time <= self.last_time {
            return;
        }

        let dt = (time - self.last_time) as f32;
        let decay = (-dt / self.params.tau).exp();

        for trace in &mut self.traces {
            *trace *= decay;
        }

        self.last_time = time;
    }

    /// Add eligibility for a synapse
    pub fn add(&mut self, synapse_id: usize, amount: f32, time: f64) {
        self.decay(time);
        if synapse_id < self.traces.len() {
            self.traces[synapse_id] = (self.traces[synapse_id] + amount).min(self.params.max_value);
        }
    }

    /// Get eligibility for a synapse
    pub fn get(&self, synapse_id: usize) -> f32 {
        self.traces.get(synapse_id).copied().unwrap_or(0.0)
    }

    /// Get all eligibility values
    pub fn all(&self) -> &[f32] {
        &self.traces
    }

    /// Apply reward to all eligible synapses
    pub fn apply_reward(&self, reward: f32, learning_rate: f32) -> Vec<(usize, f32)> {
        self.traces
            .iter()
            .enumerate()
            .filter(|(_, &e)| e.abs() > 1e-6)
            .map(|(i, &e)| (i, e * reward * learning_rate))
            .collect()
    }

    /// Reset all traces
    pub fn reset(&mut self) {
        self.traces.fill(0.0);
        self.last_time = 0.0;
    }

    /// Get number of synapses
    pub fn len(&self) -> usize {
        self.traces.len()
    }

    /// Check if empty
    pub fn is_empty(&self) -> bool {
        self.traces.is_empty()
    }
}

/// Pre-post spike pair detector
/// Detects causal and anti-causal spike pairs for STDP
#[derive(Debug, Clone)]
pub struct PairDetector {
    /// Pre-synaptic traces (one per synapse)
    pre_traces: Vec<SpikeTrace>,

    /// Post-synaptic traces (one per neuron)
    post_traces: Vec<SpikeTrace>,

    /// Trace parameters
    params: TraceParams,
}

impl PairDetector {
    /// Create new pair detector
    pub fn new(num_synapses: usize, num_neurons: usize, params: TraceParams) -> Self {
        Self {
            pre_traces: vec![SpikeTrace::new(); num_synapses],
            post_traces: vec![SpikeTrace::new(); num_neurons],
            params,
        }
    }

    /// Record pre-synaptic spike, return post trace value (for LTP)
    pub fn on_pre_spike(&mut self, synapse_id: usize, post_neuron: usize, time: f64) -> f32 {
        // Get post trace before updating pre
        let post_value = if post_neuron < self.post_traces.len() {
            self.post_traces[post_neuron].value_at(time, &self.params)
        } else {
            0.0
        };

        // Update pre trace
        if synapse_id < self.pre_traces.len() {
            self.pre_traces[synapse_id].on_spike(time, &self.params);
        }

        post_value // Positive if post fired recently (anti-causal → LTD)
    }

    /// Record post-synaptic spike, return pre trace value (for LTD)
    pub fn on_post_spike(&mut self, neuron_id: usize, synapse_ids: &[usize], time: f64) -> Vec<f32> {
        // Get pre traces before updating post
        let pre_values: Vec<f32> = synapse_ids
            .iter()
            .map(|&sid| {
                if sid < self.pre_traces.len() {
                    self.pre_traces[sid].value_at(time, &self.params)
                } else {
                    0.0
                }
            })
            .collect();

        // Update post trace
        if neuron_id < self.post_traces.len() {
            self.post_traces[neuron_id].on_spike(time, &self.params);
        }

        pre_values // Positive if pre fired recently (causal → LTP)
    }

    /// Reset all traces
    pub fn reset(&mut self) {
        for trace in &mut self.pre_traces {
            trace.reset();
        }
        for trace in &mut self.post_traces {
            trace.reset();
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_spike_trace() {
        let mut trace = SpikeTrace::new();
        let params = TraceParams::default();

        trace.on_spike(0.0, &params);
        assert!((trace.value - 1.0).abs() < 1e-6);

        // After one time constant, should be ~37% of original
        let decayed = trace.value_at(20.0, &params);
        assert!((decayed - 0.368).abs() < 0.01);
    }

    #[test]
    fn test_eligibility_trace() {
        let params = TraceParams::eligibility();
        let mut trace = EligibilityTrace::new(10, params);

        trace.add(5, 0.5, 0.0);
        assert!((trace.get(5) - 0.5).abs() < 1e-6);

        trace.decay(100.0);
        assert!(trace.get(5) < 0.5); // Should have decayed
    }

    #[test]
    fn test_pair_detector() {
        let params = TraceParams::default();
        let mut detector = PairDetector::new(10, 5, params);

        // Pre then post (causal → LTP)
        detector.on_pre_spike(0, 0, 0.0);
        let pre_values = detector.on_post_spike(0, &[0], 10.0);
        assert!(pre_values[0] > 0.5); // Pre was recent

        // Post then pre (anti-causal → LTD)
        let mut detector2 = PairDetector::new(10, 5, params);
        detector2.on_post_spike(1, &[], 0.0);
        let post_value = detector2.on_pre_spike(1, 1, 10.0);
        assert!(post_value > 0.5); // Post was recent
    }
}
