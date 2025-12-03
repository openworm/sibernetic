//! HyperPhysics STDP - Spike-Timing Dependent Plasticity
//!
//! Implements biologically-inspired learning rules for spiking neural networks.
//!
//! # Overview
//!
//! STDP (Spike-Timing Dependent Plasticity) modifies synaptic weights based on
//! the relative timing of pre- and post-synaptic spikes:
//!
//! - Pre before post (causal) → Long-Term Potentiation (LTP) → weight increase
//! - Post before pre (anti-causal) → Long-Term Depression (LTD) → weight decrease
//!
//! # Learning Rules Implemented
//!
//! - **Classical STDP**: Symmetric exponential windows
//! - **Triplet STDP**: Three-factor rule for more realistic dynamics
//! - **Reward-Modulated STDP**: Eligibility traces + reward signal
//! - **Homeostatic STDP**: Maintains target firing rates
//! - **Structural Plasticity**: Synapse creation/pruning

mod classical;
mod triplet;
mod reward;
mod homeostatic;
mod structural;
mod traces;
mod weight_bounds;

pub use classical::{ClassicalStdp, ClassicalStdpParams};
pub use triplet::{TripletStdp, TripletStdpParams};
pub use reward::{RewardModulatedStdp, RewardModulatedParams, RewardSignal};
pub use homeostatic::{HomeostaticStdp, HomeostaticParams};
pub use structural::{StructuralPlasticity, StructuralParams, SynapseCandidate};
pub use traces::{EligibilityTrace, SpikeTrace, TraceParams};
pub use weight_bounds::{WeightBounds, WeightUpdate};

use thiserror::Error;

/// STDP errors
#[derive(Debug, Error)]
pub enum StdpError {
    #[error("Invalid neuron ID: {id} (max: {max})")]
    InvalidNeuronId { id: usize, max: usize },

    #[error("Invalid synapse ID: {id} (max: {max})")]
    InvalidSynapseId { id: usize, max: usize },

    #[error("Weight out of bounds: {weight} (min: {min}, max: {max})")]
    WeightOutOfBounds { weight: f32, min: f32, max: f32 },

    #[error("Invalid learning rate: {rate}")]
    InvalidLearningRate { rate: f32 },
}

pub type Result<T> = std::result::Result<T, StdpError>;

/// Common interface for STDP rules
pub trait PlasticityRule: Send + Sync {
    /// Process a presynaptic spike
    fn on_pre_spike(&mut self, synapse_id: usize, time: f64);

    /// Process a postsynaptic spike
    fn on_post_spike(&mut self, neuron_id: usize, time: f64);

    /// Get weight updates for all synapses
    fn get_weight_updates(&self) -> Vec<WeightUpdate>;

    /// Apply weight updates
    fn apply_updates(&mut self, weights: &mut [f32]);

    /// Reset learning state
    fn reset(&mut self);

    /// Get learning statistics
    fn stats(&self) -> PlasticityStats;
}

/// Learning statistics
#[derive(Debug, Clone, Default)]
pub struct PlasticityStats {
    /// Total LTP events
    pub ltp_count: u64,

    /// Total LTD events
    pub ltd_count: u64,

    /// Average weight change
    pub avg_weight_change: f64,

    /// Maximum weight change
    pub max_weight_change: f32,

    /// Number of weights at upper bound
    pub at_upper_bound: usize,

    /// Number of weights at lower bound
    pub at_lower_bound: usize,
}

/// Composite plasticity controller
/// Combines multiple plasticity rules
pub struct PlasticityController {
    /// Active plasticity rules
    rules: Vec<Box<dyn PlasticityRule>>,

    /// Global learning rate multiplier
    learning_rate: f32,

    /// Enable/disable learning
    enabled: bool,

    /// Weight bounds
    bounds: WeightBounds,
}

impl Default for PlasticityController {
    fn default() -> Self {
        Self {
            rules: Vec::new(),
            learning_rate: 1.0,
            enabled: true,
            bounds: WeightBounds::default(),
        }
    }
}

impl PlasticityController {
    /// Create new plasticity controller
    pub fn new() -> Self {
        Self::default()
    }

    /// Add a plasticity rule
    pub fn add_rule(&mut self, rule: Box<dyn PlasticityRule>) {
        self.rules.push(rule);
    }

    /// Create with classical STDP
    pub fn with_classical_stdp(num_synapses: usize) -> Self {
        let mut controller = Self::new();
        controller.add_rule(Box::new(ClassicalStdp::new(
            num_synapses,
            ClassicalStdpParams::default(),
        )));
        controller
    }

    /// Create with reward-modulated STDP
    pub fn with_reward_stdp(num_synapses: usize) -> Self {
        let mut controller = Self::new();
        controller.add_rule(Box::new(RewardModulatedStdp::new(
            num_synapses,
            RewardModulatedParams::default(),
        )));
        controller
    }

    /// Set global learning rate
    pub fn set_learning_rate(&mut self, rate: f32) {
        self.learning_rate = rate;
    }

    /// Enable/disable learning
    pub fn set_enabled(&mut self, enabled: bool) {
        self.enabled = enabled;
    }

    /// Process presynaptic spike
    pub fn on_pre_spike(&mut self, synapse_id: usize, time: f64) {
        if !self.enabled {
            return;
        }
        for rule in &mut self.rules {
            rule.on_pre_spike(synapse_id, time);
        }
    }

    /// Process postsynaptic spike
    pub fn on_post_spike(&mut self, neuron_id: usize, time: f64) {
        if !self.enabled {
            return;
        }
        for rule in &mut self.rules {
            rule.on_post_spike(neuron_id, time);
        }
    }

    /// Apply all weight updates
    pub fn apply(&mut self, weights: &mut [f32]) {
        if !self.enabled {
            return;
        }

        for rule in &mut self.rules {
            let updates = rule.get_weight_updates();
            for update in updates {
                if update.synapse_id < weights.len() {
                    let delta = update.delta * self.learning_rate;
                    weights[update.synapse_id] =
                        self.bounds.clamp(weights[update.synapse_id] + delta);
                }
            }
        }
    }

    /// Reset all rules
    pub fn reset(&mut self) {
        for rule in &mut self.rules {
            rule.reset();
        }
    }

    /// Get combined statistics
    pub fn stats(&self) -> PlasticityStats {
        let mut combined = PlasticityStats::default();
        for rule in &self.rules {
            let s = rule.stats();
            combined.ltp_count += s.ltp_count;
            combined.ltd_count += s.ltd_count;
            combined.at_upper_bound += s.at_upper_bound;
            combined.at_lower_bound += s.at_lower_bound;
            if s.max_weight_change > combined.max_weight_change {
                combined.max_weight_change = s.max_weight_change;
            }
        }
        combined
    }
}
