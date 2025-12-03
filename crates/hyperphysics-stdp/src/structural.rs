//! Structural Plasticity
//!
//! Synapse creation and pruning based on activity patterns.
//! Enables network topology adaptation and self-organization.

#[cfg(feature = "serde")]
use serde::{Deserialize, Serialize};

use rand::Rng;
use rand_distr::{Distribution, Uniform};

use crate::weight_bounds::WeightBounds;
use crate::{PlasticityRule, PlasticityStats, WeightUpdate};

/// Structural plasticity parameters
#[derive(Debug, Clone)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct StructuralParams {
    /// Minimum weight before synapse is pruned
    pub prune_threshold: f32,

    /// Probability of creating a new synapse per timestep
    pub creation_rate: f32,

    /// Initial weight for new synapses
    pub initial_weight: f32,

    /// Maximum number of synapses per neuron
    pub max_synapses_per_neuron: usize,

    /// Time between structural updates (ms)
    pub update_interval: f32,

    /// Enable activity-dependent creation
    pub activity_dependent: bool,

    /// Weight bounds
    pub bounds: WeightBounds,
}

impl Default for StructuralParams {
    fn default() -> Self {
        Self {
            prune_threshold: 0.01,
            creation_rate: 0.001,
            initial_weight: 0.5,
            max_synapses_per_neuron: 100,
            update_interval: 1000.0, // 1 second
            activity_dependent: true,
            bounds: WeightBounds::default(),
        }
    }
}

impl StructuralParams {
    /// Conservative (less structural change)
    pub fn conservative() -> Self {
        Self {
            prune_threshold: 0.001,
            creation_rate: 0.0001,
            update_interval: 10000.0,
            ..Default::default()
        }
    }

    /// Aggressive (more structural change)
    pub fn aggressive() -> Self {
        Self {
            prune_threshold: 0.1,
            creation_rate: 0.01,
            update_interval: 100.0,
            ..Default::default()
        }
    }
}

/// Candidate for new synapse creation
#[derive(Debug, Clone, Copy)]
pub struct SynapseCandidate {
    /// Presynaptic neuron
    pub pre: usize,

    /// Postsynaptic neuron
    pub post: usize,

    /// Proposed initial weight
    pub weight: f32,

    /// Priority/score for creation
    pub priority: f32,
}

/// Structural plasticity events
#[derive(Debug, Clone)]
pub enum StructuralEvent {
    /// Synapse created
    Created {
        synapse_id: usize,
        pre: usize,
        post: usize,
        weight: f32,
    },

    /// Synapse pruned
    Pruned {
        synapse_id: usize,
        pre: usize,
        post: usize,
    },
}

/// Structural plasticity implementation
pub struct StructuralPlasticity {
    /// Parameters
    params: StructuralParams,

    /// Number of neurons
    num_neurons: usize,

    /// Synapse to pre/post mapping
    synapse_endpoints: Vec<(usize, usize)>, // (pre, post)

    /// Synapse count per postsynaptic neuron
    post_synapse_count: Vec<usize>,

    /// Activity trace per neuron (for activity-dependent creation)
    activity_trace: Vec<f32>,

    /// Pending structural events
    pending_events: Vec<StructuralEvent>,

    /// Last update time
    last_update: f64,

    /// Statistics
    synapses_created: u64,
    synapses_pruned: u64,
}

impl StructuralPlasticity {
    /// Create new structural plasticity rule
    pub fn new(num_neurons: usize, num_synapses: usize, params: StructuralParams) -> Self {
        // Default synapse endpoints
        let synapse_endpoints: Vec<(usize, usize)> = (0..num_synapses)
            .map(|i| {
                let pre = i / num_neurons.max(1);
                let post = i % num_neurons;
                (pre, post)
            })
            .collect();

        // Count synapses per post neuron
        let mut post_synapse_count = vec![0; num_neurons];
        for &(_, post) in &synapse_endpoints {
            if post < num_neurons {
                post_synapse_count[post] += 1;
            }
        }

        Self {
            params,
            num_neurons,
            synapse_endpoints,
            post_synapse_count,
            activity_trace: vec![0.0; num_neurons],
            pending_events: Vec::new(),
            last_update: 0.0,
            synapses_created: 0,
            synapses_pruned: 0,
        }
    }

    /// Check weights and prune weak synapses
    pub fn check_and_prune(&mut self, weights: &[f32]) -> Vec<usize> {
        let mut pruned = Vec::new();

        for (syn_id, &weight) in weights.iter().enumerate() {
            if weight.abs() < self.params.prune_threshold {
                if syn_id < self.synapse_endpoints.len() {
                    let (pre, post) = self.synapse_endpoints[syn_id];
                    self.pending_events.push(StructuralEvent::Pruned {
                        synapse_id: syn_id,
                        pre,
                        post,
                    });
                    pruned.push(syn_id);

                    if post < self.post_synapse_count.len() {
                        self.post_synapse_count[post] = self.post_synapse_count[post].saturating_sub(1);
                    }
                    self.synapses_pruned += 1;
                }
            }
        }

        pruned
    }

    /// Generate candidates for new synapses
    pub fn generate_candidates<R: Rng>(&self, rng: &mut R, count: usize) -> Vec<SynapseCandidate> {
        let mut candidates = Vec::new();

        if self.num_neurons < 2 {
            return candidates;
        }

        let neuron_dist = Uniform::new(0, self.num_neurons);

        for _ in 0..count {
            let pre = neuron_dist.sample(rng);
            let mut post = neuron_dist.sample(rng);

            // Avoid self-connections
            while post == pre {
                post = neuron_dist.sample(rng);
            }

            // Check if post can accept more synapses
            if post < self.post_synapse_count.len()
                && self.post_synapse_count[post] >= self.params.max_synapses_per_neuron
            {
                continue;
            }

            let priority = if self.params.activity_dependent {
                // Prefer connecting active neurons
                let pre_activity = self.activity_trace.get(pre).copied().unwrap_or(0.0);
                let post_activity = self.activity_trace.get(post).copied().unwrap_or(0.0);
                pre_activity * post_activity
            } else {
                rng.gen::<f32>()
            };

            candidates.push(SynapseCandidate {
                pre,
                post,
                weight: self.params.initial_weight,
                priority,
            });
        }

        // Sort by priority
        candidates.sort_by(|a, b| b.priority.partial_cmp(&a.priority).unwrap_or(std::cmp::Ordering::Equal));

        candidates
    }

    /// Accept a synapse candidate (create synapse)
    pub fn accept_candidate(&mut self, candidate: &SynapseCandidate, synapse_id: usize) {
        self.pending_events.push(StructuralEvent::Created {
            synapse_id,
            pre: candidate.pre,
            post: candidate.post,
            weight: candidate.weight,
        });

        if candidate.post < self.post_synapse_count.len() {
            self.post_synapse_count[candidate.post] += 1;
        }

        // Update endpoint mapping
        if synapse_id < self.synapse_endpoints.len() {
            self.synapse_endpoints[synapse_id] = (candidate.pre, candidate.post);
        } else {
            self.synapse_endpoints.push((candidate.pre, candidate.post));
        }

        self.synapses_created += 1;
    }

    /// Get pending structural events
    pub fn take_events(&mut self) -> Vec<StructuralEvent> {
        std::mem::take(&mut self.pending_events)
    }

    /// Update activity trace
    fn update_activity(&mut self, neuron_id: usize, time: f64) {
        if neuron_id < self.activity_trace.len() {
            // Simple exponential trace
            let dt = (time - self.last_update).max(0.0) as f32;
            let decay = (-dt / 1000.0).exp(); // 1 second time constant

            self.activity_trace[neuron_id] =
                self.activity_trace[neuron_id] * decay + 1.0;
        }
    }

    /// Set synapse endpoints
    pub fn set_synapse_endpoints(&mut self, synapse_id: usize, pre: usize, post: usize) {
        if synapse_id < self.synapse_endpoints.len() {
            let old_post = self.synapse_endpoints[synapse_id].1;
            if old_post < self.post_synapse_count.len() {
                self.post_synapse_count[old_post] = self.post_synapse_count[old_post].saturating_sub(1);
            }

            self.synapse_endpoints[synapse_id] = (pre, post);

            if post < self.post_synapse_count.len() {
                self.post_synapse_count[post] += 1;
            }
        }
    }
}

impl PlasticityRule for StructuralPlasticity {
    fn on_pre_spike(&mut self, synapse_id: usize, time: f64) {
        if synapse_id < self.synapse_endpoints.len() {
            let (pre, _) = self.synapse_endpoints[synapse_id];
            self.update_activity(pre, time);
        }
    }

    fn on_post_spike(&mut self, neuron_id: usize, time: f64) {
        self.update_activity(neuron_id, time);
        self.last_update = time;
    }

    fn get_weight_updates(&self) -> Vec<WeightUpdate> {
        // Structural plasticity doesn't produce weight updates directly
        Vec::new()
    }

    fn apply_updates(&mut self, _weights: &mut [f32]) {
        // Weight updates are handled through events
    }

    fn reset(&mut self) {
        self.activity_trace.fill(0.0);
        self.pending_events.clear();
        self.last_update = 0.0;
    }

    fn stats(&self) -> PlasticityStats {
        PlasticityStats {
            ltp_count: self.synapses_created,
            ltd_count: self.synapses_pruned,
            avg_weight_change: 0.0,
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
    fn test_structural_creation() {
        let sp = StructuralPlasticity::new(10, 100, StructuralParams::default());
        assert_eq!(sp.num_neurons, 10);
    }

    #[test]
    fn test_pruning() {
        let mut sp = StructuralPlasticity::new(10, 10, StructuralParams::default());
        sp.synapse_endpoints = vec![(0, 1); 10];

        let weights = vec![0.001; 10]; // All below threshold
        let pruned = sp.check_and_prune(&weights);

        assert_eq!(pruned.len(), 10);
        assert_eq!(sp.synapses_pruned, 10);
    }

    #[test]
    fn test_candidate_generation() {
        let sp = StructuralPlasticity::new(10, 0, StructuralParams::default());
        let mut rng = rand::thread_rng();

        let candidates = sp.generate_candidates(&mut rng, 5);

        assert_eq!(candidates.len(), 5);
        for c in &candidates {
            assert_ne!(c.pre, c.post); // No self-connections
        }
    }
}
