//! Mutation Operators
//!
//! Operators that modify genomes to explore the search space.

#[cfg(feature = "serde")]
use serde::{Deserialize, Serialize};

use rand::Rng;
use rand_distr::{Distribution, Normal, Uniform};

use crate::genome::{ActivationType, ConnectionGene, Genome, InnovationTracker, NodeGene, NodeType};

/// Mutation parameters
#[derive(Debug, Clone)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct MutationParams {
    /// Probability of mutating a weight
    pub weight_mutation_rate: f32,

    /// Probability of replacing weight entirely
    pub weight_replacement_rate: f32,

    /// Standard deviation for weight perturbation
    pub weight_perturbation: f32,

    /// Probability of adding a new node
    pub add_node_rate: f32,

    /// Probability of adding a new connection
    pub add_connection_rate: f32,

    /// Probability of toggling connection enabled/disabled
    pub toggle_connection_rate: f32,

    /// Probability of mutating node bias
    pub bias_mutation_rate: f32,

    /// Standard deviation for bias perturbation
    pub bias_perturbation: f32,

    /// Probability of mutating activation function
    pub activation_mutation_rate: f32,

    /// Allow recurrent connections
    pub allow_recurrent: bool,

    /// Maximum number of attempts for adding connections
    pub max_connection_attempts: usize,

    /// Weight initialization range
    pub weight_init_range: f32,
}

impl Default for MutationParams {
    fn default() -> Self {
        Self {
            weight_mutation_rate: 0.8,
            weight_replacement_rate: 0.1,
            weight_perturbation: 0.5,
            add_node_rate: 0.03,
            add_connection_rate: 0.05,
            toggle_connection_rate: 0.01,
            bias_mutation_rate: 0.4,
            bias_perturbation: 0.3,
            activation_mutation_rate: 0.1,
            allow_recurrent: true,
            max_connection_attempts: 20,
            weight_init_range: 2.0,
        }
    }
}

impl MutationParams {
    /// Exploratory mutations (more structural changes)
    pub fn exploratory() -> Self {
        Self {
            weight_mutation_rate: 0.6,
            add_node_rate: 0.1,
            add_connection_rate: 0.15,
            toggle_connection_rate: 0.05,
            activation_mutation_rate: 0.2,
            ..Default::default()
        }
    }

    /// Conservative mutations (mostly weight changes)
    pub fn conservative() -> Self {
        Self {
            weight_mutation_rate: 0.9,
            add_node_rate: 0.01,
            add_connection_rate: 0.02,
            toggle_connection_rate: 0.005,
            activation_mutation_rate: 0.05,
            ..Default::default()
        }
    }

    /// Balanced mutations
    pub fn balanced() -> Self {
        Self::default()
    }
}

/// Mutation operator
pub struct MutationOperator {
    params: MutationParams,
}

impl MutationOperator {
    pub fn new(params: MutationParams) -> Self {
        Self { params }
    }

    /// Apply all mutations to a genome
    pub fn mutate<R: Rng>(
        &self,
        genome: &mut Genome,
        tracker: &mut InnovationTracker,
        rng: &mut R,
    ) {
        // Structural mutations first
        if rng.gen::<f32>() < self.params.add_node_rate {
            self.mutate_add_node(genome, tracker, rng);
        }

        if rng.gen::<f32>() < self.params.add_connection_rate {
            self.mutate_add_connection(genome, tracker, rng);
        }

        // Weight mutations
        self.mutate_weights(genome, rng);

        // Bias mutations
        self.mutate_biases(genome, rng);

        // Toggle connections
        if rng.gen::<f32>() < self.params.toggle_connection_rate {
            self.mutate_toggle_connection(genome, rng);
        }

        // Activation function mutations
        if rng.gen::<f32>() < self.params.activation_mutation_rate {
            self.mutate_activation(genome, rng);
        }
    }

    /// Mutate weights
    pub fn mutate_weights<R: Rng>(&self, genome: &mut Genome, rng: &mut R) {
        let normal = Normal::new(0.0, self.params.weight_perturbation as f64).unwrap();
        let uniform = Uniform::new(-self.params.weight_init_range, self.params.weight_init_range);

        for conn in &mut genome.connections {
            if rng.gen::<f32>() < self.params.weight_mutation_rate {
                if rng.gen::<f32>() < self.params.weight_replacement_rate {
                    // Replace weight entirely
                    conn.weight = uniform.sample(rng);
                } else {
                    // Perturb weight
                    conn.weight += normal.sample(rng) as f32;
                }
            }
        }
    }

    /// Mutate biases
    pub fn mutate_biases<R: Rng>(&self, genome: &mut Genome, rng: &mut R) {
        let normal = Normal::new(0.0, self.params.bias_perturbation as f64).unwrap();

        for node in &mut genome.nodes {
            if node.node_type == NodeType::Hidden || node.node_type == NodeType::Output {
                if rng.gen::<f32>() < self.params.bias_mutation_rate {
                    node.bias += normal.sample(rng) as f32;
                }
            }
        }
    }

    /// Add a new node by splitting a connection
    pub fn mutate_add_node<R: Rng>(
        &self,
        genome: &mut Genome,
        tracker: &mut InnovationTracker,
        rng: &mut R,
    ) {
        // Find enabled connections to split
        let enabled: Vec<usize> = genome
            .connections
            .iter()
            .enumerate()
            .filter(|(_, c)| c.enabled)
            .map(|(i, _)| i)
            .collect();

        if enabled.is_empty() {
            return;
        }

        // Pick random connection to split
        let conn_idx = enabled[rng.gen_range(0..enabled.len())];
        let conn = &genome.connections[conn_idx];

        let from = conn.from;
        let to = conn.to;
        let old_weight = conn.weight;
        let old_innovation = conn.innovation;

        // Disable old connection
        genome.connections[conn_idx].enabled = false;

        // Create new node
        let new_node_id = tracker.get_split_node(old_innovation);

        // Determine layer for new node
        let from_layer = genome.get_node(from).map(|n| n.layer).unwrap_or(0);
        let to_layer = genome.get_node(to).map(|n| n.layer).unwrap_or(u32::MAX);
        let new_layer = if to_layer > from_layer + 1 {
            from_layer + 1
        } else {
            from_layer + 1
        };

        let new_node = NodeGene::hidden(new_node_id, new_layer);
        genome.nodes.push(new_node);

        // Create two new connections
        // Connection 1: from -> new (weight 1.0)
        let inno1 = tracker.get_connection_innovation(from, new_node_id);
        genome.connections.push(ConnectionGene::new(inno1, from, new_node_id, 1.0));

        // Connection 2: new -> to (old weight)
        let inno2 = tracker.get_connection_innovation(new_node_id, to);
        genome.connections.push(ConnectionGene::new(inno2, new_node_id, to, old_weight));
    }

    /// Add a new connection
    pub fn mutate_add_connection<R: Rng>(
        &self,
        genome: &mut Genome,
        tracker: &mut InnovationTracker,
        rng: &mut R,
    ) {
        let uniform = Uniform::new(-self.params.weight_init_range, self.params.weight_init_range);

        // Get valid source nodes (not output)
        let sources: Vec<u32> = genome
            .nodes
            .iter()
            .filter(|n| !n.node_type.is_output())
            .map(|n| n.id)
            .collect();

        // Get valid target nodes (not input/bias)
        let targets: Vec<u32> = genome
            .nodes
            .iter()
            .filter(|n| n.node_type.can_receive())
            .map(|n| n.id)
            .collect();

        if sources.is_empty() || targets.is_empty() {
            return;
        }

        // Try to find a valid connection
        for _ in 0..self.params.max_connection_attempts {
            let from = sources[rng.gen_range(0..sources.len())];
            let to = targets[rng.gen_range(0..targets.len())];

            // Skip if connection already exists
            if genome.has_connection(from, to) {
                continue;
            }

            // Check for cycles in feedforward mode
            let is_recurrent = if !self.params.allow_recurrent {
                if genome.would_create_cycle(from, to) {
                    continue;
                }
                false
            } else {
                genome.would_create_cycle(from, to)
            };

            // Create the connection
            let inno = tracker.get_connection_innovation(from, to);
            let weight = uniform.sample(rng);

            let mut conn = ConnectionGene::new(inno, from, to, weight);
            conn.recurrent = is_recurrent;
            genome.connections.push(conn);

            return;
        }
    }

    /// Toggle a connection's enabled state
    pub fn mutate_toggle_connection<R: Rng>(&self, genome: &mut Genome, rng: &mut R) {
        if genome.connections.is_empty() {
            return;
        }

        let idx = rng.gen_range(0..genome.connections.len());
        genome.connections[idx].enabled = !genome.connections[idx].enabled;
    }

    /// Mutate activation function of a hidden node
    pub fn mutate_activation<R: Rng>(&self, genome: &mut Genome, rng: &mut R) {
        let hidden: Vec<usize> = genome
            .nodes
            .iter()
            .enumerate()
            .filter(|(_, n)| n.node_type == NodeType::Hidden)
            .map(|(i, _)| i)
            .collect();

        if hidden.is_empty() {
            return;
        }

        let idx = hidden[rng.gen_range(0..hidden.len())];
        let activations = ActivationType::all();
        genome.nodes[idx].activation = activations[rng.gen_range(0..activations.len())];
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_weight_mutation() {
        let mut rng = rand::thread_rng();
        let mut tracker = InnovationTracker::new();
        let mut genome = Genome::minimal(2, 1, &mut tracker);

        let old_weights: Vec<f32> = genome.connections.iter().map(|c| c.weight).collect();

        let mutator = MutationOperator::new(MutationParams {
            weight_mutation_rate: 1.0,
            ..Default::default()
        });

        mutator.mutate_weights(&mut genome, &mut rng);

        // At least some weights should have changed
        let changed = genome
            .connections
            .iter()
            .zip(old_weights.iter())
            .any(|(c, &old)| (c.weight - old).abs() > 1e-6);

        assert!(changed);
    }

    #[test]
    fn test_add_node() {
        let mut rng = rand::thread_rng();
        let mut tracker = InnovationTracker::new();
        let mut genome = Genome::minimal(2, 1, &mut tracker);

        // Enable some connections
        for conn in &mut genome.connections {
            conn.enabled = true;
            conn.weight = 1.0;
        }

        let initial_nodes = genome.nodes.len();
        let initial_conns = genome.connections.len();

        let mutator = MutationOperator::new(MutationParams::default());
        mutator.mutate_add_node(&mut genome, &mut tracker, &mut rng);

        assert_eq!(genome.nodes.len(), initial_nodes + 1);
        assert_eq!(genome.connections.len(), initial_conns + 2);
    }
}
