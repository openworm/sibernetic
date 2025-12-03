//! Genome Encoding
//!
//! Represents neural network topology as a genome that can be
//! mutated, crossed over, and decoded into a functional network.

#[cfg(feature = "serde")]
use serde::{Deserialize, Serialize};

use std::collections::HashMap;
use hashbrown::HashSet;

/// Node type in the network
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub enum NodeType {
    /// Input node (sensory)
    Input,
    /// Hidden node (interneuron)
    Hidden,
    /// Output node (motor)
    Output,
    /// Bias node
    Bias,
}

impl NodeType {
    pub fn is_input(&self) -> bool {
        matches!(self, Self::Input | Self::Bias)
    }

    pub fn is_output(&self) -> bool {
        matches!(self, Self::Output)
    }

    pub fn can_receive(&self) -> bool {
        !matches!(self, Self::Input | Self::Bias)
    }

    pub fn can_send(&self) -> bool {
        !matches!(self, Self::Output)
    }
}

/// A node gene in the genome
#[derive(Debug, Clone)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct NodeGene {
    /// Unique node ID
    pub id: u32,

    /// Node type
    pub node_type: NodeType,

    /// Activation function type
    pub activation: ActivationType,

    /// Bias value
    pub bias: f32,

    /// Layer hint (for feedforward networks)
    pub layer: u32,
}

impl NodeGene {
    pub fn input(id: u32) -> Self {
        Self {
            id,
            node_type: NodeType::Input,
            activation: ActivationType::Identity,
            bias: 0.0,
            layer: 0,
        }
    }

    pub fn hidden(id: u32, layer: u32) -> Self {
        Self {
            id,
            node_type: NodeType::Hidden,
            activation: ActivationType::Tanh,
            bias: 0.0,
            layer,
        }
    }

    pub fn output(id: u32) -> Self {
        Self {
            id,
            node_type: NodeType::Output,
            activation: ActivationType::Tanh,
            bias: 0.0,
            layer: u32::MAX,
        }
    }

    pub fn bias(id: u32) -> Self {
        Self {
            id,
            node_type: NodeType::Bias,
            activation: ActivationType::Identity,
            bias: 1.0,
            layer: 0,
        }
    }
}

/// Activation function types
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub enum ActivationType {
    Identity,
    Sigmoid,
    Tanh,
    ReLU,
    LeakyReLU,
    Gaussian,
    Sin,
    Step,
}

impl ActivationType {
    pub fn apply(&self, x: f32) -> f32 {
        match self {
            Self::Identity => x,
            Self::Sigmoid => 1.0 / (1.0 + (-x).exp()),
            Self::Tanh => x.tanh(),
            Self::ReLU => x.max(0.0),
            Self::LeakyReLU => if x > 0.0 { x } else { 0.01 * x },
            Self::Gaussian => (-x * x).exp(),
            Self::Sin => x.sin(),
            Self::Step => if x > 0.0 { 1.0 } else { 0.0 },
        }
    }

    /// Get all activation types for random selection
    pub fn all() -> &'static [ActivationType] {
        &[
            Self::Identity,
            Self::Sigmoid,
            Self::Tanh,
            Self::ReLU,
            Self::LeakyReLU,
            Self::Gaussian,
            Self::Sin,
            Self::Step,
        ]
    }
}

/// A connection gene in the genome
#[derive(Debug, Clone)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct ConnectionGene {
    /// Innovation number (for crossover alignment)
    pub innovation: u32,

    /// Source node ID
    pub from: u32,

    /// Target node ID
    pub to: u32,

    /// Connection weight
    pub weight: f32,

    /// Is this connection enabled?
    pub enabled: bool,

    /// Is this a recurrent connection?
    pub recurrent: bool,
}

impl ConnectionGene {
    pub fn new(innovation: u32, from: u32, to: u32, weight: f32) -> Self {
        Self {
            innovation,
            from,
            to,
            weight,
            enabled: true,
            recurrent: false,
        }
    }

    pub fn recurrent(innovation: u32, from: u32, to: u32, weight: f32) -> Self {
        Self {
            innovation,
            from,
            to,
            weight,
            enabled: true,
            recurrent: true,
        }
    }
}

/// Tracks innovation numbers for new genes
#[derive(Debug, Clone, Default)]
pub struct InnovationTracker {
    /// Current innovation number
    current: u32,

    /// Map of (from, to) -> innovation number
    connection_innovations: HashMap<(u32, u32), u32>,

    /// Current node ID
    current_node: u32,

    /// Map of split connection -> new node ID
    split_innovations: HashMap<u32, u32>,
}

impl InnovationTracker {
    pub fn new() -> Self {
        Self::default()
    }

    /// Initialize with existing nodes
    pub fn with_nodes(num_nodes: u32) -> Self {
        Self {
            current: 0,
            connection_innovations: HashMap::new(),
            current_node: num_nodes,
            split_innovations: HashMap::new(),
        }
    }

    /// Get or create innovation number for a connection
    pub fn get_connection_innovation(&mut self, from: u32, to: u32) -> u32 {
        let key = (from, to);
        if let Some(&inno) = self.connection_innovations.get(&key) {
            inno
        } else {
            let inno = self.current;
            self.current += 1;
            self.connection_innovations.insert(key, inno);
            inno
        }
    }

    /// Get or create node ID for splitting a connection
    pub fn get_split_node(&mut self, connection_innovation: u32) -> u32 {
        if let Some(&node) = self.split_innovations.get(&connection_innovation) {
            node
        } else {
            let node = self.current_node;
            self.current_node += 1;
            self.split_innovations.insert(connection_innovation, node);
            node
        }
    }

    /// Get current innovation count
    pub fn innovation_count(&self) -> u32 {
        self.current
    }

    /// Get current node count
    pub fn node_count(&self) -> u32 {
        self.current_node
    }
}

/// A complete genome encoding a neural network
#[derive(Debug, Clone)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct Genome {
    /// Node genes
    pub nodes: Vec<NodeGene>,

    /// Connection genes (sorted by innovation number)
    pub connections: Vec<ConnectionGene>,

    /// Number of input nodes
    pub num_inputs: usize,

    /// Number of output nodes
    pub num_outputs: usize,

    /// Fitness value (set during evaluation)
    pub fitness: f32,

    /// Adjusted fitness (after speciation)
    pub adjusted_fitness: f32,

    /// Species ID
    pub species_id: usize,
}

impl Genome {
    /// Create minimal genome with only inputs and outputs
    pub fn minimal(num_inputs: usize, num_outputs: usize, tracker: &mut InnovationTracker) -> Self {
        let mut nodes = Vec::with_capacity(num_inputs + num_outputs + 1);

        // Input nodes
        for i in 0..num_inputs {
            nodes.push(NodeGene::input(i as u32));
        }

        // Bias node
        nodes.push(NodeGene::bias(num_inputs as u32));

        // Output nodes
        for i in 0..num_outputs {
            nodes.push(NodeGene::output((num_inputs + 1 + i) as u32));
        }

        // Create connections from all inputs (including bias) to all outputs
        let mut connections = Vec::new();
        for i in 0..=num_inputs {
            for j in 0..num_outputs {
                let from = i as u32;
                let to = (num_inputs + 1 + j) as u32;
                let inno = tracker.get_connection_innovation(from, to);
                connections.push(ConnectionGene::new(inno, from, to, 0.0));
            }
        }

        Self {
            nodes,
            connections,
            num_inputs,
            num_outputs,
            fitness: 0.0,
            adjusted_fitness: 0.0,
            species_id: 0,
        }
    }

    /// Create empty genome (for crossover)
    pub fn empty(num_inputs: usize, num_outputs: usize) -> Self {
        Self {
            nodes: Vec::new(),
            connections: Vec::new(),
            num_inputs,
            num_outputs,
            fitness: 0.0,
            adjusted_fitness: 0.0,
            species_id: 0,
        }
    }

    /// Get node by ID
    pub fn get_node(&self, id: u32) -> Option<&NodeGene> {
        self.nodes.iter().find(|n| n.id == id)
    }

    /// Get mutable node by ID
    pub fn get_node_mut(&mut self, id: u32) -> Option<&mut NodeGene> {
        self.nodes.iter_mut().find(|n| n.id == id)
    }

    /// Get connection by innovation number
    pub fn get_connection(&self, innovation: u32) -> Option<&ConnectionGene> {
        self.connections.iter().find(|c| c.innovation == innovation)
    }

    /// Check if connection exists
    pub fn has_connection(&self, from: u32, to: u32) -> bool {
        self.connections.iter().any(|c| c.from == from && c.to == to)
    }

    /// Get number of hidden nodes
    pub fn num_hidden(&self) -> usize {
        self.nodes.iter().filter(|n| n.node_type == NodeType::Hidden).count()
    }

    /// Get number of enabled connections
    pub fn num_enabled_connections(&self) -> usize {
        self.connections.iter().filter(|c| c.enabled).count()
    }

    /// Get all input node IDs
    pub fn input_ids(&self) -> Vec<u32> {
        self.nodes
            .iter()
            .filter(|n| n.node_type == NodeType::Input)
            .map(|n| n.id)
            .collect()
    }

    /// Get all output node IDs
    pub fn output_ids(&self) -> Vec<u32> {
        self.nodes
            .iter()
            .filter(|n| n.node_type == NodeType::Output)
            .map(|n| n.id)
            .collect()
    }

    /// Get all hidden node IDs
    pub fn hidden_ids(&self) -> Vec<u32> {
        self.nodes
            .iter()
            .filter(|n| n.node_type == NodeType::Hidden)
            .map(|n| n.id)
            .collect()
    }

    /// Check if adding connection would create a cycle (for feedforward networks)
    pub fn would_create_cycle(&self, from: u32, to: u32) -> bool {
        if from == to {
            return true;
        }

        // BFS to check if 'to' can reach 'from'
        let mut visited = HashSet::new();
        let mut queue = vec![to];

        while let Some(node) = queue.pop() {
            if node == from {
                return true;
            }

            if visited.contains(&node) {
                continue;
            }
            visited.insert(node);

            for conn in &self.connections {
                if conn.enabled && conn.from == node {
                    queue.push(conn.to);
                }
            }
        }

        false
    }

    /// Compute compatibility distance with another genome
    pub fn distance(&self, other: &Genome, c1: f32, c2: f32, c3: f32) -> f32 {
        let max_genes = self.connections.len().max(other.connections.len()) as f32;
        if max_genes == 0.0 {
            return 0.0;
        }

        let mut matching = 0;
        let mut disjoint = 0;
        let mut excess = 0;
        let mut weight_diff = 0.0;

        let max_inno_self = self.connections.iter().map(|c| c.innovation).max().unwrap_or(0);
        let max_inno_other = other.connections.iter().map(|c| c.innovation).max().unwrap_or(0);

        let self_innovations: HashSet<_> = self.connections.iter().map(|c| c.innovation).collect();
        let other_innovations: HashSet<_> = other.connections.iter().map(|c| c.innovation).collect();

        for conn in &self.connections {
            if other_innovations.contains(&conn.innovation) {
                matching += 1;
                if let Some(other_conn) = other.get_connection(conn.innovation) {
                    weight_diff += (conn.weight - other_conn.weight).abs();
                }
            } else if conn.innovation > max_inno_other {
                excess += 1;
            } else {
                disjoint += 1;
            }
        }

        for conn in &other.connections {
            if !self_innovations.contains(&conn.innovation) {
                if conn.innovation > max_inno_self {
                    excess += 1;
                } else {
                    disjoint += 1;
                }
            }
        }

        let avg_weight_diff = if matching > 0 {
            weight_diff / matching as f32
        } else {
            0.0
        };

        (c1 * excess as f32 + c2 * disjoint as f32) / max_genes + c3 * avg_weight_diff
    }

    /// Sort connections by innovation number
    pub fn sort_connections(&mut self) {
        self.connections.sort_by_key(|c| c.innovation);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_minimal_genome() {
        let mut tracker = InnovationTracker::new();
        let genome = Genome::minimal(3, 2, &mut tracker);

        assert_eq!(genome.nodes.len(), 3 + 1 + 2); // inputs + bias + outputs
        assert_eq!(genome.connections.len(), 4 * 2); // (inputs + bias) * outputs
    }

    #[test]
    fn test_innovation_tracking() {
        let mut tracker = InnovationTracker::new();

        let i1 = tracker.get_connection_innovation(0, 3);
        let i2 = tracker.get_connection_innovation(1, 3);
        let i3 = tracker.get_connection_innovation(0, 3); // Same as i1

        assert_eq!(i1, i3); // Same connection should have same innovation
        assert_ne!(i1, i2); // Different connections should have different innovations
    }

    #[test]
    fn test_cycle_detection() {
        let mut tracker = InnovationTracker::new();
        let mut genome = Genome::minimal(2, 1, &mut tracker);

        // Minimal genome is feedforward, no cycles
        assert!(!genome.would_create_cycle(0, 3)); // input to output

        // Adding output to hidden would create cycle if we have hidden
        // For now, just test self-loop
        assert!(genome.would_create_cycle(0, 0));
    }
}
