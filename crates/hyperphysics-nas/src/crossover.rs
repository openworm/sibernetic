//! Crossover Operators
//!
//! Combine two parent genomes to create offspring.

#[cfg(feature = "serde")]
use serde::{Deserialize, Serialize};

use rand::Rng;
use hashbrown::HashSet;

use crate::genome::{ConnectionGene, Genome, NodeGene};

/// Crossover parameters
#[derive(Debug, Clone)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct CrossoverParams {
    /// Probability of inheriting gene from less fit parent when genes match
    pub crossover_rate: f32,

    /// Probability of disabling gene if disabled in either parent
    pub disable_gene_rate: f32,
}

impl Default for CrossoverParams {
    fn default() -> Self {
        Self {
            crossover_rate: 0.5,
            disable_gene_rate: 0.75,
        }
    }
}

/// Crossover operator
pub struct CrossoverOperator {
    params: CrossoverParams,
}

impl CrossoverOperator {
    pub fn new(params: CrossoverParams) -> Self {
        Self { params }
    }

    /// Perform crossover between two parents
    ///
    /// The first parent should be the more fit one.
    /// Disjoint and excess genes are inherited from the fitter parent.
    pub fn crossover<R: Rng>(
        &self,
        parent1: &Genome,
        parent2: &Genome,
        rng: &mut R,
    ) -> Genome {
        let mut offspring = Genome::empty(parent1.num_inputs, parent1.num_outputs);

        // Collect all node IDs from both parents
        let mut node_ids: HashSet<u32> = HashSet::new();

        // Collect innovation numbers from both parents
        let p1_innovations: HashSet<u32> = parent1.connections.iter().map(|c| c.innovation).collect();
        let p2_innovations: HashSet<u32> = parent2.connections.iter().map(|c| c.innovation).collect();

        // Process connections
        for conn in &parent1.connections {
            if p2_innovations.contains(&conn.innovation) {
                // Matching gene - randomly inherit
                let inherited = if rng.gen::<f32>() < self.params.crossover_rate {
                    conn.clone()
                } else {
                    parent2.get_connection(conn.innovation).unwrap().clone()
                };

                // Check if should be disabled
                let p1_enabled = conn.enabled;
                let p2_enabled = parent2.get_connection(conn.innovation).map(|c| c.enabled).unwrap_or(true);

                let mut final_conn = inherited;
                if !p1_enabled || !p2_enabled {
                    if rng.gen::<f32>() < self.params.disable_gene_rate {
                        final_conn.enabled = false;
                    }
                }

                node_ids.insert(final_conn.from);
                node_ids.insert(final_conn.to);
                offspring.connections.push(final_conn);
            } else {
                // Disjoint or excess - inherit from fitter parent (parent1)
                node_ids.insert(conn.from);
                node_ids.insert(conn.to);
                offspring.connections.push(conn.clone());
            }
        }

        // Inherit nodes that are referenced by connections
        for &node_id in &node_ids {
            if let Some(node) = parent1.get_node(node_id) {
                offspring.nodes.push(node.clone());
            } else if let Some(node) = parent2.get_node(node_id) {
                offspring.nodes.push(node.clone());
            }
        }

        // Sort nodes by ID
        offspring.nodes.sort_by_key(|n| n.id);

        // Sort connections by innovation
        offspring.sort_connections();

        offspring
    }

    /// Crossover with equal fitness parents
    /// Uses genes from both parents for disjoint/excess
    pub fn crossover_equal<R: Rng>(
        &self,
        parent1: &Genome,
        parent2: &Genome,
        rng: &mut R,
    ) -> Genome {
        let mut offspring = Genome::empty(parent1.num_inputs, parent1.num_outputs);

        let mut node_ids: HashSet<u32> = HashSet::new();

        let p1_innovations: HashSet<u32> = parent1.connections.iter().map(|c| c.innovation).collect();
        let p2_innovations: HashSet<u32> = parent2.connections.iter().map(|c| c.innovation).collect();

        // All innovations from both parents
        let all_innovations: HashSet<u32> = p1_innovations.union(&p2_innovations).copied().collect();

        for &inno in &all_innovations {
            let has_p1 = p1_innovations.contains(&inno);
            let has_p2 = p2_innovations.contains(&inno);

            let conn = if has_p1 && has_p2 {
                // Matching - random choice
                if rng.gen::<f32>() < 0.5 {
                    parent1.get_connection(inno).unwrap().clone()
                } else {
                    parent2.get_connection(inno).unwrap().clone()
                }
            } else if has_p1 {
                // Only in p1 - 50% chance to include
                if rng.gen::<f32>() < 0.5 {
                    parent1.get_connection(inno).unwrap().clone()
                } else {
                    continue;
                }
            } else {
                // Only in p2 - 50% chance to include
                if rng.gen::<f32>() < 0.5 {
                    parent2.get_connection(inno).unwrap().clone()
                } else {
                    continue;
                }
            };

            node_ids.insert(conn.from);
            node_ids.insert(conn.to);
            offspring.connections.push(conn);
        }

        // Inherit nodes
        for &node_id in &node_ids {
            if let Some(node) = parent1.get_node(node_id) {
                offspring.nodes.push(node.clone());
            } else if let Some(node) = parent2.get_node(node_id) {
                offspring.nodes.push(node.clone());
            }
        }

        offspring.nodes.sort_by_key(|n| n.id);
        offspring.sort_connections();

        offspring
    }
}

/// Perform interspecies averaging
/// Creates a genome with averaged weights from matching genes
pub fn average_genomes(genomes: &[&Genome]) -> Option<Genome> {
    if genomes.is_empty() {
        return None;
    }

    let first = genomes[0];
    let mut result = first.clone();

    // Average weights for each connection
    for conn in &mut result.connections {
        let mut sum = conn.weight;
        let mut count = 1;

        for &genome in &genomes[1..] {
            if let Some(other_conn) = genome.get_connection(conn.innovation) {
                sum += other_conn.weight;
                count += 1;
            }
        }

        conn.weight = sum / count as f32;
    }

    // Average biases for each node
    for node in &mut result.nodes {
        let mut sum = node.bias;
        let mut count = 1;

        for &genome in &genomes[1..] {
            if let Some(other_node) = genome.get_node(node.id) {
                sum += other_node.bias;
                count += 1;
            }
        }

        node.bias = sum / count as f32;
    }

    Some(result)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::genome::InnovationTracker;

    #[test]
    fn test_crossover() {
        let mut rng = rand::thread_rng();
        let mut tracker = InnovationTracker::new();

        let parent1 = Genome::minimal(2, 1, &mut tracker);
        let parent2 = Genome::minimal(2, 1, &mut tracker);

        let operator = CrossoverOperator::new(CrossoverParams::default());
        let offspring = operator.crossover(&parent1, &parent2, &mut rng);

        // Offspring should have nodes and connections
        assert!(!offspring.nodes.is_empty());
        assert!(!offspring.connections.is_empty());
    }

    #[test]
    fn test_genome_averaging() {
        let mut tracker = InnovationTracker::new();

        let mut g1 = Genome::minimal(2, 1, &mut tracker);
        let mut g2 = Genome::minimal(2, 1, &mut tracker);

        // Set different weights
        for conn in &mut g1.connections {
            conn.weight = 1.0;
        }
        for conn in &mut g2.connections {
            conn.weight = 3.0;
        }

        let avg = average_genomes(&[&g1, &g2]).unwrap();

        // Average should be 2.0
        for conn in &avg.connections {
            assert!((conn.weight - 2.0).abs() < 1e-6);
        }
    }
}
