//! Fitness Functions
//!
//! Evaluate genome performance for evolution.

#[cfg(feature = "serde")]
use serde::{Deserialize, Serialize};

use crate::genome::Genome;

/// A fitness metric type
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub enum FitnessMetric {
    /// Task performance (primary)
    Performance,

    /// Network complexity (minimize)
    Complexity,

    /// Integrated information (Φ)
    Phi,

    /// Energy efficiency
    Energy,

    /// Behavioral diversity
    Diversity,

    /// Robustness to perturbation
    Robustness,

    /// Speed/latency
    Speed,

    /// Custom metric
    Custom(u32),
}

/// Multi-objective fitness
#[derive(Debug, Clone, Default)]
pub struct MultiFitness {
    /// Fitness values for each objective
    values: Vec<(FitnessMetric, f32)>,

    /// Pareto rank (lower is better)
    pub pareto_rank: usize,

    /// Crowding distance (for diversity)
    pub crowding_distance: f32,
}

impl MultiFitness {
    pub fn new() -> Self {
        Self::default()
    }

    /// Set a fitness value
    pub fn set(&mut self, metric: FitnessMetric, value: f32) {
        if let Some((_, v)) = self.values.iter_mut().find(|(m, _)| *m == metric) {
            *v = value;
        } else {
            self.values.push((metric, value));
        }
    }

    /// Get a fitness value
    pub fn get(&self, metric: FitnessMetric) -> Option<f32> {
        self.values.iter().find(|(m, _)| *m == metric).map(|(_, v)| *v)
    }

    /// Get primary fitness (first metric)
    pub fn primary(&self) -> f32 {
        self.values.first().map(|(_, v)| *v).unwrap_or(0.0)
    }

    /// Check if this dominates another (for Pareto)
    pub fn dominates(&self, other: &MultiFitness) -> bool {
        let mut dominated_in_one = false;

        for (metric, value) in &self.values {
            if let Some(other_value) = other.get(*metric) {
                // Assuming higher is better for all metrics
                if *value < other_value {
                    return false; // Other is better in at least one
                }
                if *value > other_value {
                    dominated_in_one = true;
                }
            }
        }

        dominated_in_one
    }

    /// Compute scalar fitness (weighted sum)
    pub fn scalar(&self, weights: &[(FitnessMetric, f32)]) -> f32 {
        let mut sum = 0.0;
        for (metric, weight) in weights {
            if let Some(value) = self.get(*metric) {
                sum += value * weight;
            }
        }
        sum
    }
}

/// Trait for fitness evaluation functions
pub trait FitnessFunction: Send + Sync {
    /// Evaluate a genome and return its fitness
    fn evaluate(&self, genome: &Genome) -> f32;

    /// Evaluate with multi-objective fitness
    fn evaluate_multi(&self, genome: &Genome) -> MultiFitness {
        let mut fitness = MultiFitness::new();
        fitness.set(FitnessMetric::Performance, self.evaluate(genome));
        fitness
    }

    /// Get the metrics this function evaluates
    fn metrics(&self) -> Vec<FitnessMetric> {
        vec![FitnessMetric::Performance]
    }
}

/// Simple fitness function wrapper
pub struct SimpleFitness<F>
where
    F: Fn(&Genome) -> f32 + Send + Sync,
{
    func: F,
}

impl<F> SimpleFitness<F>
where
    F: Fn(&Genome) -> f32 + Send + Sync,
{
    pub fn new(func: F) -> Self {
        Self { func }
    }
}

impl<F> FitnessFunction for SimpleFitness<F>
where
    F: Fn(&Genome) -> f32 + Send + Sync,
{
    fn evaluate(&self, genome: &Genome) -> f32 {
        (self.func)(genome)
    }
}

/// Complexity penalty fitness
pub struct ComplexityPenalizedFitness<F>
where
    F: FitnessFunction,
{
    inner: F,
    node_penalty: f32,
    connection_penalty: f32,
}

impl<F> ComplexityPenalizedFitness<F>
where
    F: FitnessFunction,
{
    pub fn new(inner: F, node_penalty: f32, connection_penalty: f32) -> Self {
        Self {
            inner,
            node_penalty,
            connection_penalty,
        }
    }
}

impl<F> FitnessFunction for ComplexityPenalizedFitness<F>
where
    F: FitnessFunction,
{
    fn evaluate(&self, genome: &Genome) -> f32 {
        let base_fitness = self.inner.evaluate(genome);
        let node_cost = genome.num_hidden() as f32 * self.node_penalty;
        let conn_cost = genome.num_enabled_connections() as f32 * self.connection_penalty;
        (base_fitness - node_cost - conn_cost).max(0.0)
    }

    fn evaluate_multi(&self, genome: &Genome) -> MultiFitness {
        let mut fitness = self.inner.evaluate_multi(genome);
        let complexity = genome.num_hidden() as f32 + genome.num_enabled_connections() as f32 * 0.1;
        fitness.set(FitnessMetric::Complexity, -complexity); // Negative because lower is better
        fitness
    }

    fn metrics(&self) -> Vec<FitnessMetric> {
        let mut metrics = self.inner.metrics();
        metrics.push(FitnessMetric::Complexity);
        metrics
    }
}

/// Novelty search fitness
pub struct NoveltyFitness {
    /// Archive of novel behaviors
    archive: Vec<Vec<f32>>,

    /// Behavior dimensionality
    behavior_dim: usize,

    /// k for k-nearest neighbors
    k: usize,

    /// Archive threshold
    threshold: f32,
}

impl NoveltyFitness {
    pub fn new(behavior_dim: usize) -> Self {
        Self {
            archive: Vec::new(),
            behavior_dim,
            k: 15,
            threshold: 0.1,
        }
    }

    /// Extract behavior from genome (override in subclass)
    fn extract_behavior(&self, genome: &Genome) -> Vec<f32> {
        // Default: use genome structure as behavior
        let mut behavior = vec![0.0; self.behavior_dim];

        if !behavior.is_empty() {
            behavior[0] = genome.num_hidden() as f32;
        }
        if behavior.len() > 1 {
            behavior[1] = genome.num_enabled_connections() as f32;
        }

        behavior
    }

    /// Compute novelty (distance to nearest neighbors)
    fn compute_novelty(&self, behavior: &[f32]) -> f32 {
        if self.archive.is_empty() {
            return 1.0;
        }

        let mut distances: Vec<f32> = self
            .archive
            .iter()
            .map(|b| {
                b.iter()
                    .zip(behavior.iter())
                    .map(|(a, b)| (a - b).powi(2))
                    .sum::<f32>()
                    .sqrt()
            })
            .collect();

        distances.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));

        let k = self.k.min(distances.len());
        distances[..k].iter().sum::<f32>() / k as f32
    }

    /// Add behavior to archive if novel enough
    pub fn maybe_add_to_archive(&mut self, behavior: Vec<f32>) {
        let novelty = self.compute_novelty(&behavior);
        if novelty > self.threshold {
            self.archive.push(behavior);
        }
    }
}

impl FitnessFunction for NoveltyFitness {
    fn evaluate(&self, genome: &Genome) -> f32 {
        let behavior = self.extract_behavior(genome);
        self.compute_novelty(&behavior)
    }

    fn metrics(&self) -> Vec<FitnessMetric> {
        vec![FitnessMetric::Diversity]
    }
}

/// NSGA-II Pareto ranking
pub fn pareto_rank(population: &mut [(usize, MultiFitness)]) {
    let n = population.len();
    if n == 0 {
        return;
    }

    // Count dominations
    let mut domination_count = vec![0; n];
    let mut dominated_by: Vec<Vec<usize>> = vec![Vec::new(); n];

    for i in 0..n {
        for j in 0..n {
            if i == j {
                continue;
            }

            if population[i].1.dominates(&population[j].1) {
                dominated_by[i].push(j);
            } else if population[j].1.dominates(&population[i].1) {
                domination_count[i] += 1;
            }
        }
    }

    // Assign ranks
    let mut rank = 0;
    let mut front: Vec<usize> = (0..n)
        .filter(|&i| domination_count[i] == 0)
        .collect();

    while !front.is_empty() {
        let mut next_front = Vec::new();

        for &i in &front {
            population[i].1.pareto_rank = rank;

            for &j in &dominated_by[i] {
                domination_count[j] -= 1;
                if domination_count[j] == 0 {
                    next_front.push(j);
                }
            }
        }

        rank += 1;
        front = next_front;
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::genome::InnovationTracker;

    #[test]
    fn test_multi_fitness() {
        let mut f1 = MultiFitness::new();
        f1.set(FitnessMetric::Performance, 1.0);
        f1.set(FitnessMetric::Complexity, 0.5);

        let mut f2 = MultiFitness::new();
        f2.set(FitnessMetric::Performance, 0.8);
        f2.set(FitnessMetric::Complexity, 0.5);

        assert!(f1.dominates(&f2)); // f1 is better in Performance, equal in Complexity
        assert!(!f2.dominates(&f1));
    }

    #[test]
    fn test_simple_fitness() {
        let mut tracker = InnovationTracker::new();
        let genome = Genome::minimal(2, 1, &mut tracker);

        let fitness_fn = SimpleFitness::new(|g: &Genome| g.connections.len() as f32);
        let score = fitness_fn.evaluate(&genome);

        assert!(score > 0.0);
    }
}
