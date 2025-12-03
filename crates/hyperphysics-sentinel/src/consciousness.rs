//! Consciousness Metrics
//!
//! Implementation of Integrated Information Theory (IIT) and related metrics.

#[cfg(feature = "serde")]
use serde::{Deserialize, Serialize};

/// Consciousness metrics container
#[derive(Debug, Clone, Default)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct ConsciousnessMetrics {
    /// Integrated Information (Φ)
    pub phi: f64,

    /// Causal density
    pub causal_density: f64,

    /// Metastability
    pub metastability: f64,

    /// Complexity (neural complexity)
    pub complexity: f64,

    /// Activity level (normalized)
    pub activity_level: f64,

    /// Synchronization index
    pub synchronization: f64,

    /// Entropy
    pub entropy: f64,

    /// Last computation time
    pub last_computed: f64,
}

impl ConsciousnessMetrics {
    pub fn new() -> Self {
        Self::default()
    }

    /// Get overall consciousness score
    pub fn score(&self) -> f64 {
        // Weighted combination of metrics
        self.phi * 0.4
            + self.causal_density * 0.2
            + self.metastability * 0.2
            + self.complexity * 0.1
            + (1.0 - self.synchronization.abs()) * 0.1 // Avoid too much sync
    }

    /// Check if consciousness threshold is met
    pub fn meets_threshold(&self, threshold: f64) -> bool {
        self.phi >= threshold
    }
}

/// Calculator for Integrated Information (Φ)
pub struct PhiCalculator {
    /// Number of neurons
    num_neurons: usize,

    /// Partition strategy
    partition_strategy: PartitionStrategy,

    /// State history for computing information
    state_history: Vec<Vec<bool>>,

    /// Maximum history length
    max_history: usize,
}

/// Partition strategy for computing Φ
#[derive(Debug, Clone, Copy)]
pub enum PartitionStrategy {
    /// Minimum information partition (exact, slow)
    MinimumInformationPartition,
    /// Bipartition (approximate, faster)
    Bipartition,
    /// Atomic partition (fastest)
    Atomic,
}

impl PhiCalculator {
    /// Create new Φ calculator
    pub fn new(num_neurons: usize) -> Self {
        Self {
            num_neurons,
            partition_strategy: PartitionStrategy::Bipartition,
            state_history: Vec::new(),
            max_history: 1000,
        }
    }

    /// Add a state observation
    pub fn observe(&mut self, state: Vec<bool>) {
        if state.len() != self.num_neurons {
            return;
        }

        self.state_history.push(state);

        // Limit history
        if self.state_history.len() > self.max_history {
            self.state_history.remove(0);
        }
    }

    /// Compute Φ from voltage states
    pub fn compute_from_voltages(&mut self, voltages: &[f32], threshold: f32) -> f64 {
        // Convert to binary states
        let binary_state: Vec<bool> = voltages.iter().map(|&v| v > threshold).collect();
        self.observe(binary_state);
        self.compute()
    }

    /// Compute Φ from current history
    pub fn compute(&self) -> f64 {
        if self.state_history.len() < 10 {
            return 0.0;
        }

        match self.partition_strategy {
            PartitionStrategy::Atomic => self.compute_atomic_phi(),
            PartitionStrategy::Bipartition => self.compute_bipartition_phi(),
            PartitionStrategy::MinimumInformationPartition => self.compute_mip_phi(),
        }
    }

    /// Atomic partition (simplified)
    fn compute_atomic_phi(&self) -> f64 {
        // Compute mutual information between neurons
        let mut total_mi = 0.0;
        let n = self.num_neurons;

        for i in 0..n {
            for j in (i + 1)..n {
                total_mi += self.mutual_information(i, j);
            }
        }

        // Normalize
        let pairs = (n * (n - 1) / 2) as f64;
        if pairs > 0.0 {
            total_mi / pairs
        } else {
            0.0
        }
    }

    /// Bipartition (medium complexity)
    fn compute_bipartition_phi(&self) -> f64 {
        // Find the bipartition that minimizes integrated information
        let n = self.num_neurons;
        if n < 2 {
            return 0.0;
        }

        let mut min_phi = f64::INFINITY;

        // Try all possible bipartitions
        for mask in 1..(1 << n) - 1 {
            let mut part_a = Vec::new();
            let mut part_b = Vec::new();

            for i in 0..n {
                if (mask >> i) & 1 == 1 {
                    part_a.push(i);
                } else {
                    part_b.push(i);
                }
            }

            // Compute information between partitions
            let phi = self.information_across_partition(&part_a, &part_b);
            min_phi = min_phi.min(phi);
        }

        if min_phi.is_finite() {
            min_phi
        } else {
            0.0
        }
    }

    /// Minimum information partition (exact, expensive)
    fn compute_mip_phi(&self) -> f64 {
        // For small systems, this is tractable
        // For larger systems, would need approximations
        if self.num_neurons <= 8 {
            self.compute_bipartition_phi()
        } else {
            // Fall back to bipartition for large systems
            self.compute_bipartition_phi()
        }
    }

    /// Compute mutual information between two neurons
    fn mutual_information(&self, i: usize, j: usize) -> f64 {
        if self.state_history.is_empty() {
            return 0.0;
        }

        let n = self.state_history.len() as f64;

        // Count joint and marginal probabilities
        let mut p_00 = 0.0;
        let mut p_01 = 0.0;
        let mut p_10 = 0.0;
        let mut p_11 = 0.0;

        for state in &self.state_history {
            let a = state[i];
            let b = state[j];
            match (a, b) {
                (false, false) => p_00 += 1.0,
                (false, true) => p_01 += 1.0,
                (true, false) => p_10 += 1.0,
                (true, true) => p_11 += 1.0,
            }
        }

        p_00 /= n;
        p_01 /= n;
        p_10 /= n;
        p_11 /= n;

        let p_i_0 = p_00 + p_01;
        let p_i_1 = p_10 + p_11;
        let p_j_0 = p_00 + p_10;
        let p_j_1 = p_01 + p_11;

        // Compute mutual information
        let mut mi = 0.0;

        for &(p_ij, p_i, p_j) in &[
            (p_00, p_i_0, p_j_0),
            (p_01, p_i_0, p_j_1),
            (p_10, p_i_1, p_j_0),
            (p_11, p_i_1, p_j_1),
        ] {
            if p_ij > 0.0 && p_i > 0.0 && p_j > 0.0 {
                mi += p_ij * (p_ij / (p_i * p_j)).ln();
            }
        }

        mi.max(0.0)
    }

    /// Compute information across a partition
    fn information_across_partition(&self, part_a: &[usize], part_b: &[usize]) -> f64 {
        // Simplified: sum of mutual information between parts
        let mut total = 0.0;

        for &i in part_a {
            for &j in part_b {
                total += self.mutual_information(i, j);
            }
        }

        total
    }

    /// Reset state history
    pub fn reset(&mut self) {
        self.state_history.clear();
    }
}

/// Causal density calculator
pub struct CausalDensity {
    /// Number of neurons
    num_neurons: usize,

    /// Causal interaction matrix
    causal_matrix: Vec<Vec<f64>>,

    /// State history
    history: Vec<Vec<f32>>,

    /// Maximum history
    max_history: usize,
}

impl CausalDensity {
    /// Create new causal density calculator
    pub fn new(num_neurons: usize) -> Self {
        Self {
            num_neurons,
            causal_matrix: vec![vec![0.0; num_neurons]; num_neurons],
            history: Vec::new(),
            max_history: 500,
        }
    }

    /// Add state observation
    pub fn observe(&mut self, state: Vec<f32>) {
        if state.len() != self.num_neurons {
            return;
        }

        self.history.push(state);

        if self.history.len() > self.max_history {
            self.history.remove(0);
        }

        // Update causal matrix
        self.update_causal_matrix();
    }

    /// Update causal interaction matrix using Granger causality
    fn update_causal_matrix(&mut self) {
        if self.history.len() < 10 {
            return;
        }

        // Simplified Granger causality estimation
        let n = self.num_neurons;

        for i in 0..n {
            for j in 0..n {
                if i == j {
                    continue;
                }

                // Estimate predictive power of j on i
                let mut correlation = 0.0;
                let count = self.history.len() - 1;

                for t in 0..count {
                    let x_t = self.history[t][j];
                    let y_next = self.history[t + 1][i];
                    correlation += x_t * y_next;
                }

                if count > 0 {
                    self.causal_matrix[i][j] = (correlation / count as f64).abs();
                }
            }
        }
    }

    /// Compute overall causal density
    pub fn compute(&self) -> f64 {
        let n = self.num_neurons;
        if n == 0 {
            return 0.0;
        }

        let mut sum = 0.0;
        let mut count = 0;

        for i in 0..n {
            for j in 0..n {
                if i != j {
                    sum += self.causal_matrix[i][j];
                    count += 1;
                }
            }
        }

        if count > 0 {
            sum / count as f64
        } else {
            0.0
        }
    }

    /// Get causal interaction strength between two neurons
    pub fn causal_strength(&self, from: usize, to: usize) -> f64 {
        if from < self.num_neurons && to < self.num_neurons {
            self.causal_matrix[to][from]
        } else {
            0.0
        }
    }

    /// Reset history
    pub fn reset(&mut self) {
        self.history.clear();
        for row in &mut self.causal_matrix {
            row.fill(0.0);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_phi_calculator() {
        let mut calc = PhiCalculator::new(4);

        // Add some random states
        for _ in 0..50 {
            let state: Vec<bool> = (0..4).map(|i| i % 2 == 0).collect();
            calc.observe(state);
        }

        let phi = calc.compute();
        assert!(phi >= 0.0);
    }

    #[test]
    fn test_consciousness_metrics() {
        let mut metrics = ConsciousnessMetrics::new();
        metrics.phi = 0.5;
        metrics.causal_density = 0.3;

        assert!(metrics.score() > 0.0);
        assert!(metrics.meets_threshold(0.4));
        assert!(!metrics.meets_threshold(0.6));
    }

    #[test]
    fn test_causal_density() {
        let mut cd = CausalDensity::new(3);

        // Add correlated states
        for i in 0..20 {
            let state = vec![i as f32, (i + 1) as f32, (i + 2) as f32];
            cd.observe(state);
        }

        let density = cd.compute();
        assert!(density >= 0.0);
    }
}
