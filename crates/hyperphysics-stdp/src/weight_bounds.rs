//! Weight Bounds and Update Tracking
//!
//! Manages synaptic weight constraints and accumulates updates.

#[cfg(feature = "serde")]
use serde::{Deserialize, Serialize};

/// Weight bounds configuration
#[derive(Debug, Clone, Copy)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct WeightBounds {
    /// Minimum weight
    pub min: f32,

    /// Maximum weight
    pub max: f32,

    /// Soft bounds: apply multiplicative scaling near limits
    pub soft_bounds: bool,

    /// Soft bound steepness (higher = sharper transition)
    pub soft_steepness: f32,
}

impl Default for WeightBounds {
    fn default() -> Self {
        Self {
            min: 0.0,
            max: 10.0,
            soft_bounds: true,
            soft_steepness: 2.0,
        }
    }
}

impl WeightBounds {
    /// Hard bounds only
    pub fn hard(min: f32, max: f32) -> Self {
        Self {
            min,
            max,
            soft_bounds: false,
            soft_steepness: 0.0,
        }
    }

    /// Clamp weight to bounds
    pub fn clamp(&self, weight: f32) -> f32 {
        weight.clamp(self.min, self.max)
    }

    /// Apply soft bounds to weight update
    /// Returns scaled update that naturally respects bounds
    pub fn apply_soft_bounds(&self, weight: f32, delta: f32) -> f32 {
        if !self.soft_bounds {
            return delta;
        }

        let range = self.max - self.min;
        if range <= 0.0 {
            return 0.0;
        }

        // Normalize weight to 0-1 range
        let w_norm = (weight - self.min) / range;

        // Apply multiplicative soft bounds
        let scale = if delta > 0.0 {
            // Potentiation: reduce as approaching max
            (1.0 - w_norm).powf(self.soft_steepness)
        } else {
            // Depression: reduce as approaching min
            w_norm.powf(self.soft_steepness)
        };

        delta * scale
    }

    /// Check if weight is at upper bound
    pub fn at_upper(&self, weight: f32) -> bool {
        (weight - self.max).abs() < 1e-6
    }

    /// Check if weight is at lower bound
    pub fn at_lower(&self, weight: f32) -> bool {
        (weight - self.min).abs() < 1e-6
    }
}

/// A single weight update
#[derive(Debug, Clone, Copy)]
pub struct WeightUpdate {
    /// Synapse index
    pub synapse_id: usize,

    /// Weight change (can be positive or negative)
    pub delta: f32,
}

impl WeightUpdate {
    /// Create new update
    pub fn new(synapse_id: usize, delta: f32) -> Self {
        Self { synapse_id, delta }
    }

    /// Check if this is LTP (weight increase)
    pub fn is_ltp(&self) -> bool {
        self.delta > 0.0
    }

    /// Check if this is LTD (weight decrease)
    pub fn is_ltd(&self) -> bool {
        self.delta < 0.0
    }
}

/// Accumulator for weight updates
/// Collects updates before batch application
#[derive(Debug, Clone, Default)]
pub struct UpdateAccumulator {
    /// Pending updates per synapse
    updates: Vec<f32>,

    /// Count of updates per synapse
    counts: Vec<u32>,

    /// Total LTP events
    ltp_count: u64,

    /// Total LTD events
    ltd_count: u64,
}

impl UpdateAccumulator {
    /// Create new accumulator
    pub fn new(num_synapses: usize) -> Self {
        Self {
            updates: vec![0.0; num_synapses],
            counts: vec![0; num_synapses],
            ltp_count: 0,
            ltd_count: 0,
        }
    }

    /// Add an update
    pub fn add(&mut self, synapse_id: usize, delta: f32) {
        if synapse_id < self.updates.len() {
            self.updates[synapse_id] += delta;
            self.counts[synapse_id] += 1;

            if delta > 0.0 {
                self.ltp_count += 1;
            } else if delta < 0.0 {
                self.ltd_count += 1;
            }
        }
    }

    /// Get accumulated update for a synapse
    pub fn get(&self, synapse_id: usize) -> f32 {
        self.updates.get(synapse_id).copied().unwrap_or(0.0)
    }

    /// Get all non-zero updates
    pub fn get_updates(&self) -> Vec<WeightUpdate> {
        self.updates
            .iter()
            .enumerate()
            .filter(|(_, &delta)| delta.abs() > 1e-9)
            .map(|(id, &delta)| WeightUpdate::new(id, delta))
            .collect()
    }

    /// Apply all updates to weights
    pub fn apply(&self, weights: &mut [f32], bounds: &WeightBounds) {
        for (i, &delta) in self.updates.iter().enumerate() {
            if i < weights.len() && delta.abs() > 1e-9 {
                let scaled_delta = bounds.apply_soft_bounds(weights[i], delta);
                weights[i] = bounds.clamp(weights[i] + scaled_delta);
            }
        }
    }

    /// Clear accumulated updates
    pub fn clear(&mut self) {
        self.updates.fill(0.0);
        self.counts.fill(0);
    }

    /// Reset including statistics
    pub fn reset(&mut self) {
        self.clear();
        self.ltp_count = 0;
        self.ltd_count = 0;
    }

    /// Get LTP count
    pub fn ltp_count(&self) -> u64 {
        self.ltp_count
    }

    /// Get LTD count
    pub fn ltd_count(&self) -> u64 {
        self.ltd_count
    }

    /// Get average update magnitude
    pub fn average_magnitude(&self) -> f32 {
        let non_zero: Vec<_> = self.updates.iter().filter(|&&x| x.abs() > 1e-9).collect();
        if non_zero.is_empty() {
            return 0.0;
        }
        let sum: f32 = non_zero.iter().map(|x| x.abs()).sum();
        sum / non_zero.len() as f32
    }
}

/// Weight normalization methods
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub enum NormalizationMethod {
    /// No normalization
    None,

    /// Multiplicative normalization (total weight constant)
    Multiplicative,

    /// Subtractive normalization (mean weight constant)
    Subtractive,

    /// Both multiplicative and subtractive
    Both,
}

/// Normalize weights to maintain homeostasis
pub fn normalize_weights(
    weights: &mut [f32],
    method: NormalizationMethod,
    target_sum: f32,
    target_mean: f32,
) {
    if weights.is_empty() {
        return;
    }

    match method {
        NormalizationMethod::None => {}

        NormalizationMethod::Multiplicative => {
            let sum: f32 = weights.iter().sum();
            if sum > 1e-6 {
                let scale = target_sum / sum;
                for w in weights {
                    *w *= scale;
                }
            }
        }

        NormalizationMethod::Subtractive => {
            let mean: f32 = weights.iter().sum::<f32>() / weights.len() as f32;
            let delta = target_mean - mean;
            for w in weights {
                *w += delta;
            }
        }

        NormalizationMethod::Both => {
            // First subtractive (center), then multiplicative (scale)
            let mean: f32 = weights.iter().sum::<f32>() / weights.len() as f32;
            let delta = target_mean - mean;
            for w in weights {
                *w += delta;
            }

            let sum: f32 = weights.iter().sum();
            if sum > 1e-6 {
                let scale = target_sum / sum;
                for w in weights {
                    *w *= scale;
                }
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_weight_bounds() {
        let bounds = WeightBounds::default();

        assert!((bounds.clamp(-1.0) - 0.0).abs() < 1e-6);
        assert!((bounds.clamp(15.0) - 10.0).abs() < 1e-6);
        assert!((bounds.clamp(5.0) - 5.0).abs() < 1e-6);
    }

    #[test]
    fn test_soft_bounds() {
        let bounds = WeightBounds::default();

        // Near max, potentiation should be reduced
        let delta_near_max = bounds.apply_soft_bounds(9.0, 1.0);
        assert!(delta_near_max < 1.0);

        // Near min, depression should be reduced
        let delta_near_min = bounds.apply_soft_bounds(1.0, -1.0);
        assert!(delta_near_min > -1.0);
    }

    #[test]
    fn test_update_accumulator() {
        let mut acc = UpdateAccumulator::new(10);

        acc.add(0, 0.5);
        acc.add(0, 0.3);
        acc.add(1, -0.2);

        let updates = acc.get_updates();
        assert_eq!(updates.len(), 2);

        assert_eq!(acc.ltp_count(), 2);
        assert_eq!(acc.ltd_count(), 1);
    }

    #[test]
    fn test_normalization() {
        let mut weights = vec![1.0, 2.0, 3.0, 4.0];

        normalize_weights(&mut weights, NormalizationMethod::Multiplicative, 20.0, 0.0);
        let sum: f32 = weights.iter().sum();
        assert!((sum - 20.0).abs() < 1e-6);
    }
}
