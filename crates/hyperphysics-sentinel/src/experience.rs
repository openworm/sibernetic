//! Experience Buffer
//!
//! Stores agent experiences for learning and replay.

#[cfg(feature = "serde")]
use serde::{Deserialize, Serialize};

use crate::SimTime;

/// A single experience tuple
#[derive(Debug, Clone)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct Experience {
    /// Time of experience
    pub time: SimTime,

    /// Observation/state
    pub state: Vec<f32>,

    /// Action taken
    pub action: Vec<f32>,

    /// Reward received
    pub reward: f64,

    /// Episode terminal flag
    pub done: bool,
}

impl Experience {
    pub fn new(time: SimTime, state: Vec<f32>, action: Vec<f32>) -> Self {
        Self {
            time,
            state,
            action,
            reward: 0.0,
            done: false,
        }
    }
}

/// Reward signal
#[derive(Debug, Clone, Copy)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct Reward {
    /// Reward value
    pub value: f64,

    /// Time of reward
    pub time: SimTime,

    /// Reward source/type
    pub source: RewardSource,
}

impl Reward {
    pub fn new(value: f64, time: SimTime) -> Self {
        Self {
            value,
            time,
            source: RewardSource::Environment,
        }
    }

    pub fn with_source(value: f64, time: SimTime, source: RewardSource) -> Self {
        Self { value, time, source }
    }
}

/// Source of reward signal
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub enum RewardSource {
    /// External environment reward
    Environment,
    /// Intrinsic curiosity reward
    Curiosity,
    /// Social reward from other agents
    Social,
    /// Survival/homeostasis reward
    Survival,
    /// Trading P&L
    Trading,
    /// Custom reward
    Custom(u32),
}

/// Circular buffer for storing experiences
#[derive(Debug, Clone)]
pub struct ExperienceBuffer {
    /// Experience storage
    buffer: Vec<Experience>,

    /// Maximum capacity
    capacity: usize,

    /// Write position
    write_pos: usize,

    /// Number of elements
    size: usize,

    /// Episode boundaries
    episode_starts: Vec<usize>,
}

impl ExperienceBuffer {
    /// Create new buffer with given capacity
    pub fn new(capacity: usize) -> Self {
        Self {
            buffer: Vec::with_capacity(capacity),
            capacity,
            write_pos: 0,
            size: 0,
            episode_starts: vec![0],
        }
    }

    /// Push new experience
    pub fn push(&mut self, experience: Experience) {
        if experience.done {
            self.episode_starts.push(self.size + 1);
        }

        if self.buffer.len() < self.capacity {
            self.buffer.push(experience);
        } else {
            self.buffer[self.write_pos] = experience;
        }

        self.write_pos = (self.write_pos + 1) % self.capacity;
        self.size = self.size.saturating_add(1).min(self.capacity);
    }

    /// Get experience at index
    pub fn get(&self, index: usize) -> Option<&Experience> {
        if index < self.size {
            let actual_idx = if self.size < self.capacity {
                index
            } else {
                (self.write_pos + index) % self.capacity
            };
            self.buffer.get(actual_idx)
        } else {
            None
        }
    }

    /// Get last experience (mutable)
    pub fn last_mut(&mut self) -> Option<&mut Experience> {
        if self.size == 0 {
            return None;
        }

        let idx = if self.size < self.capacity {
            self.size - 1
        } else {
            (self.write_pos + self.capacity - 1) % self.capacity
        };

        self.buffer.get_mut(idx)
    }

    /// Sample random batch
    pub fn sample(&self, batch_size: usize) -> Vec<&Experience> {
        use rand::seq::SliceRandom;
        let mut rng = rand::thread_rng();

        let indices: Vec<usize> = (0..self.size).collect();
        let sampled: Vec<usize> = indices
            .choose_multiple(&mut rng, batch_size.min(self.size))
            .copied()
            .collect();

        sampled
            .iter()
            .filter_map(|&i| self.get(i))
            .collect()
    }

    /// Sample sequence of experiences
    pub fn sample_sequence(&self, length: usize) -> Option<Vec<&Experience>> {
        if self.size < length {
            return None;
        }

        use rand::Rng;
        let mut rng = rand::thread_rng();

        let start = rng.gen_range(0..(self.size - length));

        let sequence: Vec<&Experience> = (start..(start + length))
            .filter_map(|i| self.get(i))
            .collect();

        if sequence.len() == length {
            Some(sequence)
        } else {
            None
        }
    }

    /// Get total reward in buffer
    pub fn total_reward(&self) -> f64 {
        (0..self.size)
            .filter_map(|i| self.get(i))
            .map(|e| e.reward)
            .sum()
    }

    /// Get average reward
    pub fn average_reward(&self) -> f64 {
        if self.size == 0 {
            return 0.0;
        }
        self.total_reward() / self.size as f64
    }

    /// Number of experiences stored
    pub fn len(&self) -> usize {
        self.size
    }

    /// Check if empty
    pub fn is_empty(&self) -> bool {
        self.size == 0
    }

    /// Check if full
    pub fn is_full(&self) -> bool {
        self.size >= self.capacity
    }

    /// Clear buffer
    pub fn clear(&mut self) {
        self.buffer.clear();
        self.write_pos = 0;
        self.size = 0;
        self.episode_starts = vec![0];
    }

    /// Get capacity
    pub fn capacity(&self) -> usize {
        self.capacity
    }

    /// Get number of episodes
    pub fn num_episodes(&self) -> usize {
        self.episode_starts.len()
    }

    /// Iterate over all experiences
    pub fn iter(&self) -> impl Iterator<Item = &Experience> {
        (0..self.size).filter_map(move |i| self.get(i))
    }
}

/// Prioritized experience replay
#[derive(Debug, Clone)]
pub struct PrioritizedBuffer {
    /// Base buffer
    buffer: ExperienceBuffer,

    /// Priority values
    priorities: Vec<f64>,

    /// Alpha for prioritization
    alpha: f64,

    /// Beta for importance sampling
    beta: f64,
}

impl PrioritizedBuffer {
    /// Create new prioritized buffer
    pub fn new(capacity: usize, alpha: f64, beta: f64) -> Self {
        Self {
            buffer: ExperienceBuffer::new(capacity),
            priorities: Vec::with_capacity(capacity),
            alpha,
            beta,
        }
    }

    /// Push with priority
    pub fn push(&mut self, experience: Experience, priority: f64) {
        self.buffer.push(experience);

        if self.priorities.len() < self.buffer.capacity() {
            self.priorities.push(priority.powf(self.alpha));
        } else {
            let idx = (self.buffer.write_pos + self.buffer.capacity() - 1) % self.buffer.capacity();
            self.priorities[idx] = priority.powf(self.alpha);
        }
    }

    /// Sample based on priorities
    pub fn sample(&self, batch_size: usize) -> Vec<(usize, &Experience, f64)> {
        use rand::distributions::{Distribution, WeightedIndex};
        let mut rng = rand::thread_rng();

        let total: f64 = self.priorities.iter().take(self.buffer.size).sum();
        if total <= 0.0 {
            return Vec::new();
        }

        let weights: Vec<f64> = self.priorities.iter().take(self.buffer.size).copied().collect();
        let dist = match WeightedIndex::new(&weights) {
            Ok(d) => d,
            Err(_) => return Vec::new(),
        };

        let mut samples = Vec::new();
        for _ in 0..batch_size {
            let idx = dist.sample(&mut rng);
            if let Some(exp) = self.buffer.get(idx) {
                // Importance sampling weight
                let prob = weights[idx] / total;
                let weight = (self.buffer.size as f64 * prob).powf(-self.beta);
                samples.push((idx, exp, weight));
            }
        }

        // Normalize weights
        let max_weight = samples.iter().map(|&(_, _, w)| w).fold(0.0, f64::max);
        samples
            .into_iter()
            .map(|(i, e, w)| (i, e, w / max_weight))
            .collect()
    }

    /// Update priority for an experience
    pub fn update_priority(&mut self, index: usize, priority: f64) {
        if index < self.priorities.len() {
            self.priorities[index] = priority.powf(self.alpha);
        }
    }

    /// Get buffer size
    pub fn len(&self) -> usize {
        self.buffer.len()
    }

    /// Check if empty
    pub fn is_empty(&self) -> bool {
        self.buffer.is_empty()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_experience_buffer() {
        let mut buffer = ExperienceBuffer::new(10);

        for i in 0..5 {
            buffer.push(Experience::new(i as f64, vec![i as f32], vec![]));
        }

        assert_eq!(buffer.len(), 5);
        assert!(!buffer.is_full());

        let exp = buffer.get(0).unwrap();
        assert_eq!(exp.time, 0.0);
    }

    #[test]
    fn test_buffer_overflow() {
        let mut buffer = ExperienceBuffer::new(5);

        for i in 0..10 {
            buffer.push(Experience::new(i as f64, vec![i as f32], vec![]));
        }

        assert_eq!(buffer.len(), 5);
        assert!(buffer.is_full());

        // Should have experiences 5-9
        let exp = buffer.get(0).unwrap();
        assert_eq!(exp.time, 5.0);
    }

    #[test]
    fn test_sampling() {
        let mut buffer = ExperienceBuffer::new(100);

        for i in 0..50 {
            buffer.push(Experience::new(i as f64, vec![i as f32], vec![]));
        }

        let samples = buffer.sample(10);
        assert_eq!(samples.len(), 10);
    }
}
