//! Time Synchronization
//!
//! Handles multi-rate integration between neural and physics simulations.

#[cfg(feature = "serde")]
use serde::{Deserialize, Serialize};

/// Integration strategy for multi-rate simulation
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub enum IntegrationStrategy {
    /// Run neural steps, then physics step
    Sequential,

    /// Interleave neural and physics with interpolation
    Interleaved,

    /// Adaptive stepping based on activity
    Adaptive,
}

/// Time synchronization controller
#[derive(Debug, Clone)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct TimeSync {
    /// Neural simulation time step (ms)
    pub neural_dt: f32,

    /// Physics simulation time step (ms)
    pub physics_dt: f32,

    /// Integration strategy
    pub strategy: IntegrationStrategy,

    /// Current neural time (ms)
    neural_time: f64,

    /// Current physics time (ms)
    physics_time: f64,

    /// Neural steps per physics step (computed)
    steps_per_physics: u32,

    /// Accumulated neural time since last physics step
    accumulated_neural: f64,
}

impl Default for TimeSync {
    fn default() -> Self {
        Self::new(0.025, 0.5) // 0.025 ms neural, 0.5 ms physics
    }
}

impl TimeSync {
    /// Create new time synchronizer
    pub fn new(neural_dt: f32, physics_dt: f32) -> Self {
        let steps = (physics_dt / neural_dt).round() as u32;
        let steps = steps.max(1);

        Self {
            neural_dt,
            physics_dt,
            strategy: IntegrationStrategy::Sequential,
            neural_time: 0.0,
            physics_time: 0.0,
            steps_per_physics: steps,
            accumulated_neural: 0.0,
        }
    }

    /// Create for specific model level
    pub fn for_model_level(level: hyperphysics_connectome::ModelLevel, physics_dt: f32) -> Self {
        let neural_dt = level.recommended_dt();
        Self::new(neural_dt, physics_dt)
    }

    /// Get number of neural steps per physics step
    pub fn neural_steps_per_physics(&self) -> u32 {
        self.steps_per_physics
    }

    /// Get current neural time
    pub fn neural_time(&self) -> f64 {
        self.neural_time
    }

    /// Get current physics time
    pub fn physics_time(&self) -> f64 {
        self.physics_time
    }

    /// Advance neural time by one step
    pub fn advance_neural(&mut self) {
        self.neural_time += self.neural_dt as f64;
        self.accumulated_neural += self.neural_dt as f64;
    }

    /// Advance physics time by one step
    pub fn advance_physics(&mut self) {
        self.physics_time += self.physics_dt as f64;
    }

    /// Check if physics step should occur
    pub fn should_step_physics(&self) -> bool {
        self.accumulated_neural >= self.physics_dt as f64
    }

    /// Reset accumulated neural time (after physics step)
    pub fn reset_accumulated(&mut self) {
        self.accumulated_neural = 0.0;
    }

    /// Synchronize times (call after physics step)
    pub fn synchronize(&mut self) {
        // Align neural time to physics time
        self.neural_time = self.physics_time;
        self.accumulated_neural = 0.0;
    }

    /// Reset to initial state
    pub fn reset(&mut self) {
        self.neural_time = 0.0;
        self.physics_time = 0.0;
        self.accumulated_neural = 0.0;
    }

    /// Get interpolation factor for smooth muscle activation
    /// Returns value between 0 and 1 indicating position within physics step
    pub fn interpolation_factor(&self) -> f32 {
        (self.accumulated_neural / self.physics_dt as f64) as f32
    }

    /// Compute number of neural steps needed to reach next physics time
    pub fn steps_until_physics(&self) -> u32 {
        let remaining = self.physics_dt as f64 - self.accumulated_neural;
        (remaining / self.neural_dt as f64).ceil() as u32
    }
}

/// Schedule for one integration cycle
#[derive(Debug, Clone)]
pub struct IntegrationSchedule {
    /// Number of neural steps before physics
    pub neural_steps_before: u32,

    /// Run physics step
    pub physics_step: bool,

    /// Number of neural steps after physics
    pub neural_steps_after: u32,

    /// Whether to update proprioception
    pub update_proprioception: bool,
}

impl Default for IntegrationSchedule {
    fn default() -> Self {
        Self {
            neural_steps_before: 10,
            physics_step: true,
            neural_steps_after: 0,
            update_proprioception: true,
        }
    }
}

/// Compute integration schedule for one physics step
pub fn compute_schedule(sync: &TimeSync) -> IntegrationSchedule {
    match sync.strategy {
        IntegrationStrategy::Sequential => IntegrationSchedule {
            neural_steps_before: sync.steps_per_physics,
            physics_step: true,
            neural_steps_after: 0,
            update_proprioception: true,
        },

        IntegrationStrategy::Interleaved => {
            // Half before, half after
            let half = sync.steps_per_physics / 2;
            IntegrationSchedule {
                neural_steps_before: half,
                physics_step: true,
                neural_steps_after: sync.steps_per_physics - half,
                update_proprioception: true,
            }
        }

        IntegrationStrategy::Adaptive => {
            // For now, same as sequential
            // In full implementation, would adapt based on neural activity
            IntegrationSchedule {
                neural_steps_before: sync.steps_per_physics,
                physics_step: true,
                neural_steps_after: 0,
                update_proprioception: true,
            }
        }
    }
}

/// Multi-rate integration state machine
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum IntegrationPhase {
    /// Running neural simulation
    Neural,

    /// Transferring muscle activation
    TransferToBody,

    /// Running physics simulation
    Physics,

    /// Transferring proprioception
    TransferToNeural,

    /// Cycle complete
    Complete,
}

/// Integration cycle controller
#[derive(Debug, Clone)]
pub struct IntegrationCycle {
    /// Current phase
    pub phase: IntegrationPhase,

    /// Neural steps completed in current cycle
    neural_steps_done: u32,

    /// Target neural steps for current cycle
    neural_steps_target: u32,

    /// Schedule for this cycle
    schedule: IntegrationSchedule,
}

impl IntegrationCycle {
    /// Start new integration cycle
    pub fn new(schedule: IntegrationSchedule) -> Self {
        Self {
            phase: IntegrationPhase::Neural,
            neural_steps_done: 0,
            neural_steps_target: schedule.neural_steps_before,
            schedule,
        }
    }

    /// Advance to next step in cycle
    pub fn advance(&mut self) -> IntegrationPhase {
        match self.phase {
            IntegrationPhase::Neural => {
                self.neural_steps_done += 1;
                if self.neural_steps_done >= self.neural_steps_target {
                    if self.schedule.physics_step && self.neural_steps_target == self.schedule.neural_steps_before {
                        self.phase = IntegrationPhase::TransferToBody;
                    } else {
                        self.phase = IntegrationPhase::Complete;
                    }
                }
            }

            IntegrationPhase::TransferToBody => {
                self.phase = IntegrationPhase::Physics;
            }

            IntegrationPhase::Physics => {
                if self.schedule.update_proprioception {
                    self.phase = IntegrationPhase::TransferToNeural;
                } else if self.schedule.neural_steps_after > 0 {
                    self.phase = IntegrationPhase::Neural;
                    self.neural_steps_done = 0;
                    self.neural_steps_target = self.schedule.neural_steps_after;
                } else {
                    self.phase = IntegrationPhase::Complete;
                }
            }

            IntegrationPhase::TransferToNeural => {
                if self.schedule.neural_steps_after > 0 {
                    self.phase = IntegrationPhase::Neural;
                    self.neural_steps_done = 0;
                    self.neural_steps_target = self.schedule.neural_steps_after;
                } else {
                    self.phase = IntegrationPhase::Complete;
                }
            }

            IntegrationPhase::Complete => {
                // Stay in complete
            }
        }

        self.phase
    }

    /// Check if cycle is complete
    pub fn is_complete(&self) -> bool {
        self.phase == IntegrationPhase::Complete
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_time_sync() {
        let sync = TimeSync::new(0.025, 0.5);
        assert_eq!(sync.neural_steps_per_physics(), 20);
    }

    #[test]
    fn test_integration_cycle() {
        let schedule = IntegrationSchedule {
            neural_steps_before: 3,
            physics_step: true,
            neural_steps_after: 0,
            update_proprioception: true,
        };

        let mut cycle = IntegrationCycle::new(schedule);

        // 3 neural steps
        assert_eq!(cycle.phase, IntegrationPhase::Neural);
        cycle.advance();
        cycle.advance();
        cycle.advance();

        assert_eq!(cycle.phase, IntegrationPhase::TransferToBody);
        cycle.advance();

        assert_eq!(cycle.phase, IntegrationPhase::Physics);
        cycle.advance();

        assert_eq!(cycle.phase, IntegrationPhase::TransferToNeural);
        cycle.advance();

        assert_eq!(cycle.phase, IntegrationPhase::Complete);
    }
}
