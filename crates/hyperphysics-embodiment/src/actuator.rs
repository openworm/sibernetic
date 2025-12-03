//! Muscle Actuator Model
//!
//! Converts neural muscle activation signals to physical forces
//! applied to SPH particles.

#[cfg(feature = "serde")]
use serde::{Deserialize, Serialize};

use hyperphysics_sph::MuscleActivation;

/// Actuator configuration
#[derive(Debug, Clone)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct ActuatorConfig {
    /// Maximum muscle force (N)
    pub max_force: f32,

    /// Force rise time constant (ms)
    pub tau_rise: f32,

    /// Force decay time constant (ms)
    pub tau_decay: f32,

    /// Activation to force nonlinearity exponent
    pub force_exponent: f32,

    /// Fatigue time constant (ms) - longer = slower fatigue
    pub tau_fatigue: f32,

    /// Recovery time constant (ms)
    pub tau_recovery: f32,

    /// Enable fatigue modeling
    pub fatigue_enabled: bool,
}

impl Default for ActuatorConfig {
    fn default() -> Self {
        Self {
            max_force: 10.0,
            tau_rise: 10.0,
            tau_decay: 30.0,
            force_exponent: 2.0,
            tau_fatigue: 5000.0,
            tau_recovery: 2000.0,
            fatigue_enabled: false,
        }
    }
}

impl ActuatorConfig {
    /// Fast twitch muscles
    pub fn fast_twitch() -> Self {
        Self {
            max_force: 15.0,
            tau_rise: 5.0,
            tau_decay: 15.0,
            force_exponent: 2.0,
            tau_fatigue: 2000.0,
            tau_recovery: 1000.0,
            fatigue_enabled: true,
        }
    }

    /// Slow twitch muscles (more like C. elegans body wall)
    pub fn slow_twitch() -> Self {
        Self {
            max_force: 8.0,
            tau_rise: 20.0,
            tau_decay: 50.0,
            force_exponent: 1.5,
            tau_fatigue: 10000.0,
            tau_recovery: 3000.0,
            fatigue_enabled: true,
        }
    }
}

/// Single muscle actuator state
#[derive(Debug, Clone, Copy, Default)]
pub struct ActuatorState {
    /// Current force output (0-1, normalized)
    pub force: f32,

    /// Current activation level
    pub activation: f32,

    /// Fatigue level (0 = fresh, 1 = fully fatigued)
    pub fatigue: f32,

    /// Internal calcium-like activation variable
    calcium: f32,
}

impl ActuatorState {
    /// Get effective force (accounting for fatigue)
    pub fn effective_force(&self) -> f32 {
        self.force * (1.0 - self.fatigue * 0.5)
    }
}

/// Muscle actuator system for the whole body
#[derive(Debug, Clone)]
pub struct Actuator {
    /// Configuration
    config: ActuatorConfig,

    /// State for each muscle (96 for C. elegans)
    states: Vec<ActuatorState>,

    /// Smoothed activation buffer
    smoothed_activation: Vec<f32>,
}

impl Actuator {
    /// Create actuator system for given number of muscles
    pub fn new(num_muscles: usize, config: ActuatorConfig) -> Self {
        Self {
            config,
            states: vec![ActuatorState::default(); num_muscles],
            smoothed_activation: vec![0.0; num_muscles],
        }
    }

    /// Create for C. elegans (96 muscles)
    pub fn celegans() -> Self {
        Self::new(96, ActuatorConfig::slow_twitch())
    }

    /// Update actuators with new neural activation
    pub fn update(&mut self, neural_activation: &[f32; 96], dt: f32) {
        for (i, &activation) in neural_activation.iter().enumerate() {
            if i >= self.states.len() {
                break;
            }

            // Smooth activation (low-pass filter)
            let alpha = if activation > self.smoothed_activation[i] {
                dt / (self.config.tau_rise + dt)
            } else {
                dt / (self.config.tau_decay + dt)
            };
            self.smoothed_activation[i] += alpha * (activation - self.smoothed_activation[i]);

            let state = &mut self.states[i];
            state.activation = self.smoothed_activation[i];

            // Calcium dynamics (simplified)
            let target_calcium = state.activation;
            state.calcium += (target_calcium - state.calcium) * dt / 5.0;

            // Force from calcium with nonlinearity
            state.force = state.calcium.powf(self.config.force_exponent);
            state.force = state.force.clamp(0.0, 1.0);

            // Fatigue dynamics
            if self.config.fatigue_enabled {
                if state.force > 0.1 {
                    // Fatigue accumulates during activity
                    state.fatigue += state.force * dt / self.config.tau_fatigue;
                } else {
                    // Recovery during rest
                    state.fatigue -= state.fatigue * dt / self.config.tau_recovery;
                }
                state.fatigue = state.fatigue.clamp(0.0, 1.0);
            }
        }
    }

    /// Get force output for all muscles
    pub fn get_forces(&self) -> [f32; 96] {
        let mut forces = [0.0_f32; 96];
        for (i, state) in self.states.iter().enumerate() {
            if i < 96 {
                forces[i] = state.effective_force() * self.config.max_force;
            }
        }
        forces
    }

    /// Get normalized activation suitable for SPH
    pub fn get_activation(&self) -> MuscleActivation {
        let mut activation = MuscleActivation::new();
        for (i, state) in self.states.iter().enumerate() {
            if i < 96 {
                activation.set_by_index(i, state.effective_force());
            }
        }
        activation
    }

    /// Get state of specific muscle
    pub fn get_muscle_state(&self, index: usize) -> Option<&ActuatorState> {
        self.states.get(index)
    }

    /// Reset all actuators
    pub fn reset(&mut self) {
        for state in &mut self.states {
            *state = ActuatorState::default();
        }
        self.smoothed_activation.fill(0.0);
    }

    /// Get total energy expenditure (arbitrary units)
    pub fn energy_expenditure(&self) -> f32 {
        self.states.iter().map(|s| s.force * s.force).sum()
    }

    /// Get fatigue level (average across all muscles)
    pub fn average_fatigue(&self) -> f32 {
        if self.states.is_empty() {
            return 0.0;
        }
        let sum: f32 = self.states.iter().map(|s| s.fatigue).sum();
        sum / self.states.len() as f32
    }
}

/// Convert neural activation pattern to muscle forces
pub fn neural_to_muscle_forces(
    neural_activation: &[f32; 96],
    actuator: &mut Actuator,
    dt: f32,
) -> MuscleActivation {
    actuator.update(neural_activation, dt);
    actuator.get_activation()
}

/// Compute antagonist activation pattern
/// (for smooth locomotion, opposing muscles should be anti-phase)
pub fn compute_antagonist_pattern(activation: &[f32; 96]) -> [f32; 96] {
    let mut result = [0.0_f32; 96];

    for seg in 0..24 {
        // Dorsal and ventral are antagonists
        let dr = activation[seg * 4];     // Dorsal right
        let vr = activation[seg * 4 + 1]; // Ventral right
        let vl = activation[seg * 4 + 2]; // Ventral left
        let dl = activation[seg * 4 + 3]; // Dorsal left

        // Anti-phase pattern
        result[seg * 4] = vr.max(vl) * 0.5;     // DR gets anti-ventral
        result[seg * 4 + 1] = dr.max(dl) * 0.5; // VR gets anti-dorsal
        result[seg * 4 + 2] = dr.max(dl) * 0.5; // VL gets anti-dorsal
        result[seg * 4 + 3] = vr.max(vl) * 0.5; // DL gets anti-ventral
    }

    result
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_actuator_creation() {
        let actuator = Actuator::celegans();
        assert_eq!(actuator.states.len(), 96);
    }

    #[test]
    fn test_actuator_update() {
        let mut actuator = Actuator::new(96, ActuatorConfig::default());

        let mut activation = [0.0_f32; 96];
        activation[0] = 1.0;

        // Run for a few steps
        for _ in 0..100 {
            actuator.update(&activation, 0.1);
        }

        // Muscle 0 should have significant force
        assert!(actuator.states[0].force > 0.5);

        // Other muscles should be near zero
        assert!(actuator.states[1].force < 0.01);
    }

    #[test]
    fn test_fatigue() {
        let config = ActuatorConfig {
            fatigue_enabled: true,
            tau_fatigue: 100.0, // Fast fatigue for testing
            ..Default::default()
        };
        let mut actuator = Actuator::new(1, config);

        let activation = [1.0_f32; 96];

        // Run until fatigued
        for _ in 0..1000 {
            actuator.update(&activation, 0.1);
        }

        assert!(actuator.states[0].fatigue > 0.3);
    }
}
