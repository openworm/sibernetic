//! Main Worm Simulation
//!
//! Unified simulation controller for C. elegans.

use crate::config::SimulationConfig;
use crate::{Result, SimTime, SimulationError};

use hyperphysics_connectome::ModelLevel;
use hyperphysics_embodiment::EmbodiedWorm;
use hyperphysics_sph::MuscleActivation;

/// Simulation mode
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SimulationMode {
    /// Standard embodied simulation
    Embodied,
    /// Neural network only (no physics)
    NeuralOnly,
    /// Physics only (no neural)
    PhysicsOnly,
    /// Playback mode (replay recorded data)
    Playback,
}

/// Main worm simulation
pub struct WormSimulation {
    /// Configuration
    config: SimulationConfig,

    /// Simulation mode
    mode: SimulationMode,

    /// Embodied worm (neural + physics)
    embodiment: EmbodiedWorm,

    /// Current simulation time (ms)
    time: SimTime,

    /// Step counter
    step_count: u64,

    /// Recording buffer
    recordings: SimulationRecordings,

    /// Paused state
    paused: bool,

    /// Initialized flag
    initialized: bool,
}

impl WormSimulation {
    /// Create new simulation
    pub fn new(config: SimulationConfig) -> Self {
        let embodiment = EmbodiedWorm::new(config.model_level, config.coupling_config.clone());

        Self {
            config,
            mode: SimulationMode::Embodied,
            embodiment,
            time: 0.0,
            step_count: 0,
            recordings: SimulationRecordings::new(),
            paused: false,
            initialized: false,
        }
    }

    /// Create with preset configuration
    pub fn fast() -> Self {
        Self::new(SimulationConfig::fast())
    }

    /// Create high-fidelity simulation
    pub fn high_fidelity() -> Self {
        Self::new(SimulationConfig::high_fidelity())
    }

    /// Create demo simulation
    pub fn demo() -> Self {
        Self::new(SimulationConfig::demo())
    }

    /// Initialize simulation
    pub fn initialize(&mut self) -> Result<()> {
        // Validate configuration
        self.config.validate()
            .map_err(|e| SimulationError::Configuration(e))?;

        // Initialize embodied worm
        self.embodiment.initialize_body()?;

        self.initialized = true;
        Ok(())
    }

    /// Set simulation mode
    pub fn set_mode(&mut self, mode: SimulationMode) {
        self.mode = mode;
    }

    /// Perform one simulation step
    pub fn step(&mut self) {
        if self.paused || !self.initialized {
            return;
        }

        match self.mode {
            SimulationMode::Embodied => {
                self.embodiment.step();
            }
            SimulationMode::NeuralOnly => {
                // Step neural network without physics
                self.embodiment.neural_mut().step(self.config.dt);
            }
            SimulationMode::PhysicsOnly => {
                // Step physics without neural
                self.embodiment.body_mut().step();
            }
            SimulationMode::Playback => {
                // Playback mode - advance time only
            }
        }

        // Record data
        self.record();

        // Update time
        self.time += self.config.dt as f64;
        self.step_count += 1;
    }

    /// Run simulation for specified duration
    pub fn run(&mut self, duration_ms: SimTime) {
        let steps = (duration_ms / self.config.dt as f64).ceil() as u64;
        for _ in 0..steps {
            self.step();

            // Check max time
            if let Some(max_time) = self.config.max_time {
                if self.time >= max_time {
                    break;
                }
            }
        }
    }

    /// Run until condition is met
    pub fn run_until<F: Fn(&Self) -> bool>(&mut self, condition: F) {
        while !condition(self) {
            self.step();

            // Safety limit
            if self.step_count > 10_000_000 {
                break;
            }
        }
    }

    /// Apply stimulus to sensory neurons
    pub fn stimulate(&mut self, neuron_name: &str, current: f32) {
        self.embodiment.stimulate(neuron_name, current);
    }

    /// Apply touch stimulus
    pub fn touch(&mut self, position: [f32; 3], strength: f32) {
        self.embodiment.touch(position, strength);
    }

    /// Clear all stimuli
    pub fn clear_stimuli(&mut self) {
        self.embodiment.neural_mut().clear_inputs();
    }

    /// Record current state
    fn record(&mut self) {
        if self.step_count % self.config.output.recording_interval as u64 != 0 {
            return;
        }

        if self.config.output.record_spikes {
            let spikes = self.embodiment.neural().get_spikes();
            if !spikes.is_empty() {
                self.recordings.spikes.push((self.time, spikes));
            }
        }

        if self.config.output.record_muscles {
            let muscles = *self.embodiment.muscle_activations();
            self.recordings.muscles.push((self.time, muscles));
        }

        if self.config.output.record_positions {
            let com = self.embodiment.body().center_of_mass();
            self.recordings.positions.push((self.time, com));
        }

        if self.config.output.record_voltages {
            let voltages = self.embodiment.neural().get_voltages();
            self.recordings.voltages.push((self.time, voltages));
        }
    }

    // ========== Control ==========

    /// Pause simulation
    pub fn pause(&mut self) {
        self.paused = true;
    }

    /// Resume simulation
    pub fn resume(&mut self) {
        self.paused = false;
    }

    /// Check if paused
    pub fn is_paused(&self) -> bool {
        self.paused
    }

    /// Reset simulation
    pub fn reset(&mut self) {
        self.embodiment.reset();
        self.time = 0.0;
        self.step_count = 0;
        self.recordings.clear();
    }

    // ========== Getters ==========

    /// Get current time
    pub fn time(&self) -> SimTime {
        self.time
    }

    /// Get step count
    pub fn step_count(&self) -> u64 {
        self.step_count
    }

    /// Get muscle activations
    pub fn get_muscle_activations(&self) -> &[f32; 96] {
        self.embodiment.muscle_activations()
    }

    /// Get membrane voltages
    pub fn get_voltages(&self) -> Vec<f32> {
        self.embodiment.neural().get_voltages()
    }

    /// Get spike events
    pub fn get_spikes(&self) -> Vec<u16> {
        self.embodiment.neural().get_spikes()
    }

    /// Get center of mass
    pub fn get_center_of_mass(&self) -> [f64; 3] {
        self.embodiment.body().center_of_mass()
    }

    /// Get kinetic energy
    pub fn get_kinetic_energy(&self) -> f64 {
        self.embodiment.body().kinetic_energy()
    }

    /// Get number of particles
    pub fn get_num_particles(&self) -> usize {
        self.embodiment.body().num_particles()
    }

    /// Get particle positions
    pub fn get_positions(&self) -> &[f32] {
        self.embodiment.body().positions()
    }

    /// Get embodied worm reference
    pub fn embodiment(&self) -> &EmbodiedWorm {
        &self.embodiment
    }

    /// Get mutable embodied worm reference
    pub fn embodiment_mut(&mut self) -> &mut EmbodiedWorm {
        &mut self.embodiment
    }

    /// Get recordings
    pub fn recordings(&self) -> &SimulationRecordings {
        &self.recordings
    }

    /// Get configuration
    pub fn config(&self) -> &SimulationConfig {
        &self.config
    }

    /// Check if initialized
    pub fn is_initialized(&self) -> bool {
        self.initialized
    }
}

/// Recorded simulation data
#[derive(Debug, Clone, Default)]
pub struct SimulationRecordings {
    /// Spike times: (time, neuron_ids)
    pub spikes: Vec<(SimTime, Vec<u16>)>,

    /// Muscle activations: (time, 96 values)
    pub muscles: Vec<(SimTime, [f32; 96])>,

    /// Center of mass positions: (time, xyz)
    pub positions: Vec<(SimTime, [f64; 3])>,

    /// Membrane voltages: (time, all neurons)
    pub voltages: Vec<(SimTime, Vec<f32>)>,
}

impl SimulationRecordings {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn clear(&mut self) {
        self.spikes.clear();
        self.muscles.clear();
        self.positions.clear();
        self.voltages.clear();
    }

    /// Get total recording time
    pub fn duration(&self) -> SimTime {
        let times = [
            self.spikes.last().map(|(t, _)| *t),
            self.muscles.last().map(|(t, _)| *t),
            self.positions.last().map(|(t, _)| *t),
            self.voltages.last().map(|(t, _)| *t),
        ];

        times.into_iter().flatten().fold(0.0, f64::max)
    }

    /// Get total spike count
    pub fn total_spikes(&self) -> usize {
        self.spikes.iter().map(|(_, s)| s.len()).sum()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_simulation_creation() {
        let sim = WormSimulation::fast();
        assert!(!sim.is_initialized());
        assert_eq!(sim.time(), 0.0);
    }

    #[test]
    fn test_simulation_init() {
        let mut sim = WormSimulation::fast();
        assert!(sim.initialize().is_ok());
        assert!(sim.is_initialized());
    }

    #[test]
    fn test_simulation_step() {
        let mut sim = WormSimulation::fast();
        sim.initialize().unwrap();

        sim.step();
        assert!(sim.time() > 0.0);
        assert_eq!(sim.step_count(), 1);
    }

    #[test]
    fn test_simulation_run() {
        let mut sim = WormSimulation::fast();
        sim.initialize().unwrap();

        sim.run(10.0); // Run for 10 ms
        assert!(sim.time() >= 10.0);
    }
}
