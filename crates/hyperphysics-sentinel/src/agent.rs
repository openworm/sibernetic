//! Sentinel Agent
//!
//! Core agent implementation combining neural network, body, and learning.

#[cfg(feature = "serde")]
use serde::{Deserialize, Serialize};

use hyperphysics_connectome::{ModelLevel, SpikingNetwork};
use hyperphysics_embodiment::{EmbodiedWorm, CouplingConfig};
use hyperphysics_stdp::{PlasticityController, RewardModulatedParams, RewardModulatedStdp, RewardSignal};
use hyperphysics_nas::{Genome, InnovationTracker};

use crate::consciousness::ConsciousnessMetrics;
use crate::experience::{Experience, ExperienceBuffer, Reward};
use crate::lifecycle::LifecycleStage;
use crate::{AgentId, Result, SimTime};

/// Agent state
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub enum AgentState {
    /// Just spawned, initializing
    Initializing,
    /// Active and running
    Active,
    /// Learning/training mode
    Learning,
    /// Evaluating fitness
    Evaluating,
    /// Reproducing/spawning offspring
    Reproducing,
    /// Paused
    Paused,
    /// Terminated
    Dead,
}

impl Default for AgentState {
    fn default() -> Self {
        Self::Initializing
    }
}

/// Sentinel configuration
#[derive(Debug, Clone)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct SentinelConfig {
    /// Model complexity level
    pub model_level: ModelLevel,

    /// Enable online learning
    pub learning_enabled: bool,

    /// Learning rate
    pub learning_rate: f32,

    /// Experience buffer size
    pub experience_buffer_size: usize,

    /// Target consciousness level (Φ)
    pub target_phi: f64,

    /// Enable consciousness monitoring
    pub consciousness_monitoring: bool,

    /// Reward discount factor (γ)
    pub gamma: f32,

    /// Maximum age (simulation time)
    pub max_age: Option<SimTime>,

    /// Neural coupling configuration
    pub coupling: CouplingConfig,
}

impl Default for SentinelConfig {
    fn default() -> Self {
        Self {
            model_level: ModelLevel::B,
            learning_enabled: true,
            learning_rate: 0.01,
            experience_buffer_size: 10000,
            target_phi: 0.1,
            consciousness_monitoring: true,
            gamma: 0.99,
            max_age: None,
            coupling: CouplingConfig::default(),
        }
    }
}

impl SentinelConfig {
    /// High-fidelity configuration
    pub fn high_fidelity() -> Self {
        Self {
            model_level: ModelLevel::C,
            learning_enabled: true,
            learning_rate: 0.005,
            consciousness_monitoring: true,
            coupling: CouplingConfig::high_fidelity(),
            ..Default::default()
        }
    }

    /// Fast/simple configuration
    pub fn fast() -> Self {
        Self {
            model_level: ModelLevel::A,
            learning_enabled: false,
            consciousness_monitoring: false,
            coupling: CouplingConfig::open_loop(),
            ..Default::default()
        }
    }

    /// Trading-optimized configuration
    pub fn trading() -> Self {
        Self {
            model_level: ModelLevel::A,
            learning_enabled: true,
            learning_rate: 0.1, // Fast learning
            experience_buffer_size: 1000,
            gamma: 0.95,
            consciousness_monitoring: false, // Speed over monitoring
            coupling: CouplingConfig::open_loop(),
            ..Default::default()
        }
    }
}

/// A Sentinel agent
pub struct Sentinel {
    /// Unique identifier
    id: AgentId,

    /// Configuration
    config: SentinelConfig,

    /// Current state
    state: AgentState,

    /// Lifecycle stage
    lifecycle: LifecycleStage,

    /// Embodied simulation
    embodiment: EmbodiedWorm,

    /// Plasticity controller
    plasticity: PlasticityController,

    /// Experience buffer
    experiences: ExperienceBuffer,

    /// Consciousness metrics
    consciousness: ConsciousnessMetrics,

    /// Genome (for evolution)
    genome: Option<Genome>,

    /// Birth time
    birth_time: SimTime,

    /// Current simulation time
    current_time: SimTime,

    /// Total reward accumulated
    total_reward: f64,

    /// Step count
    step_count: u64,

    /// Parent IDs (for lineage tracking)
    parent_ids: Vec<AgentId>,

    /// Generation number
    generation: u32,
}

impl Sentinel {
    /// Create new sentinel agent
    pub fn new(id: AgentId, config: SentinelConfig) -> Result<Self> {
        let embodiment = EmbodiedWorm::new(config.model_level, config.coupling.clone());

        let num_synapses = embodiment.neural().num_synapses();
        let plasticity = if config.learning_enabled {
            PlasticityController::with_reward_stdp(num_synapses)
        } else {
            PlasticityController::new()
        };

        Ok(Self {
            id,
            config: config.clone(),
            state: AgentState::Initializing,
            lifecycle: LifecycleStage::Embryo,
            embodiment,
            plasticity,
            experiences: ExperienceBuffer::new(config.experience_buffer_size),
            consciousness: ConsciousnessMetrics::new(),
            genome: None,
            birth_time: 0.0,
            current_time: 0.0,
            total_reward: 0.0,
            step_count: 0,
            parent_ids: Vec::new(),
            generation: 0,
        })
    }

    /// Create from genome
    pub fn from_genome(id: AgentId, genome: Genome, config: SentinelConfig) -> Result<Self> {
        let mut agent = Self::new(id, config)?;
        agent.genome = Some(genome);
        Ok(agent)
    }

    /// Initialize the agent
    pub fn initialize(&mut self) -> Result<()> {
        self.embodiment.initialize_body()?;
        self.birth_time = self.current_time;
        self.state = AgentState::Active;
        self.lifecycle = LifecycleStage::Juvenile;
        Ok(())
    }

    /// Perform one simulation step
    pub fn step(&mut self, dt: f32) {
        if self.state != AgentState::Active && self.state != AgentState::Learning {
            return;
        }

        // Step embodied simulation
        self.embodiment.step();

        // Record experience
        let state = self.capture_state();
        self.experiences.push(Experience {
            time: self.current_time,
            state,
            action: self.embodiment.muscle_activations().to_vec(),
            reward: 0.0, // Will be set by reward signal
            done: false,
        });

        // Update consciousness metrics periodically
        if self.config.consciousness_monitoring && self.step_count % 100 == 0 {
            self.update_consciousness_metrics();
        }

        // Update time
        self.current_time += dt as f64;
        self.step_count += 1;

        // Check age limit
        if let Some(max_age) = self.config.max_age {
            if self.age() > max_age {
                self.state = AgentState::Dead;
            }
        }
    }

    /// Apply reward signal
    pub fn receive_reward(&mut self, reward: Reward) {
        self.total_reward += reward.value;

        // Update last experience
        if let Some(exp) = self.experiences.last_mut() {
            exp.reward = reward.value;
        }

        // Apply to plasticity
        if self.config.learning_enabled {
            let signal = RewardSignal::reward(reward.value as f32, self.current_time);
            // Would need to access the reward STDP rule here
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

    /// Capture current state as observation
    fn capture_state(&self) -> Vec<f32> {
        let mut state = Vec::new();

        // Neural state (voltage summary)
        let voltages = self.embodiment.neural().get_voltages();
        let avg_v: f32 = voltages.iter().sum::<f32>() / voltages.len().max(1) as f32;
        state.push(avg_v);

        // Body state (center of mass, velocity)
        let body_state = self.embodiment.state();
        state.push(body_state.center_of_mass[0] as f32);
        state.push(body_state.center_of_mass[1] as f32);
        state.push(body_state.center_of_mass[2] as f32);

        // Energy
        state.push(body_state.kinetic_energy as f32);

        state
    }

    /// Update consciousness metrics
    fn update_consciousness_metrics(&mut self) {
        // Simplified consciousness calculation
        // In full implementation, would compute Φ, causal density, etc.
        let voltages = self.embodiment.neural().get_voltages();
        let activity: f32 = voltages.iter().map(|v| (v + 65.0).max(0.0)).sum();
        let normalized_activity = activity / (voltages.len() as f32 * 100.0);

        self.consciousness.activity_level = normalized_activity as f64;

        // Estimate Φ from activity patterns (simplified)
        self.consciousness.phi = normalized_activity as f64 * 0.5;
    }

    // ========== Getters ==========

    /// Get agent ID
    pub fn id(&self) -> AgentId {
        self.id
    }

    /// Get current state
    pub fn state(&self) -> AgentState {
        self.state
    }

    /// Set state
    pub fn set_state(&mut self, state: AgentState) {
        self.state = state;
    }

    /// Get lifecycle stage
    pub fn lifecycle(&self) -> LifecycleStage {
        self.lifecycle
    }

    /// Set lifecycle stage
    pub fn set_lifecycle(&mut self, stage: LifecycleStage) {
        self.lifecycle = stage;
    }

    /// Get current simulation time
    pub fn time(&self) -> SimTime {
        self.current_time
    }

    /// Get age (time since birth)
    pub fn age(&self) -> SimTime {
        self.current_time - self.birth_time
    }

    /// Get total reward
    pub fn total_reward(&self) -> f64 {
        self.total_reward
    }

    /// Get step count
    pub fn step_count(&self) -> u64 {
        self.step_count
    }

    /// Get consciousness metrics
    pub fn consciousness(&self) -> &ConsciousnessMetrics {
        &self.consciousness
    }

    /// Get embodiment reference
    pub fn embodiment(&self) -> &EmbodiedWorm {
        &self.embodiment
    }

    /// Get mutable embodiment reference
    pub fn embodiment_mut(&mut self) -> &mut EmbodiedWorm {
        &mut self.embodiment
    }

    /// Get genome
    pub fn genome(&self) -> Option<&Genome> {
        self.genome.as_ref()
    }

    /// Set genome
    pub fn set_genome(&mut self, genome: Genome) {
        self.genome = Some(genome);
    }

    /// Get generation
    pub fn generation(&self) -> u32 {
        self.generation
    }

    /// Set generation
    pub fn set_generation(&mut self, gen: u32) {
        self.generation = gen;
    }

    /// Get parent IDs
    pub fn parents(&self) -> &[AgentId] {
        &self.parent_ids
    }

    /// Add parent ID
    pub fn add_parent(&mut self, parent_id: AgentId) {
        self.parent_ids.push(parent_id);
    }

    /// Get experience buffer
    pub fn experiences(&self) -> &ExperienceBuffer {
        &self.experiences
    }

    /// Get fitness (for evolution)
    pub fn fitness(&self) -> f32 {
        // Combine reward with consciousness
        let reward_component = self.total_reward as f32;
        let phi_component = self.consciousness.phi as f32;
        let efficiency = self.step_count as f32 / (self.total_reward.abs() as f32 + 1.0);

        reward_component + phi_component * 10.0 - efficiency * 0.01
    }

    /// Check if alive
    pub fn is_alive(&self) -> bool {
        self.state != AgentState::Dead
    }

    /// Kill the agent
    pub fn kill(&mut self) {
        self.state = AgentState::Dead;
        self.lifecycle = LifecycleStage::Dead;
    }

    /// Reset agent state
    pub fn reset(&mut self) {
        self.embodiment.reset();
        self.experiences.clear();
        self.total_reward = 0.0;
        self.step_count = 0;
        self.current_time = self.birth_time;
        self.state = AgentState::Active;
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_agent_creation() {
        let agent = Sentinel::new(1, SentinelConfig::fast()).unwrap();
        assert_eq!(agent.id(), 1);
        assert_eq!(agent.state(), AgentState::Initializing);
    }

    #[test]
    fn test_agent_initialization() {
        let mut agent = Sentinel::new(1, SentinelConfig::fast()).unwrap();
        agent.initialize().unwrap();
        assert_eq!(agent.state(), AgentState::Active);
    }
}
