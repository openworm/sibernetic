//! Agent Lifecycle Management
//!
//! Manages agent birth, growth, reproduction, and death.

#[cfg(feature = "serde")]
use serde::{Deserialize, Serialize};

use crate::agent::{Sentinel, SentinelConfig};
use crate::{AgentId, Result, SimTime};

/// Lifecycle stages
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub enum LifecycleStage {
    /// Just created, not yet initialized
    Embryo,
    /// Learning and developing
    Juvenile,
    /// Fully mature, can reproduce
    Adult,
    /// Declining capabilities
    Senescent,
    /// No longer active
    Dead,
}

impl Default for LifecycleStage {
    fn default() -> Self {
        Self::Embryo
    }
}

impl LifecycleStage {
    /// Check if can reproduce
    pub fn can_reproduce(&self) -> bool {
        matches!(self, Self::Adult)
    }

    /// Check if alive
    pub fn is_alive(&self) -> bool {
        !matches!(self, Self::Dead)
    }

    /// Get next stage
    pub fn next(&self) -> Option<Self> {
        match self {
            Self::Embryo => Some(Self::Juvenile),
            Self::Juvenile => Some(Self::Adult),
            Self::Adult => Some(Self::Senescent),
            Self::Senescent => Some(Self::Dead),
            Self::Dead => None,
        }
    }
}

/// Spawn configuration
#[derive(Debug, Clone)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct SpawnConfig {
    /// Agent configuration
    pub agent_config: SentinelConfig,

    /// Initial position (if applicable)
    pub position: Option<[f32; 3]>,

    /// Parent agent IDs
    pub parents: Vec<AgentId>,

    /// Generation number
    pub generation: u32,

    /// Mutation rate for offspring
    pub mutation_rate: f32,
}

impl Default for SpawnConfig {
    fn default() -> Self {
        Self {
            agent_config: SentinelConfig::default(),
            position: None,
            parents: Vec::new(),
            generation: 0,
            mutation_rate: 0.1,
        }
    }
}

/// Lifecycle manager
pub struct LifecycleManager {
    /// Age thresholds for stage transitions
    juvenile_age: SimTime,
    adult_age: SimTime,
    senescent_age: SimTime,
    max_age: SimTime,

    /// Reproduction cooldown
    reproduction_cooldown: SimTime,

    /// Next agent ID
    next_id: AgentId,

    /// Birth registry
    births: Vec<(AgentId, SimTime)>,

    /// Death registry
    deaths: Vec<(AgentId, SimTime)>,
}

impl Default for LifecycleManager {
    fn default() -> Self {
        Self {
            juvenile_age: 100.0,
            adult_age: 500.0,
            senescent_age: 2000.0,
            max_age: 5000.0,
            reproduction_cooldown: 200.0,
            next_id: 1,
            births: Vec::new(),
            deaths: Vec::new(),
        }
    }
}

impl LifecycleManager {
    /// Create new lifecycle manager
    pub fn new() -> Self {
        Self::default()
    }

    /// Configure age thresholds
    pub fn with_ages(
        juvenile_age: SimTime,
        adult_age: SimTime,
        senescent_age: SimTime,
        max_age: SimTime,
    ) -> Self {
        Self {
            juvenile_age,
            adult_age,
            senescent_age,
            max_age,
            ..Default::default()
        }
    }

    /// Spawn new agent
    pub fn spawn(&mut self, config: SpawnConfig, current_time: SimTime) -> Result<Sentinel> {
        let id = self.next_id;
        self.next_id += 1;

        let mut agent = Sentinel::new(id, config.agent_config)?;

        // Set lineage info
        for &parent_id in &config.parents {
            agent.add_parent(parent_id);
        }
        agent.set_generation(config.generation);

        // Register birth
        self.births.push((id, current_time));

        Ok(agent)
    }

    /// Spawn from parents (sexual reproduction)
    pub fn spawn_from_parents(
        &mut self,
        parent1: &Sentinel,
        parent2: &Sentinel,
        current_time: SimTime,
    ) -> Result<Sentinel> {
        let config = SpawnConfig {
            agent_config: SentinelConfig::default(), // Would blend parent configs
            parents: vec![parent1.id(), parent2.id()],
            generation: parent1.generation().max(parent2.generation()) + 1,
            ..Default::default()
        };

        self.spawn(config, current_time)
    }

    /// Spawn from single parent (asexual reproduction)
    pub fn spawn_from_parent(
        &mut self,
        parent: &Sentinel,
        current_time: SimTime,
    ) -> Result<Sentinel> {
        let config = SpawnConfig {
            agent_config: SentinelConfig::default(), // Would clone parent config
            parents: vec![parent.id()],
            generation: parent.generation() + 1,
            ..Default::default()
        };

        self.spawn(config, current_time)
    }

    /// Update agent lifecycle stage based on age
    pub fn update_stage(&self, agent: &mut Sentinel) {
        let age = agent.age();
        let current_stage = agent.lifecycle();

        let new_stage = if age >= self.max_age {
            LifecycleStage::Dead
        } else if age >= self.senescent_age {
            LifecycleStage::Senescent
        } else if age >= self.adult_age {
            LifecycleStage::Adult
        } else if age >= self.juvenile_age {
            LifecycleStage::Juvenile
        } else {
            LifecycleStage::Embryo
        };

        if new_stage != current_stage {
            agent.set_lifecycle(new_stage);

            if new_stage == LifecycleStage::Dead {
                agent.kill();
            }
        }
    }

    /// Check if agent can reproduce
    pub fn can_reproduce(&self, agent: &Sentinel, current_time: SimTime) -> bool {
        if !agent.lifecycle().can_reproduce() {
            return false;
        }

        // Check fitness threshold
        if agent.fitness() < 0.0 {
            return false;
        }

        true
    }

    /// Register death
    pub fn register_death(&mut self, agent_id: AgentId, time: SimTime) {
        self.deaths.push((agent_id, time));
    }

    /// Get birth count
    pub fn birth_count(&self) -> usize {
        self.births.len()
    }

    /// Get death count
    pub fn death_count(&self) -> usize {
        self.deaths.len()
    }

    /// Get next ID
    pub fn next_id(&self) -> AgentId {
        self.next_id
    }

    /// Calculate population growth rate
    pub fn growth_rate(&self, time_window: SimTime, current_time: SimTime) -> f64 {
        let start_time = (current_time - time_window).max(0.0);

        let births_in_window = self
            .births
            .iter()
            .filter(|(_, t)| *t >= start_time)
            .count();

        let deaths_in_window = self
            .deaths
            .iter()
            .filter(|(_, t)| *t >= start_time)
            .count();

        (births_in_window as f64 - deaths_in_window as f64) / time_window
    }
}

/// Population statistics
#[derive(Debug, Clone, Default)]
pub struct PopulationStats {
    /// Total population
    pub total: usize,

    /// Population by stage
    pub by_stage: Vec<(LifecycleStage, usize)>,

    /// Average age
    pub avg_age: SimTime,

    /// Average fitness
    pub avg_fitness: f64,

    /// Average generation
    pub avg_generation: f64,

    /// Diversity (unique genomes)
    pub diversity: f64,
}

/// Compute population statistics
pub fn compute_population_stats(agents: &[Sentinel]) -> PopulationStats {
    if agents.is_empty() {
        return PopulationStats::default();
    }

    let total = agents.len();

    // Count by stage
    let mut stage_counts = std::collections::HashMap::new();
    for agent in agents {
        *stage_counts.entry(agent.lifecycle()).or_insert(0) += 1;
    }
    let by_stage: Vec<_> = stage_counts.into_iter().collect();

    // Compute averages
    let sum_age: SimTime = agents.iter().map(|a| a.age()).sum();
    let sum_fitness: f64 = agents.iter().map(|a| a.fitness() as f64).sum();
    let sum_gen: u32 = agents.iter().map(|a| a.generation()).sum();

    PopulationStats {
        total,
        by_stage,
        avg_age: sum_age / total as f64,
        avg_fitness: sum_fitness / total as f64,
        avg_generation: sum_gen as f64 / total as f64,
        diversity: 1.0, // Would compute actual diversity
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_lifecycle_stages() {
        assert!(LifecycleStage::Adult.can_reproduce());
        assert!(!LifecycleStage::Juvenile.can_reproduce());
        assert!(LifecycleStage::Embryo.is_alive());
        assert!(!LifecycleStage::Dead.is_alive());
    }

    #[test]
    fn test_lifecycle_manager() {
        let mut manager = LifecycleManager::new();

        let config = SpawnConfig::default();
        let agent = manager.spawn(config, 0.0).unwrap();

        assert_eq!(agent.id(), 1);
        assert_eq!(manager.birth_count(), 1);
    }
}
