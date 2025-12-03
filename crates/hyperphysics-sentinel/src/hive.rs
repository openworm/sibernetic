//! Hive - Multi-Agent Coordination
//!
//! Manages populations of Sentinel agents with swarm behavior.

#[cfg(feature = "serde")]
use serde::{Deserialize, Serialize};

use hashbrown::HashMap;

use crate::agent::{AgentState, Sentinel, SentinelConfig};
use crate::lifecycle::{LifecycleManager, PopulationStats, SpawnConfig, compute_population_stats};
use crate::{AgentId, Result, SimTime};

/// Swarm behavior type
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub enum SwarmBehavior {
    /// Independent agents
    Independent,
    /// Cooperative behavior
    Cooperative,
    /// Competitive behavior
    Competitive,
    /// Hierarchical (leader-follower)
    Hierarchical,
    /// Stigmergic (environment-mediated)
    Stigmergic,
}

/// Hive configuration
#[derive(Debug, Clone)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct HiveConfig {
    /// Maximum population
    pub max_population: usize,

    /// Minimum population
    pub min_population: usize,

    /// Default agent configuration
    pub agent_config: SentinelConfig,

    /// Swarm behavior
    pub swarm_behavior: SwarmBehavior,

    /// Enable evolution
    pub evolution_enabled: bool,

    /// Selection pressure (0-1)
    pub selection_pressure: f32,

    /// Reproduction rate
    pub reproduction_rate: f32,

    /// Death rate
    pub death_rate: f32,
}

impl Default for HiveConfig {
    fn default() -> Self {
        Self {
            max_population: 100,
            min_population: 10,
            agent_config: SentinelConfig::default(),
            swarm_behavior: SwarmBehavior::Cooperative,
            evolution_enabled: true,
            selection_pressure: 0.5,
            reproduction_rate: 0.1,
            death_rate: 0.05,
        }
    }
}

/// A hive of Sentinel agents
pub struct Hive {
    /// Configuration
    config: HiveConfig,

    /// Agents
    agents: HashMap<AgentId, Sentinel>,

    /// Lifecycle manager
    lifecycle: LifecycleManager,

    /// Current simulation time
    time: SimTime,

    /// Generation counter
    generation: u32,

    /// Global reward buffer
    global_reward: f64,

    /// Leader agent (for hierarchical behavior)
    leader_id: Option<AgentId>,
}

impl Hive {
    /// Create new hive
    pub fn new(config: HiveConfig) -> Self {
        Self {
            config,
            agents: HashMap::new(),
            lifecycle: LifecycleManager::new(),
            time: 0.0,
            generation: 0,
            global_reward: 0.0,
            leader_id: None,
        }
    }

    /// Initialize hive with random agents
    pub fn initialize(&mut self, count: usize) -> Result<()> {
        for _ in 0..count.min(self.config.max_population) {
            self.spawn_agent(SpawnConfig {
                agent_config: self.config.agent_config.clone(),
                generation: 0,
                ..Default::default()
            })?;
        }
        Ok(())
    }

    /// Spawn new agent
    pub fn spawn_agent(&mut self, config: SpawnConfig) -> Result<AgentId> {
        if self.agents.len() >= self.config.max_population {
            return Ok(0); // No room
        }

        let mut agent = self.lifecycle.spawn(config, self.time)?;
        agent.initialize()?;
        let id = agent.id();
        self.agents.insert(id, agent);
        Ok(id)
    }

    /// Remove agent
    pub fn remove_agent(&mut self, id: AgentId) {
        if let Some(agent) = self.agents.remove(&id) {
            self.lifecycle.register_death(agent.id(), self.time);

            if self.leader_id == Some(id) {
                self.elect_leader();
            }
        }
    }

    /// Step all agents
    pub fn step(&mut self, dt: f32) {
        // Step each agent
        let ids: Vec<AgentId> = self.agents.keys().copied().collect();

        for id in &ids {
            if let Some(agent) = self.agents.get_mut(id) {
                agent.step(dt);
                self.lifecycle.update_stage(agent);
            }
        }

        // Remove dead agents
        let dead_ids: Vec<AgentId> = self
            .agents
            .iter()
            .filter(|(_, a)| !a.is_alive())
            .map(|(id, _)| *id)
            .collect();

        for id in dead_ids {
            self.remove_agent(id);
        }

        // Maintain minimum population
        while self.agents.len() < self.config.min_population {
            let _ = self.spawn_agent(SpawnConfig {
                agent_config: self.config.agent_config.clone(),
                generation: self.generation,
                ..Default::default()
            });
        }

        // Swarm behavior
        self.apply_swarm_behavior();

        // Evolution
        if self.config.evolution_enabled {
            self.evolve();
        }

        self.time += dt as f64;
    }

    /// Run simulation for duration
    pub fn run(&mut self, duration: SimTime, dt: f32) {
        let steps = (duration / dt as f64).ceil() as usize;
        for _ in 0..steps {
            self.step(dt);
        }
    }

    /// Apply swarm behavior
    fn apply_swarm_behavior(&mut self) {
        match self.config.swarm_behavior {
            SwarmBehavior::Independent => {}
            SwarmBehavior::Cooperative => self.apply_cooperative_behavior(),
            SwarmBehavior::Competitive => self.apply_competitive_behavior(),
            SwarmBehavior::Hierarchical => self.apply_hierarchical_behavior(),
            SwarmBehavior::Stigmergic => self.apply_stigmergic_behavior(),
        }
    }

    /// Cooperative behavior: share rewards
    fn apply_cooperative_behavior(&mut self) {
        if self.agents.is_empty() {
            return;
        }

        // Average reward distribution
        let total_reward: f64 = self.agents.values().map(|a| a.total_reward()).sum();
        let avg_reward = total_reward / self.agents.len() as f64;

        self.global_reward = avg_reward;
    }

    /// Competitive behavior: winner takes all
    fn apply_competitive_behavior(&mut self) {
        // Find best performing agent
        if let Some((&best_id, _)) = self.agents.iter().max_by(|a, b| {
            a.1.fitness().partial_cmp(&b.1.fitness()).unwrap_or(std::cmp::Ordering::Equal)
        }) {
            self.leader_id = Some(best_id);
        }
    }

    /// Hierarchical behavior
    fn apply_hierarchical_behavior(&mut self) {
        if self.leader_id.is_none() {
            self.elect_leader();
        }

        // Leader influences followers
        if let Some(leader_id) = self.leader_id {
            if let Some(leader) = self.agents.get(&leader_id) {
                let leader_action = leader.embodiment().muscle_activations().clone();
                // Could modify follower behavior based on leader
            }
        }
    }

    /// Stigmergic behavior (environment-based communication)
    fn apply_stigmergic_behavior(&mut self) {
        // Agents modify and respond to environment
        // Would integrate with environment state
    }

    /// Elect leader (for hierarchical behavior)
    fn elect_leader(&mut self) {
        if let Some((&best_id, _)) = self.agents.iter().max_by(|a, b| {
            a.1.fitness().partial_cmp(&b.1.fitness()).unwrap_or(std::cmp::Ordering::Equal)
        }) {
            self.leader_id = Some(best_id);
        } else {
            self.leader_id = None;
        }
    }

    /// Perform evolutionary step
    fn evolve(&mut self) {
        use rand::Rng;
        let mut rng = rand::thread_rng();

        // Selection based on fitness
        let mut agents_by_fitness: Vec<(&AgentId, &Sentinel)> = self.agents.iter().collect();
        agents_by_fitness.sort_by(|a, b| {
            b.1.fitness().partial_cmp(&a.1.fitness()).unwrap_or(std::cmp::Ordering::Equal)
        });

        // Reproduction
        if self.agents.len() < self.config.max_population {
            let reproduction_prob = self.config.reproduction_rate;

            for (id, agent) in &agents_by_fitness[..agents_by_fitness.len().min(10)] {
                if self.lifecycle.can_reproduce(agent, self.time)
                    && rng.gen::<f32>() < reproduction_prob
                {
                    let config = SpawnConfig {
                        agent_config: self.config.agent_config.clone(),
                        parents: vec![**id],
                        generation: agent.generation() + 1,
                        ..Default::default()
                    };

                    let _ = self.spawn_agent(config);
                    break;
                }
            }
        }

        // Selection (death of unfit)
        let death_prob = self.config.death_rate * self.config.selection_pressure;

        let weak_agents: Vec<AgentId> = agents_by_fitness
            .iter()
            .rev()
            .take(agents_by_fitness.len() / 4)
            .filter(|(_, a)| a.fitness() < 0.0 && rng.gen::<f32>() < death_prob)
            .map(|(id, _)| **id)
            .collect();

        for id in weak_agents {
            self.remove_agent(id);
        }

        self.generation += 1;
    }

    /// Distribute reward to all agents
    pub fn distribute_reward(&mut self, reward: f64) {
        use crate::experience::Reward;

        self.global_reward += reward;

        let per_agent = reward / self.agents.len().max(1) as f64;

        for agent in self.agents.values_mut() {
            agent.receive_reward(Reward::new(per_agent, self.time));
        }
    }

    // ========== Getters ==========

    /// Get agent by ID
    pub fn get_agent(&self, id: AgentId) -> Option<&Sentinel> {
        self.agents.get(&id)
    }

    /// Get mutable agent by ID
    pub fn get_agent_mut(&mut self, id: AgentId) -> Option<&mut Sentinel> {
        self.agents.get_mut(&id)
    }

    /// Get all agents
    pub fn agents(&self) -> impl Iterator<Item = &Sentinel> {
        self.agents.values()
    }

    /// Get population size
    pub fn population(&self) -> usize {
        self.agents.len()
    }

    /// Get current time
    pub fn time(&self) -> SimTime {
        self.time
    }

    /// Get generation
    pub fn generation(&self) -> u32 {
        self.generation
    }

    /// Get leader ID
    pub fn leader(&self) -> Option<AgentId> {
        self.leader_id
    }

    /// Get population statistics
    pub fn stats(&self) -> PopulationStats {
        let agents: Vec<_> = self.agents.values().collect();
        // Note: compute_population_stats expects &[Sentinel], not &[&Sentinel]
        // This is simplified for the example
        PopulationStats {
            total: self.agents.len(),
            avg_fitness: self.agents.values().map(|a| a.fitness() as f64).sum::<f64>()
                / self.agents.len().max(1) as f64,
            avg_generation: self.agents.values().map(|a| a.generation() as f64).sum::<f64>()
                / self.agents.len().max(1) as f64,
            ..Default::default()
        }
    }

    /// Get best agent
    pub fn best_agent(&self) -> Option<&Sentinel> {
        self.agents.values().max_by(|a, b| {
            a.fitness().partial_cmp(&b.fitness()).unwrap_or(std::cmp::Ordering::Equal)
        })
    }

    /// Get global reward
    pub fn global_reward(&self) -> f64 {
        self.global_reward
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_hive_creation() {
        let config = HiveConfig {
            max_population: 10,
            min_population: 5,
            ..Default::default()
        };

        let mut hive = Hive::new(config);
        hive.initialize(5).unwrap();

        assert_eq!(hive.population(), 5);
    }

    #[test]
    fn test_hive_step() {
        let config = HiveConfig {
            max_population: 10,
            min_population: 2,
            agent_config: SentinelConfig::fast(),
            ..Default::default()
        };

        let mut hive = Hive::new(config);
        hive.initialize(3).unwrap();

        hive.step(0.5);

        assert!(hive.time() > 0.0);
    }
}
