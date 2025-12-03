//! Performance Metrics
//!
//! Tracking and analysis of agent performance.

#[cfg(feature = "serde")]
use serde::{Deserialize, Serialize};

use crate::{AgentId, SimTime};

/// Agent performance metrics
#[derive(Debug, Clone, Default)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct AgentStats {
    /// Agent ID
    pub id: AgentId,

    /// Total reward accumulated
    pub total_reward: f64,

    /// Average reward per step
    pub avg_reward: f64,

    /// Maximum reward in single step
    pub max_reward: f64,

    /// Minimum reward in single step
    pub min_reward: f64,

    /// Total steps taken
    pub total_steps: u64,

    /// Survival time
    pub survival_time: SimTime,

    /// Number of offspring
    pub offspring_count: u32,

    /// Consciousness metrics
    pub avg_phi: f64,
    pub max_phi: f64,

    /// Energy efficiency (reward per energy)
    pub efficiency: f64,

    /// Behavior diversity
    pub diversity: f64,
}

impl AgentStats {
    pub fn new(id: AgentId) -> Self {
        Self {
            id,
            min_reward: f64::INFINITY,
            ..Default::default()
        }
    }

    /// Update with new reward
    pub fn update_reward(&mut self, reward: f64) {
        self.total_reward += reward;
        self.max_reward = self.max_reward.max(reward);
        self.min_reward = self.min_reward.min(reward);

        if self.total_steps > 0 {
            self.avg_reward = self.total_reward / self.total_steps as f64;
        }
    }

    /// Increment step count
    pub fn increment_steps(&mut self) {
        self.total_steps += 1;
    }

    /// Update survival time
    pub fn update_survival(&mut self, time: SimTime) {
        self.survival_time = time;
    }

    /// Update consciousness metrics
    pub fn update_consciousness(&mut self, phi: f64) {
        self.max_phi = self.max_phi.max(phi);

        // Running average
        if self.total_steps > 0 {
            let n = self.total_steps as f64;
            self.avg_phi = (self.avg_phi * (n - 1.0) + phi) / n;
        } else {
            self.avg_phi = phi;
        }
    }
}

/// Performance metrics tracker
#[derive(Debug, Clone, Default)]
pub struct PerformanceMetrics {
    /// Per-agent statistics
    agent_stats: Vec<AgentStats>,

    /// Global metrics
    pub global_reward: f64,
    pub global_steps: u64,
    pub total_agents: u64,
    pub active_agents: usize,

    /// Time series data
    reward_history: Vec<(SimTime, f64)>,
    population_history: Vec<(SimTime, usize)>,
    phi_history: Vec<(SimTime, f64)>,

    /// Moving averages
    reward_ma: f64,
    phi_ma: f64,

    /// Best ever metrics
    pub best_reward: f64,
    pub best_phi: f64,
    pub best_agent_id: Option<AgentId>,
}

impl PerformanceMetrics {
    pub fn new() -> Self {
        Self::default()
    }

    /// Record reward for an agent
    pub fn record_reward(&mut self, agent_id: AgentId, reward: f64, time: SimTime) {
        // Update agent stats
        if let Some(stats) = self.agent_stats.iter_mut().find(|s| s.id == agent_id) {
            stats.update_reward(reward);
        } else {
            let mut stats = AgentStats::new(agent_id);
            stats.update_reward(reward);
            self.agent_stats.push(stats);
        }

        // Update global metrics
        self.global_reward += reward;

        // Update history
        self.reward_history.push((time, reward));

        // Update moving average
        let alpha = 0.01;
        self.reward_ma = alpha * reward + (1.0 - alpha) * self.reward_ma;

        // Update best
        if reward > self.best_reward {
            self.best_reward = reward;
            self.best_agent_id = Some(agent_id);
        }
    }

    /// Record step for an agent
    pub fn record_step(&mut self, agent_id: AgentId) {
        if let Some(stats) = self.agent_stats.iter_mut().find(|s| s.id == agent_id) {
            stats.increment_steps();
        }
        self.global_steps += 1;
    }

    /// Record consciousness metric
    pub fn record_consciousness(&mut self, agent_id: AgentId, phi: f64, time: SimTime) {
        if let Some(stats) = self.agent_stats.iter_mut().find(|s| s.id == agent_id) {
            stats.update_consciousness(phi);
        }

        self.phi_history.push((time, phi));

        let alpha = 0.01;
        self.phi_ma = alpha * phi + (1.0 - alpha) * self.phi_ma;

        if phi > self.best_phi {
            self.best_phi = phi;
        }
    }

    /// Record population count
    pub fn record_population(&mut self, count: usize, time: SimTime) {
        self.active_agents = count;
        self.population_history.push((time, count));
    }

    /// Get agent statistics
    pub fn get_agent_stats(&self, agent_id: AgentId) -> Option<&AgentStats> {
        self.agent_stats.iter().find(|s| s.id == agent_id)
    }

    /// Get reward moving average
    pub fn reward_moving_average(&self) -> f64 {
        self.reward_ma
    }

    /// Get phi moving average
    pub fn phi_moving_average(&self) -> f64 {
        self.phi_ma
    }

    /// Get reward history
    pub fn reward_history(&self) -> &[(SimTime, f64)] {
        &self.reward_history
    }

    /// Get population history
    pub fn population_history(&self) -> &[(SimTime, usize)] {
        &self.population_history
    }

    /// Get phi history
    pub fn phi_history(&self) -> &[(SimTime, f64)] {
        &self.phi_history
    }

    /// Compute summary statistics
    pub fn summary(&self) -> MetricsSummary {
        let n = self.agent_stats.len();

        if n == 0 {
            return MetricsSummary::default();
        }

        let total_reward: f64 = self.agent_stats.iter().map(|s| s.total_reward).sum();
        let total_steps: u64 = self.agent_stats.iter().map(|s| s.total_steps).sum();
        let total_offspring: u32 = self.agent_stats.iter().map(|s| s.offspring_count).sum();
        let avg_phi: f64 = self.agent_stats.iter().map(|s| s.avg_phi).sum::<f64>() / n as f64;

        MetricsSummary {
            num_agents: n,
            total_reward,
            avg_reward_per_agent: total_reward / n as f64,
            total_steps,
            avg_steps_per_agent: total_steps as f64 / n as f64,
            total_offspring,
            avg_phi,
            best_phi: self.best_phi,
            reward_ma: self.reward_ma,
            phi_ma: self.phi_ma,
        }
    }

    /// Clear history (keep statistics)
    pub fn clear_history(&mut self) {
        self.reward_history.clear();
        self.population_history.clear();
        self.phi_history.clear();
    }

    /// Reset all metrics
    pub fn reset(&mut self) {
        *self = Self::new();
    }
}

/// Summary of all metrics
#[derive(Debug, Clone, Default)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct MetricsSummary {
    pub num_agents: usize,
    pub total_reward: f64,
    pub avg_reward_per_agent: f64,
    pub total_steps: u64,
    pub avg_steps_per_agent: f64,
    pub total_offspring: u32,
    pub avg_phi: f64,
    pub best_phi: f64,
    pub reward_ma: f64,
    pub phi_ma: f64,
}

/// Trading-specific metrics
#[derive(Debug, Clone, Default)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct TradingMetrics {
    /// Total P&L
    pub total_pnl: f64,

    /// Number of trades
    pub num_trades: u64,

    /// Win rate
    pub win_rate: f64,

    /// Average trade size
    pub avg_trade_size: f64,

    /// Maximum drawdown
    pub max_drawdown: f64,

    /// Sharpe ratio
    pub sharpe_ratio: f64,

    /// Decision latency (μs)
    pub avg_latency_us: f64,
    pub max_latency_us: f64,

    /// Signals generated
    pub signals_generated: u64,

    /// Market states observed
    pub states_observed: u64,
}

impl TradingMetrics {
    pub fn new() -> Self {
        Self::default()
    }

    /// Record trade result
    pub fn record_trade(&mut self, pnl: f64, size: f64) {
        self.total_pnl += pnl;
        self.num_trades += 1;

        if pnl > 0.0 {
            self.win_rate = (self.win_rate * (self.num_trades - 1) as f64 + 1.0)
                / self.num_trades as f64;
        } else {
            self.win_rate = (self.win_rate * (self.num_trades - 1) as f64)
                / self.num_trades as f64;
        }

        self.avg_trade_size = (self.avg_trade_size * (self.num_trades - 1) as f64 + size)
            / self.num_trades as f64;
    }

    /// Record decision latency
    pub fn record_latency(&mut self, latency_us: f64) {
        self.max_latency_us = self.max_latency_us.max(latency_us);

        let n = (self.signals_generated + 1) as f64;
        self.avg_latency_us = (self.avg_latency_us * (n - 1.0) + latency_us) / n;

        self.signals_generated += 1;
    }

    /// Update drawdown
    pub fn update_drawdown(&mut self, current_equity: f64, peak_equity: f64) {
        if peak_equity > 0.0 {
            let drawdown = (peak_equity - current_equity) / peak_equity;
            self.max_drawdown = self.max_drawdown.max(drawdown);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_agent_stats() {
        let mut stats = AgentStats::new(1);

        stats.update_reward(1.0);
        stats.update_reward(2.0);
        stats.increment_steps();
        stats.increment_steps();

        assert_eq!(stats.total_reward, 3.0);
        assert_eq!(stats.max_reward, 2.0);
        assert_eq!(stats.total_steps, 2);
    }

    #[test]
    fn test_performance_metrics() {
        let mut metrics = PerformanceMetrics::new();

        metrics.record_reward(1, 1.0, 0.0);
        metrics.record_reward(1, 2.0, 1.0);
        metrics.record_step(1);
        metrics.record_step(1);

        assert!(metrics.global_reward > 0.0);
        assert_eq!(metrics.global_steps, 2);
    }

    #[test]
    fn test_trading_metrics() {
        let mut metrics = TradingMetrics::new();

        metrics.record_trade(100.0, 1.0);
        metrics.record_trade(-50.0, 1.0);

        assert_eq!(metrics.num_trades, 2);
        assert_eq!(metrics.total_pnl, 50.0);
        assert!((metrics.win_rate - 0.5).abs() < 0.01);
    }
}
