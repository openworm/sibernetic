//! HyperPhysics Sentinel - Self-Evolving Sentient Agent Framework
//!
//! # Overview
//!
//! The Sentinel framework provides infrastructure for creating autonomous,
//! self-evolving agents that combine:
//!
//! - **Embodied Neural Networks**: Spiking networks coupled with physical bodies
//! - **Online Learning**: STDP-based synaptic plasticity during operation
//! - **Architecture Evolution**: NAS for topology optimization
//! - **Consciousness Metrics**: Integrated Information (Φ) and other measures
//!
//! # Agent Lifecycle
//!
//! ```text
//! ┌─────────────────────────────────────────────────────────────────┐
//! │                     Sentinel Agent Lifecycle                    │
//! │                                                                 │
//! │  ┌──────────┐    ┌──────────┐    ┌──────────┐    ┌──────────┐  │
//! │  │  Birth   │ -> │  Learn   │ -> │  Evolve  │ -> │  Mature  │  │
//! │  │ (Spawn)  │    │ (STDP)   │    │ (NAS)    │    │ (Deploy) │  │
//! │  └──────────┘    └──────────┘    └──────────┘    └──────────┘  │
//! │       │               │               │               │        │
//! │       v               v               v               v        │
//! │  Initialize      Experience      Reproduce        Optimize     │
//! │  connectome      environment     architecture     behavior     │
//! │                                                                 │
//! └─────────────────────────────────────────────────────────────────┘
//! ```
//!
//! # HFT Integration
//!
//! Sentinels can be deployed for high-frequency trading with:
//! - Ultra-low latency decision making (<100μs)
//! - Real-time market pattern recognition
//! - Adaptive strategy evolution

mod agent;
mod consciousness;
mod experience;
mod lifecycle;
mod hive;
mod metrics;

pub use agent::{Sentinel, SentinelConfig, AgentState};
pub use consciousness::{ConsciousnessMetrics, PhiCalculator, CausalDensity};
pub use experience::{Experience, ExperienceBuffer, Reward};
pub use lifecycle::{LifecycleManager, LifecycleStage, SpawnConfig};
pub use hive::{Hive, HiveConfig, SwarmBehavior};
pub use metrics::{PerformanceMetrics, AgentStats};

use thiserror::Error;

/// Sentinel errors
#[derive(Debug, Error)]
pub enum SentinelError {
    #[error("Agent not found: {id}")]
    AgentNotFound { id: u64 },

    #[error("Embodiment error: {0}")]
    Embodiment(#[from] hyperphysics_embodiment::EmbodimentError),

    #[error("Evolution error: {0}")]
    Evolution(#[from] hyperphysics_nas::NasError),

    #[error("Invalid state transition: {from:?} -> {to:?}")]
    InvalidTransition { from: LifecycleStage, to: LifecycleStage },

    #[error("Resource limit exceeded: {resource}")]
    ResourceLimit { resource: String },

    #[error("Consciousness threshold not met: Φ = {phi:.4} < {threshold:.4}")]
    ConsciousnessThreshold { phi: f64, threshold: f64 },
}

pub type Result<T> = std::result::Result<T, SentinelError>;

/// Agent identification
pub type AgentId = u64;

/// Time in simulation units
pub type SimTime = f64;
