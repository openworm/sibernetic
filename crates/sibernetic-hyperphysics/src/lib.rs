//! Sibernetic-HyperPhysics Backend Adapter
//!
//! This crate integrates the Sibernetic C. elegans simulation with the
//! HyperPhysics ecosystem, providing a unified interface for embodied
//! neural simulation.
//!
//! # Architecture
//!
//! ```text
//! ┌─────────────────────────────────────────────────────────────────┐
//! │                    HyperPhysics Ecosystem                       │
//! │                                                                 │
//! │  ┌─────────────┐  ┌─────────────┐  ┌─────────────┐             │
//! │  │   Rapier    │  │    Jolt     │  │   MuJoCo    │             │
//! │  │  Backend    │  │  Backend    │  │  Backend    │  ...        │
//! │  └─────────────┘  └─────────────┘  └─────────────┘             │
//! │         │               │               │                       │
//! │         └───────────────┴───────────────┘                       │
//! │                         │                                       │
//! │              ┌──────────────────┐                               │
//! │              │ Physics Backend  │                               │
//! │              │    Interface     │                               │
//! │              └──────────────────┘                               │
//! │                         │                                       │
//! │  ┌──────────────────────┴──────────────────────┐               │
//! │  │         Sibernetic-HyperPhysics             │               │
//! │  │  ┌─────────────┐  ┌─────────────────────┐   │               │
//! │  │  │    SPH      │  │    Connectome       │   │               │
//! │  │  │  Simulator  │  │    (c302 port)      │   │               │
//! │  │  └─────────────┘  └─────────────────────┘   │               │
//! │  │         │                   │               │               │
//! │  │         └───────┬───────────┘               │               │
//! │  │                 │                           │               │
//! │  │         ┌───────────────┐                   │               │
//! │  │         │  Embodiment   │                   │               │
//! │  │         └───────────────┘                   │               │
//! │  └─────────────────────────────────────────────┘               │
//! │                                                                 │
//! └─────────────────────────────────────────────────────────────────┘
//! ```
//!
//! # Usage
//!
//! ```rust,ignore
//! use sibernetic_hyperphysics::{WormSimulation, SimulationConfig};
//!
//! let mut sim = WormSimulation::new(SimulationConfig::default());
//! sim.initialize()?;
//!
//! // Run simulation
//! for _ in 0..1000 {
//!     sim.step();
//! }
//!
//! // Get muscle output for visualization
//! let muscles = sim.get_muscle_activations();
//! ```

mod simulation;
mod config;
mod bridge;
mod visualization;
#[cfg(feature = "python")]
mod python;

pub use simulation::{WormSimulation, SimulationMode};
pub use config::{SimulationConfig, PhysicsBackend, NeuralBackend};
pub use bridge::{BackendBridge, BackendCapabilities};
pub use visualization::{VisualizationState, RenderData};

// Re-export key types from sub-crates
pub use hyperphysics_sph::{SphWorld, MuscleActivation, ParticleType};
pub use hyperphysics_connectome::{SpikingNetwork, ModelLevel, Connectome};
pub use hyperphysics_embodiment::EmbodiedWorm;
pub use hyperphysics_sentinel::{Sentinel, Hive, SentinelConfig, HiveConfig};

use thiserror::Error;

/// Simulation errors
#[derive(Debug, Error)]
pub enum SimulationError {
    #[error("SPH error: {0}")]
    Sph(#[from] hyperphysics_sph::SphError),

    #[error("Embodiment error: {0}")]
    Embodiment(#[from] hyperphysics_embodiment::EmbodimentError),

    #[error("Sentinel error: {0}")]
    Sentinel(#[from] hyperphysics_sentinel::SentinelError),

    #[error("Backend not available: {backend}")]
    BackendNotAvailable { backend: String },

    #[error("Configuration error: {0}")]
    Configuration(String),

    #[error("Initialization failed: {0}")]
    Initialization(String),
}

pub type Result<T> = std::result::Result<T, SimulationError>;

/// Simulation time in milliseconds
pub type SimTime = f64;

/// Library version
pub const VERSION: &str = env!("CARGO_PKG_VERSION");
