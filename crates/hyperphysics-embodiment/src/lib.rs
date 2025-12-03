//! HyperPhysics Embodiment - Neural-Body Coupling
//!
//! This crate bridges the neural network simulation (connectome) with the
//! body physics simulation (SPH), enabling closed-loop embodied simulation.
//!
//! # Architecture
//!
//! ```text
//! ┌─────────────────────────────────────────────────────────────────┐
//! │                    Embodied Simulation                          │
//! │                                                                 │
//! │  ┌──────────────┐      Muscle       ┌──────────────────┐       │
//! │  │   Neural     │   Activation      │      Body        │       │
//! │  │  Network     │ ───────────────►  │    Physics       │       │
//! │  │  (302 neurons)│                   │   (SPH particles)│       │
//! │  │              │ ◄───────────────  │                  │       │
//! │  └──────────────┘  Proprioception   └──────────────────┘       │
//! │                                                                 │
//! └─────────────────────────────────────────────────────────────────┘
//! ```
//!
//! # Multi-Rate Integration
//!
//! Neural dynamics typically require finer time steps than body physics:
//! - Neural: dt = 0.025 - 0.1 ms (depending on model level)
//! - Physics: dt = 0.5 - 1.0 ms
//!
//! The embodiment layer handles this by running multiple neural steps
//! per physics step with proper interpolation.

mod coupling;
mod proprioception;
mod actuator;
mod embodied;
mod time_sync;

pub use coupling::{CouplingConfig, CouplingMode};
pub use proprioception::{Proprioceptor, ProprioceptiveState, StretchReceptor};
pub use actuator::{Actuator, ActuatorConfig};
pub use embodied::{EmbodiedWorm, SimulationState};
pub use time_sync::{TimeSync, IntegrationStrategy};

/// Re-export key types from sub-crates
pub use hyperphysics_sph::{SphWorld, MuscleActivation};
pub use hyperphysics_connectome::{SpikingNetwork, ModelLevel, Connectome};

use thiserror::Error;

/// Embodiment errors
#[derive(Debug, Error)]
pub enum EmbodimentError {
    #[error("SPH error: {0}")]
    Sph(#[from] hyperphysics_sph::SphError),

    #[error("Muscle mapping mismatch: expected {expected} muscles, got {actual}")]
    MuscleCountMismatch { expected: usize, actual: usize },

    #[error("Invalid segment index: {index} (max: {max})")]
    InvalidSegment { index: usize, max: usize },

    #[error("Time synchronization error: neural dt ({neural_dt}) > physics dt ({physics_dt})")]
    TimeSyncError { neural_dt: f64, physics_dt: f64 },
}

pub type Result<T> = std::result::Result<T, EmbodimentError>;
