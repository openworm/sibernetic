//! # HyperPhysics Connectome
//!
//! C. elegans connectome data and spiking neural network models,
//! ported from the OpenWorm c302 project.
//!
//! ## Features
//!
//! - Complete C. elegans connectome (302 neurons, ~7000 synapses)
//! - Multiple neuron models (LIF, Izhikevich, Hodgkin-Huxley)
//! - Synaptic dynamics (chemical and electrical)
//! - Neuromuscular junction mapping
//!
//! ## Model Levels (from c302)
//!
//! - **Level A**: Integrate-and-fire neurons (fast, for behavior studies)
//! - **Level B**: Izhikevich neurons (captures spiking dynamics)
//! - **Level C**: Conductance-based (includes basic ion channels)
//! - **Level D**: Full Hodgkin-Huxley (biophysically detailed)
//!
//! ## Example
//!
//! ```rust,no_run
//! use hyperphysics_connectome::{Connectome, SpikingNetwork, ModelLevel};
//!
//! // Load the C. elegans connectome
//! let connectome = Connectome::celegans();
//!
//! // Create a spiking network with Izhikevich neurons
//! let mut network = SpikingNetwork::from_connectome(&connectome, ModelLevel::B);
//!
//! // Run simulation
//! let dt = 0.1; // ms
//! for _ in 0..1000 {
//!     network.step(dt);
//! }
//!
//! // Get motor neuron outputs for muscle control
//! let muscle_signals = network.get_muscle_output();
//! ```

#![warn(missing_docs)]

pub mod neuron;
pub mod synapse;
pub mod connectome;
pub mod network;
pub mod models;
pub mod muscle_map;

// Re-exports
pub use neuron::{Neuron, NeuronId, NeuronClass, NeuronState};
pub use synapse::{Synapse, SynapseType, SynapticState};
pub use connectome::Connectome;
pub use network::SpikingNetwork;
pub use models::{ModelLevel, ModelParams, NeuronModel};
pub use muscle_map::{MuscleMap, NeuromuscularJunction};

/// Result type for connectome operations
pub type Result<T> = std::result::Result<T, ConnectomeError>;

/// Errors that can occur in connectome operations
#[derive(Debug, thiserror::Error)]
pub enum ConnectomeError {
    /// Neuron not found
    #[error("Neuron not found: {0}")]
    NeuronNotFound(String),

    /// Invalid model level
    #[error("Invalid model level: {0}")]
    InvalidModelLevel(String),

    /// Simulation error
    #[error("Simulation error: {0}")]
    Simulation(String),

    /// IO error
    #[error("IO error: {0}")]
    Io(#[from] std::io::Error),

    /// Serialization error
    #[cfg(feature = "serde")]
    #[error("Serialization error: {0}")]
    Serialization(String),
}
