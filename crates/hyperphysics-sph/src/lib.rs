//! # HyperPhysics SPH
//!
//! Smoothed Particle Hydrodynamics (SPH) engine for soft body and fluid simulation.
//! Ported from OpenWorm's Sibernetic project with enhancements for HyperPhysics integration.
//!
//! ## Features
//!
//! - **PCISPH**: Predictive-Corrective Incompressible SPH for stable fluid simulation
//! - **Elastic Bodies**: Spring-based connections for deformable solids
//! - **Membranes**: Surface tension and boundary handling
//! - **Muscle Simulation**: Contractile forces for biomechanical modeling
//! - **Multi-backend**: CPU (SIMD), wgpu, OpenCL support
//!
//! ## References
//!
//! - Solenthaler & Pajarola (2009): Predictive-Corrective Incompressible SPH
//! - Müller et al. (2003): Particle-based fluid simulation for interactive applications
//! - Ihmsen et al. (2010): Boundary handling and adaptive time-stepping for PCISPH
//!
//! ## Example
//!
//! ```rust,no_run
//! use hyperphysics_sph::{SphWorld, SphConfig, ParticleType};
//!
//! let config = SphConfig::default();
//! let mut world = SphWorld::new(config);
//!
//! // Add liquid particles
//! world.add_particle([0.0, 1.0, 0.0], [0.0, 0.0, 0.0], ParticleType::Liquid);
//!
//! // Simulation loop
//! for _ in 0..1000 {
//!     world.step();
//! }
//! ```

#![warn(missing_docs)]
#![allow(clippy::excessive_precision)]

pub mod config;
pub mod particle;
pub mod spatial_hash;
pub mod kernels;
pub mod pcisph;
pub mod elastic;
pub mod membrane;
pub mod muscle;
pub mod world;

#[cfg(feature = "gpu-wgpu")]
pub mod gpu_wgpu;

#[cfg(feature = "gpu-opencl")]
pub mod gpu_opencl;

// Re-exports
pub use config::{SphConfig, PhysicsConstants};
pub use particle::{Particle, ParticleType, ParticleBuffer};
pub use spatial_hash::SpatialHash;
pub use elastic::{ElasticConnection, ElasticNetwork};
pub use membrane::{Membrane, MembraneTriangle};
pub use muscle::{MuscleActivation, MuscleSegment};
pub use world::SphWorld;

/// Result type for SPH operations
pub type Result<T> = std::result::Result<T, SphError>;

/// Errors that can occur in SPH simulation
#[derive(Debug, thiserror::Error)]
pub enum SphError {
    /// Configuration error
    #[error("Invalid configuration: {0}")]
    Config(String),

    /// Particle limit exceeded
    #[error("Maximum particle count exceeded: {count} > {max}")]
    ParticleLimit { count: usize, max: usize },

    /// GPU error
    #[error("GPU error: {0}")]
    Gpu(String),

    /// IO error
    #[error("IO error: {0}")]
    Io(#[from] std::io::Error),

    /// Serialization error
    #[error("Serialization error: {0}")]
    Serialization(String),
}
