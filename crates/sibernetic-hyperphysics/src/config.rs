//! Simulation Configuration
//!
//! Configuration types for the Sibernetic-HyperPhysics simulation.

#[cfg(feature = "serde")]
use serde::{Deserialize, Serialize};

use hyperphysics_connectome::ModelLevel;
use hyperphysics_embodiment::CouplingConfig;
use hyperphysics_sph::SphConfig;

/// Physics backend selection
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub enum PhysicsBackend {
    /// Native SPH solver (Sibernetic port)
    NativeSph,
    /// Rapier physics engine
    Rapier,
    /// Jolt physics engine
    Jolt,
    /// MuJoCo physics engine
    MuJoCo,
    /// PhysX physics engine
    PhysX,
    /// Custom backend
    Custom(u32),
}

impl Default for PhysicsBackend {
    fn default() -> Self {
        Self::NativeSph
    }
}

/// Neural simulation backend
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub enum NeuralBackend {
    /// Native spiking network (c302 port)
    NativeSpiking,
    /// NEURON simulator integration
    Neuron,
    /// Brian2 integration
    Brian2,
    /// PyNN integration
    PyNN,
    /// Custom backend
    Custom(u32),
}

impl Default for NeuralBackend {
    fn default() -> Self {
        Self::NativeSpiking
    }
}

/// Main simulation configuration
#[derive(Debug, Clone)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct SimulationConfig {
    /// Physics backend
    pub physics_backend: PhysicsBackend,

    /// Neural backend
    pub neural_backend: NeuralBackend,

    /// Model complexity level
    pub model_level: ModelLevel,

    /// SPH physics configuration
    pub sph_config: SphConfig,

    /// Neural-body coupling configuration
    pub coupling_config: CouplingConfig,

    /// Simulation time step (ms)
    pub dt: f32,

    /// Enable real-time mode
    pub real_time: bool,

    /// Enable visualization
    pub visualization: bool,

    /// Enable learning/plasticity
    pub learning_enabled: bool,

    /// Enable evolution
    pub evolution_enabled: bool,

    /// Random seed
    pub seed: Option<u64>,

    /// Maximum simulation time (None = unlimited)
    pub max_time: Option<f64>,

    /// Output configuration
    pub output: OutputConfig,
}

impl Default for SimulationConfig {
    fn default() -> Self {
        Self {
            physics_backend: PhysicsBackend::NativeSph,
            neural_backend: NeuralBackend::NativeSpiking,
            model_level: ModelLevel::B,
            sph_config: SphConfig::celegans(),
            coupling_config: CouplingConfig::default(),
            dt: 0.5,
            real_time: false,
            visualization: false,
            learning_enabled: false,
            evolution_enabled: false,
            seed: None,
            max_time: None,
            output: OutputConfig::default(),
        }
    }
}

impl SimulationConfig {
    /// Fast simulation configuration
    pub fn fast() -> Self {
        Self {
            model_level: ModelLevel::A,
            dt: 1.0,
            coupling_config: CouplingConfig::open_loop(),
            ..Default::default()
        }
    }

    /// High-fidelity configuration
    pub fn high_fidelity() -> Self {
        Self {
            model_level: ModelLevel::C,
            dt: 0.1,
            coupling_config: CouplingConfig::high_fidelity(),
            ..Default::default()
        }
    }

    /// Research configuration (full features)
    pub fn research() -> Self {
        Self {
            model_level: ModelLevel::C,
            dt: 0.25,
            coupling_config: CouplingConfig::high_fidelity(),
            learning_enabled: true,
            evolution_enabled: false,
            visualization: true,
            output: OutputConfig::full(),
            ..Default::default()
        }
    }

    /// Demo configuration
    pub fn demo() -> Self {
        Self {
            model_level: ModelLevel::B,
            dt: 0.5,
            real_time: true,
            visualization: true,
            ..Default::default()
        }
    }

    /// Validate configuration
    pub fn validate(&self) -> Result<(), String> {
        if self.dt <= 0.0 {
            return Err("Time step must be positive".to_string());
        }

        if self.dt > 10.0 {
            return Err("Time step too large (max 10 ms)".to_string());
        }

        Ok(())
    }
}

/// Output configuration
#[derive(Debug, Clone)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct OutputConfig {
    /// Record membrane potentials
    pub record_voltages: bool,

    /// Record spike times
    pub record_spikes: bool,

    /// Record muscle activations
    pub record_muscles: bool,

    /// Record body positions
    pub record_positions: bool,

    /// Recording interval (steps)
    pub recording_interval: u32,

    /// Output file path
    pub output_path: Option<String>,
}

impl Default for OutputConfig {
    fn default() -> Self {
        Self {
            record_voltages: false,
            record_spikes: true,
            record_muscles: true,
            record_positions: true,
            recording_interval: 10,
            output_path: None,
        }
    }
}

impl OutputConfig {
    /// Minimal output
    pub fn minimal() -> Self {
        Self {
            record_voltages: false,
            record_spikes: false,
            record_muscles: true,
            record_positions: false,
            recording_interval: 100,
            output_path: None,
        }
    }

    /// Full output
    pub fn full() -> Self {
        Self {
            record_voltages: true,
            record_spikes: true,
            record_muscles: true,
            record_positions: true,
            recording_interval: 1,
            output_path: None,
        }
    }
}

/// Worm body configuration
#[derive(Debug, Clone)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct WormBodyConfig {
    /// Number of body segments
    pub num_segments: usize,

    /// Particles per segment
    pub particles_per_segment: usize,

    /// Body length (mm)
    pub length: f32,

    /// Body radius (mm)
    pub radius: f32,

    /// Muscle stiffness
    pub muscle_stiffness: f32,

    /// Cuticle stiffness
    pub cuticle_stiffness: f32,
}

impl Default for WormBodyConfig {
    fn default() -> Self {
        Self {
            num_segments: 24,
            particles_per_segment: 8,
            length: 1.0,  // 1 mm
            radius: 0.04, // 40 μm
            muscle_stiffness: 200.0,
            cuticle_stiffness: 100.0,
        }
    }
}

/// Environment configuration
#[derive(Debug, Clone)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct EnvironmentConfig {
    /// Environment type
    pub env_type: EnvironmentType,

    /// Fluid viscosity
    pub viscosity: f32,

    /// Gravity vector
    pub gravity: [f32; 3],

    /// Boundary conditions
    pub boundaries: BoundaryConditions,

    /// Add food sources
    pub food_sources: Vec<FoodSource>,

    /// Add obstacles
    pub obstacles: Vec<Obstacle>,
}

impl Default for EnvironmentConfig {
    fn default() -> Self {
        Self {
            env_type: EnvironmentType::Water,
            viscosity: 0.001,
            gravity: [0.0, -9.81, 0.0],
            boundaries: BoundaryConditions::default(),
            food_sources: Vec::new(),
            obstacles: Vec::new(),
        }
    }
}

/// Environment type
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub enum EnvironmentType {
    /// Swimming in water
    Water,
    /// Crawling on agar
    Agar,
    /// Crawling in soil
    Soil,
    /// Custom environment
    Custom,
}

/// Boundary conditions
#[derive(Debug, Clone)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct BoundaryConditions {
    /// Minimum coordinates
    pub min: [f32; 3],
    /// Maximum coordinates
    pub max: [f32; 3],
    /// Boundary type
    pub boundary_type: BoundaryType,
}

impl Default for BoundaryConditions {
    fn default() -> Self {
        Self {
            min: [-10.0, -10.0, -10.0],
            max: [10.0, 10.0, 10.0],
            boundary_type: BoundaryType::Reflective,
        }
    }
}

/// Boundary type
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub enum BoundaryType {
    /// Particles bounce off
    Reflective,
    /// Particles wrap around
    Periodic,
    /// Particles are absorbed
    Absorbing,
    /// No boundary
    Open,
}

/// Food source in environment
#[derive(Debug, Clone)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct FoodSource {
    pub position: [f32; 3],
    pub concentration: f32,
    pub radius: f32,
}

/// Obstacle in environment
#[derive(Debug, Clone)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct Obstacle {
    pub position: [f32; 3],
    pub size: [f32; 3],
    pub shape: ObstacleShape,
}

/// Obstacle shape
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub enum ObstacleShape {
    Box,
    Sphere,
    Cylinder,
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_config_validation() {
        let config = SimulationConfig::default();
        assert!(config.validate().is_ok());

        let mut bad_config = SimulationConfig::default();
        bad_config.dt = -1.0;
        assert!(bad_config.validate().is_err());
    }

    #[test]
    fn test_presets() {
        let fast = SimulationConfig::fast();
        assert!(fast.dt >= 1.0);

        let hifi = SimulationConfig::high_fidelity();
        assert!(hifi.dt <= 0.25);
    }
}
