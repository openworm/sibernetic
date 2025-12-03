//! Backend Bridge
//!
//! Abstract interface for different physics and neural backends.

use crate::config::{PhysicsBackend, NeuralBackend};

/// Backend capabilities
#[derive(Debug, Clone, Default)]
pub struct BackendCapabilities {
    /// Supports GPU acceleration
    pub gpu: bool,

    /// Supports multi-threading
    pub multithreaded: bool,

    /// Supports SIMD
    pub simd: bool,

    /// Maximum particles (physics)
    pub max_particles: usize,

    /// Maximum neurons (neural)
    pub max_neurons: usize,

    /// Real-time capable
    pub real_time: bool,

    /// Supports deterministic execution
    pub deterministic: bool,
}

impl BackendCapabilities {
    /// Native SPH capabilities
    pub fn native_sph() -> Self {
        Self {
            gpu: false, // CPU only for now
            multithreaded: false, // Single-threaded
            simd: true, // Uses SIMD where available
            max_particles: 100_000,
            max_neurons: 0,
            real_time: true,
            deterministic: true,
        }
    }

    /// Native spiking network capabilities
    pub fn native_spiking() -> Self {
        Self {
            gpu: false,
            multithreaded: false,
            simd: true,
            max_particles: 0,
            max_neurons: 10_000,
            real_time: true,
            deterministic: true,
        }
    }
}

/// Backend bridge interface
pub struct BackendBridge {
    /// Physics backend
    physics: PhysicsBackend,

    /// Neural backend
    neural: NeuralBackend,

    /// Capabilities
    capabilities: BackendCapabilities,
}

impl BackendBridge {
    /// Create native backend bridge
    pub fn native() -> Self {
        Self {
            physics: PhysicsBackend::NativeSph,
            neural: NeuralBackend::NativeSpiking,
            capabilities: BackendCapabilities {
                gpu: false,
                multithreaded: false,
                simd: true,
                max_particles: 100_000,
                max_neurons: 10_000,
                real_time: true,
                deterministic: true,
            },
        }
    }

    /// Create with specific backends
    pub fn new(physics: PhysicsBackend, neural: NeuralBackend) -> Self {
        let capabilities = match (&physics, &neural) {
            (PhysicsBackend::NativeSph, NeuralBackend::NativeSpiking) => {
                BackendCapabilities {
                    gpu: false,
                    multithreaded: false,
                    simd: true,
                    max_particles: 100_000,
                    max_neurons: 10_000,
                    real_time: true,
                    deterministic: true,
                }
            }
            (PhysicsBackend::Rapier, _) => {
                BackendCapabilities {
                    gpu: false,
                    multithreaded: true,
                    simd: true,
                    max_particles: 1_000_000,
                    max_neurons: 10_000,
                    real_time: true,
                    deterministic: true,
                }
            }
            (PhysicsBackend::MuJoCo, _) => {
                BackendCapabilities {
                    gpu: true,
                    multithreaded: true,
                    simd: true,
                    max_particles: 10_000_000,
                    max_neurons: 10_000,
                    real_time: true,
                    deterministic: false, // GPU float non-determinism
                }
            }
            _ => BackendCapabilities::default(),
        };

        Self {
            physics,
            neural,
            capabilities,
        }
    }

    /// Get physics backend
    pub fn physics(&self) -> PhysicsBackend {
        self.physics
    }

    /// Get neural backend
    pub fn neural(&self) -> NeuralBackend {
        self.neural
    }

    /// Get capabilities
    pub fn capabilities(&self) -> &BackendCapabilities {
        &self.capabilities
    }

    /// Check if backend combination is supported
    pub fn is_supported(&self) -> bool {
        match (&self.physics, &self.neural) {
            (PhysicsBackend::NativeSph, NeuralBackend::NativeSpiking) => true,
            // Other combinations would require integration
            _ => false,
        }
    }

    /// Get backend name
    pub fn name(&self) -> String {
        format!("{:?}+{:?}", self.physics, self.neural)
    }
}

/// Trait for physics backend implementations
pub trait PhysicsBackendImpl: Send + Sync {
    /// Step physics simulation
    fn step(&mut self, dt: f32);

    /// Get particle positions
    fn positions(&self) -> &[f32];

    /// Get particle velocities
    fn velocities(&self) -> &[f32];

    /// Apply muscle forces
    fn apply_muscle_forces(&mut self, activations: &[f32; 96]);

    /// Get center of mass
    fn center_of_mass(&self) -> [f64; 3];

    /// Reset simulation
    fn reset(&mut self);
}

/// Trait for neural backend implementations
pub trait NeuralBackendImpl: Send + Sync {
    /// Step neural simulation
    fn step(&mut self, dt: f32);

    /// Get membrane potentials
    fn voltages(&self) -> &[f32];

    /// Set input current
    fn set_input(&mut self, neuron_id: u16, current: f32);

    /// Get spike events
    fn get_spikes(&self) -> Vec<u16>;

    /// Get muscle output
    fn get_muscle_output(&self) -> [f32; 96];

    /// Reset simulation
    fn reset(&mut self);
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_native_bridge() {
        let bridge = BackendBridge::native();
        assert!(bridge.is_supported());
        assert!(bridge.capabilities().real_time);
    }

    #[test]
    fn test_unsupported_bridge() {
        let bridge = BackendBridge::new(PhysicsBackend::MuJoCo, NeuralBackend::Brian2);
        assert!(!bridge.is_supported());
    }
}
