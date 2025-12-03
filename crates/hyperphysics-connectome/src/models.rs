//! Model Levels and Parameters
//!
//! Defines the different abstraction levels from c302 (A through D)
//! and their associated parameter sets.

#[cfg(feature = "serde")]
use serde::{Deserialize, Serialize};

use crate::neuron::{LIFParams, IzhikevichParams, HHParams, NeuronState};
use crate::synapse::ChemicalSynapseParams;

/// Model complexity level (from c302)
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub enum ModelLevel {
    /// Level A: Simple integrate-and-fire
    /// - Fast simulation (~1ms timescale)
    /// - Good for behavioral studies
    /// - Captures basic spike timing
    A,

    /// Level B: Izhikevich neurons
    /// - Medium complexity
    /// - Captures various spiking patterns
    /// - Good balance of speed and realism
    B,

    /// Level C: Conductance-based (simplified)
    /// - Includes basic ion channels
    /// - More realistic dynamics
    /// - Slower than A/B
    C,

    /// Level D: Full Hodgkin-Huxley
    /// - Biophysically detailed
    /// - All major ion channels
    /// - Slowest but most accurate
    D,
}

impl ModelLevel {
    /// Get computational cost multiplier relative to Level A
    pub fn cost_factor(&self) -> f32 {
        match self {
            Self::A => 1.0,
            Self::B => 5.0,
            Self::C => 20.0,
            Self::D => 100.0,
        }
    }

    /// Get recommended time step (ms)
    pub fn recommended_dt(&self) -> f32 {
        match self {
            Self::A => 0.5,
            Self::B => 0.1,
            Self::C => 0.05,
            Self::D => 0.01,
        }
    }

    /// Get default parameters for this level
    pub fn default_params(&self) -> ModelParams {
        match self {
            Self::A => ModelParams::level_a(),
            Self::B => ModelParams::level_b(),
            Self::C => ModelParams::level_c(),
            Self::D => ModelParams::level_d(),
        }
    }
}

impl Default for ModelLevel {
    fn default() -> Self {
        Self::B
    }
}

/// Complete parameter set for a model level
#[derive(Debug, Clone)]
pub struct ModelParams {
    /// Model level
    pub level: ModelLevel,
    /// Neuron model implementation
    pub neuron_model: NeuronModel,
    /// Synapse parameters
    pub synapse: ChemicalSynapseParams,
    /// Default synaptic weight scaling
    pub weight_scale: f32,
    /// Default gap junction conductance (nS)
    pub gap_junction_g: f32,
}

impl ModelParams {
    /// Level A parameters (integrate-and-fire)
    pub fn level_a() -> Self {
        Self {
            level: ModelLevel::A,
            neuron_model: NeuronModel::LIF(LIFParams::default()),
            synapse: ChemicalSynapseParams {
                tau_rise: 0.5,
                tau_decay: 5.0,
                g_max: 1.0,
                delay: 1.0,
            },
            weight_scale: 1.0,
            gap_junction_g: 0.5,
        }
    }

    /// Level B parameters (Izhikevich)
    pub fn level_b() -> Self {
        Self {
            level: ModelLevel::B,
            neuron_model: NeuronModel::Izhikevich(IzhikevichParams::default()),
            synapse: ChemicalSynapseParams {
                tau_rise: 0.2,
                tau_decay: 3.0,
                g_max: 2.0,
                delay: 0.5,
            },
            weight_scale: 5.0,
            gap_junction_g: 1.0,
        }
    }

    /// Level C parameters (conductance-based)
    pub fn level_c() -> Self {
        Self {
            level: ModelLevel::C,
            neuron_model: NeuronModel::HH(HHParams::default()),
            synapse: ChemicalSynapseParams::ampa(),
            weight_scale: 0.1,
            gap_junction_g: 0.1,
        }
    }

    /// Level D parameters (full Hodgkin-Huxley)
    pub fn level_d() -> Self {
        Self {
            level: ModelLevel::D,
            neuron_model: NeuronModel::HH(HHParams::default()),
            synapse: ChemicalSynapseParams {
                tau_rise: 0.1,
                tau_decay: 2.0,
                g_max: 0.5,
                delay: 0.5,
            },
            weight_scale: 0.05,
            gap_junction_g: 0.05,
        }
    }
}

/// Neuron model variants
#[derive(Debug, Clone)]
pub enum NeuronModel {
    /// Leaky integrate-and-fire
    LIF(LIFParams),
    /// Izhikevich model
    Izhikevich(IzhikevichParams),
    /// Hodgkin-Huxley model
    HH(HHParams),
}

impl NeuronModel {
    /// Get resting potential
    pub fn v_rest(&self) -> f32 {
        match self {
            Self::LIF(p) => p.v_rest,
            Self::Izhikevich(_) => -65.0,
            Self::HH(p) => p.e_l,
        }
    }

    /// Get spike threshold
    pub fn v_thresh(&self) -> f32 {
        match self {
            Self::LIF(p) => p.v_thresh,
            Self::Izhikevich(p) => p.v_peak,
            Self::HH(p) => p.v_thresh,
        }
    }

    /// Update neuron state for one timestep
    pub fn step(&self, state: &mut NeuronState, dt: f32) {
        match self {
            Self::LIF(p) => p.step(state, dt),
            Self::Izhikevich(p) => p.step(state, dt),
            Self::HH(p) => p.step(state, dt),
        }
    }

    /// Create initial state
    pub fn initial_state(&self) -> NeuronState {
        NeuronState::resting(self.v_rest())
    }
}

impl Default for NeuronModel {
    fn default() -> Self {
        Self::Izhikevich(IzhikevichParams::default())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_model_levels() {
        assert!(ModelLevel::A.cost_factor() < ModelLevel::D.cost_factor());
        assert!(ModelLevel::D.recommended_dt() < ModelLevel::A.recommended_dt());
    }

    #[test]
    fn test_model_params() {
        let params_a = ModelParams::level_a();
        let params_d = ModelParams::level_d();

        assert!(matches!(params_a.neuron_model, NeuronModel::LIF(_)));
        assert!(matches!(params_d.neuron_model, NeuronModel::HH(_)));
    }

    #[test]
    fn test_neuron_model_step() {
        let model = NeuronModel::Izhikevich(IzhikevichParams::default());
        let mut state = model.initial_state();

        state.set_external(10.0);

        for _ in 0..100 {
            model.step(&mut state, 0.1);
        }

        // Voltage should have changed
        assert!((state.v - model.v_rest()).abs() > 1.0);
    }
}
