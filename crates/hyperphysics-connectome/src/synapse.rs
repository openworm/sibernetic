//! Synaptic Connections and Dynamics
//!
//! Models chemical synapses and electrical gap junctions.

#[cfg(feature = "serde")]
use serde::{Deserialize, Serialize};

use crate::neuron::NeuronId;

/// Synapse type
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub enum SynapseType {
    /// Chemical synapse (one-directional)
    Chemical,
    /// Electrical gap junction (bidirectional)
    GapJunction,
}

/// A synaptic connection between two neurons
#[derive(Debug, Clone)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct Synapse {
    /// Presynaptic neuron
    pub pre: NeuronId,
    /// Postsynaptic neuron
    pub post: NeuronId,
    /// Synapse type
    pub synapse_type: SynapseType,
    /// Synaptic weight (conductance in nS)
    pub weight: f32,
    /// Number of synaptic contacts (from electron microscopy)
    pub num_contacts: u32,
    /// Reversal potential for chemical synapses (mV)
    pub e_rev: f32,
}

impl Synapse {
    /// Create a new chemical synapse
    pub fn chemical(pre: NeuronId, post: NeuronId, weight: f32) -> Self {
        Self {
            pre,
            post,
            synapse_type: SynapseType::Chemical,
            weight,
            num_contacts: 1,
            e_rev: 0.0, // Excitatory default
        }
    }

    /// Create a new gap junction
    pub fn gap_junction(pre: NeuronId, post: NeuronId, weight: f32) -> Self {
        Self {
            pre,
            post,
            synapse_type: SynapseType::GapJunction,
            weight,
            num_contacts: 1,
            e_rev: 0.0, // Not used for gap junctions
        }
    }

    /// Set as excitatory (E_rev = 0 mV)
    pub fn excitatory(mut self) -> Self {
        self.e_rev = 0.0;
        self
    }

    /// Set as inhibitory (E_rev = -80 mV)
    pub fn inhibitory(mut self) -> Self {
        self.e_rev = -80.0;
        self
    }

    /// Set reversal potential
    pub fn with_reversal(mut self, e_rev: f32) -> Self {
        self.e_rev = e_rev;
        self
    }

    /// Set number of contacts
    pub fn with_contacts(mut self, num_contacts: u32) -> Self {
        self.num_contacts = num_contacts;
        self
    }

    /// Check if synapse is excitatory
    pub fn is_excitatory(&self) -> bool {
        self.e_rev > -40.0
    }

    /// Check if synapse is inhibitory
    pub fn is_inhibitory(&self) -> bool {
        self.e_rev <= -40.0
    }

    /// Check if this is a gap junction
    pub fn is_gap_junction(&self) -> bool {
        matches!(self.synapse_type, SynapseType::GapJunction)
    }
}

/// Dynamic state of a synapse during simulation
#[derive(Debug, Clone, Copy, Default)]
pub struct SynapticState {
    /// Synaptic conductance (fraction of max, 0-1)
    pub g: f32,
    /// Presynaptic calcium concentration (for plasticity)
    pub ca_pre: f32,
    /// Postsynaptic calcium concentration (for plasticity)
    pub ca_post: f32,
    /// Eligibility trace (for reward-modulated plasticity)
    pub eligibility: f32,
    /// Short-term plasticity: facilitation variable
    pub facilitation: f32,
    /// Short-term plasticity: depression variable
    pub depression: f32,
}

impl SynapticState {
    /// Create new synaptic state
    pub fn new() -> Self {
        Self {
            g: 0.0,
            ca_pre: 0.0,
            ca_post: 0.0,
            eligibility: 0.0,
            facilitation: 1.0,
            depression: 1.0,
        }
    }
}

/// Parameters for chemical synapse dynamics
#[derive(Debug, Clone, Copy)]
pub struct ChemicalSynapseParams {
    /// Rise time constant (ms)
    pub tau_rise: f32,
    /// Decay time constant (ms)
    pub tau_decay: f32,
    /// Maximum conductance (nS)
    pub g_max: f32,
    /// Delay (ms)
    pub delay: f32,
}

impl Default for ChemicalSynapseParams {
    fn default() -> Self {
        Self {
            tau_rise: 0.5,
            tau_decay: 5.0,
            g_max: 1.0,
            delay: 1.0,
        }
    }
}

impl ChemicalSynapseParams {
    /// AMPA-type (fast excitatory)
    pub fn ampa() -> Self {
        Self {
            tau_rise: 0.2,
            tau_decay: 2.0,
            g_max: 1.0,
            delay: 1.0,
        }
    }

    /// NMDA-type (slow excitatory)
    pub fn nmda() -> Self {
        Self {
            tau_rise: 2.0,
            tau_decay: 100.0,
            g_max: 0.5,
            delay: 1.0,
        }
    }

    /// GABA-A type (fast inhibitory)
    pub fn gaba_a() -> Self {
        Self {
            tau_rise: 0.5,
            tau_decay: 10.0,
            g_max: 1.0,
            delay: 1.0,
        }
    }

    /// Update conductance after presynaptic spike
    pub fn on_spike(&self, state: &mut SynapticState) {
        // Instantaneous rise approximation
        state.g += self.g_max * state.facilitation * state.depression;
    }

    /// Update conductance decay
    pub fn step(&self, state: &mut SynapticState, dt: f32) {
        // Exponential decay
        state.g *= (-dt / self.tau_decay).exp();

        // Short-term plasticity recovery
        state.facilitation += (1.0 - state.facilitation) * dt / 500.0;
        state.depression += (1.0 - state.depression) * dt / 200.0;
    }
}

/// Parameters for gap junction dynamics
#[derive(Debug, Clone, Copy)]
pub struct GapJunctionParams {
    /// Conductance (nS)
    pub g: f32,
    /// Rectification factor (1.0 = no rectification)
    pub rectification: f32,
}

impl Default for GapJunctionParams {
    fn default() -> Self {
        Self {
            g: 1.0,
            rectification: 1.0,
        }
    }
}

impl GapJunctionParams {
    /// Calculate current through gap junction
    ///
    /// I = g * (V_pre - V_post)
    pub fn current(&self, v_pre: f32, v_post: f32) -> f32 {
        let dv = v_pre - v_post;

        // Apply rectification if needed
        let g_eff = if self.rectification != 1.0 {
            if dv > 0.0 {
                self.g * self.rectification
            } else {
                self.g / self.rectification
            }
        } else {
            self.g
        };

        g_eff * dv
    }
}

/// Short-term plasticity parameters (Tsodyks-Markram model)
#[derive(Debug, Clone, Copy)]
pub struct ShortTermPlasticityParams {
    /// Use parameter (fraction released per spike)
    pub u: f32,
    /// Time constant for facilitation recovery (ms)
    pub tau_f: f32,
    /// Time constant for depression recovery (ms)
    pub tau_d: f32,
}

impl Default for ShortTermPlasticityParams {
    fn default() -> Self {
        Self {
            u: 0.5,
            tau_f: 500.0,
            tau_d: 200.0,
        }
    }
}

impl ShortTermPlasticityParams {
    /// Facilitating synapse
    pub fn facilitating() -> Self {
        Self {
            u: 0.1,
            tau_f: 1000.0,
            tau_d: 100.0,
        }
    }

    /// Depressing synapse
    pub fn depressing() -> Self {
        Self {
            u: 0.8,
            tau_f: 100.0,
            tau_d: 800.0,
        }
    }

    /// Update on presynaptic spike
    pub fn on_spike(&self, state: &mut SynapticState) {
        // Facilitation increases
        state.facilitation += self.u * (1.0 - state.facilitation);
        // Depression decreases
        state.depression *= 1.0 - state.facilitation;
    }

    /// Update recovery between spikes
    pub fn step(&self, state: &mut SynapticState, dt: f32) {
        state.facilitation += (1.0 - state.facilitation) * dt / self.tau_f;
        state.depression += (1.0 - state.depression) * dt / self.tau_d;
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_synapse_creation() {
        let syn = Synapse::chemical(0, 1, 1.5).excitatory();
        assert_eq!(syn.pre, 0);
        assert_eq!(syn.post, 1);
        assert!(syn.is_excitatory());

        let inh = Synapse::chemical(0, 1, 1.0).inhibitory();
        assert!(inh.is_inhibitory());

        let gap = Synapse::gap_junction(0, 1, 0.5);
        assert!(gap.is_gap_junction());
    }

    #[test]
    fn test_synaptic_dynamics() {
        let params = ChemicalSynapseParams::default();
        let mut state = SynapticState::new();

        // Trigger spike
        params.on_spike(&mut state);
        let g_peak = state.g;

        // Decay
        for _ in 0..100 {
            params.step(&mut state, 0.1);
        }

        assert!(state.g < g_peak, "Conductance should decay");
    }

    #[test]
    fn test_gap_junction_current() {
        let params = GapJunctionParams::default();

        // Current should flow from higher to lower potential
        let i = params.current(-50.0, -70.0);
        assert!(i > 0.0, "Current should be positive when V_pre > V_post");

        let i2 = params.current(-70.0, -50.0);
        assert!(i2 < 0.0, "Current should be negative when V_pre < V_post");
    }
}
