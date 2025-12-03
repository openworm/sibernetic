//! Neuron Types and Models
//!
//! Defines neuron structures and multiple biophysical models from simple
//! integrate-and-fire to detailed Hodgkin-Huxley.

#[cfg(feature = "serde")]
use serde::{Deserialize, Serialize};

/// Unique identifier for a neuron
pub type NeuronId = u32;

/// Neuron functional class
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub enum NeuronClass {
    /// Sensory neuron - receives external input
    Sensory,
    /// Interneuron - processes information
    Interneuron,
    /// Motor neuron - controls muscles
    Motor,
}

impl NeuronClass {
    /// Check if this neuron can receive sensory input
    pub fn is_sensory(&self) -> bool {
        matches!(self, Self::Sensory)
    }

    /// Check if this neuron controls muscles
    pub fn is_motor(&self) -> bool {
        matches!(self, Self::Motor)
    }
}

/// Neurotransmitter type
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub enum Neurotransmitter {
    /// Acetylcholine (excitatory)
    Acetylcholine,
    /// GABA (inhibitory)
    GABA,
    /// Glutamate (excitatory)
    Glutamate,
    /// Dopamine (modulatory)
    Dopamine,
    /// Serotonin (modulatory)
    Serotonin,
    /// Unknown/unspecified
    Unknown,
}

impl Neurotransmitter {
    /// Check if this neurotransmitter is typically excitatory
    pub fn is_excitatory(&self) -> bool {
        matches!(self, Self::Acetylcholine | Self::Glutamate)
    }

    /// Check if this neurotransmitter is typically inhibitory
    pub fn is_inhibitory(&self) -> bool {
        matches!(self, Self::GABA)
    }
}

/// A neuron in the connectome
#[derive(Debug, Clone)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct Neuron {
    /// Unique identifier
    pub id: NeuronId,
    /// Neuron name (e.g., "AVAL", "VB01")
    pub name: String,
    /// Functional class
    pub class: NeuronClass,
    /// 3D position in the worm body (μm)
    pub position: [f32; 3],
    /// Primary neurotransmitter
    pub neurotransmitter: Neurotransmitter,
    /// Bilateral pair partner (if any)
    pub partner: Option<NeuronId>,
}

impl Neuron {
    /// Create a new neuron
    pub fn new(id: NeuronId, name: &str, class: NeuronClass) -> Self {
        Self {
            id,
            name: name.to_string(),
            class,
            position: [0.0, 0.0, 0.0],
            neurotransmitter: Neurotransmitter::Unknown,
            partner: None,
        }
    }

    /// Set position
    pub fn with_position(mut self, position: [f32; 3]) -> Self {
        self.position = position;
        self
    }

    /// Set neurotransmitter
    pub fn with_neurotransmitter(mut self, nt: Neurotransmitter) -> Self {
        self.neurotransmitter = nt;
        self
    }

    /// Set bilateral partner
    pub fn with_partner(mut self, partner: NeuronId) -> Self {
        self.partner = Some(partner);
        self
    }

    /// Check if this is a left-side neuron (name ends with 'L')
    pub fn is_left(&self) -> bool {
        self.name.ends_with('L')
    }

    /// Check if this is a right-side neuron (name ends with 'R')
    pub fn is_right(&self) -> bool {
        self.name.ends_with('R')
    }
}

/// Dynamic state of a neuron during simulation
#[derive(Debug, Clone, Copy, Default)]
pub struct NeuronState {
    /// Membrane potential (mV)
    pub v: f32,
    /// Recovery variable (for Izhikevich model)
    pub u: f32,
    /// Sodium channel activation (for HH model)
    pub m: f32,
    /// Sodium channel inactivation (for HH model)
    pub h: f32,
    /// Potassium channel activation (for HH model)
    pub n: f32,
    /// Calcium concentration (mM)
    pub ca: f32,
    /// Total synaptic input current (nA)
    pub i_syn: f32,
    /// External input current (nA)
    pub i_ext: f32,
    /// Time since last spike (ms)
    pub t_since_spike: f32,
    /// Did neuron spike this timestep?
    pub spiked: bool,
}

impl NeuronState {
    /// Create state at resting potential
    pub fn resting(v_rest: f32) -> Self {
        Self {
            v: v_rest,
            u: 0.0,
            m: 0.05,
            h: 0.6,
            n: 0.32,
            ca: 0.0001, // 100 nM
            i_syn: 0.0,
            i_ext: 0.0,
            t_since_spike: f32::INFINITY,
            spiked: false,
        }
    }

    /// Reset spike flag and update time since spike
    pub fn post_step(&mut self, dt: f32) {
        if self.spiked {
            self.t_since_spike = 0.0;
        } else {
            self.t_since_spike += dt;
        }
        self.spiked = false;
        self.i_syn = 0.0; // Reset synaptic input for next step
    }

    /// Add synaptic input
    pub fn add_input(&mut self, current: f32) {
        self.i_syn += current;
    }

    /// Set external input
    pub fn set_external(&mut self, current: f32) {
        self.i_ext = current;
    }

    /// Get total input current
    pub fn total_current(&self) -> f32 {
        self.i_syn + self.i_ext
    }
}

/// Leaky Integrate-and-Fire neuron model (Level A)
#[derive(Debug, Clone, Copy)]
pub struct LIFParams {
    /// Membrane time constant (ms)
    pub tau_m: f32,
    /// Resting potential (mV)
    pub v_rest: f32,
    /// Threshold potential (mV)
    pub v_thresh: f32,
    /// Reset potential (mV)
    pub v_reset: f32,
    /// Membrane resistance (MΩ)
    pub r_m: f32,
    /// Refractory period (ms)
    pub t_ref: f32,
}

impl Default for LIFParams {
    fn default() -> Self {
        Self {
            tau_m: 10.0,      // 10 ms
            v_rest: -65.0,    // -65 mV
            v_thresh: -50.0,  // -50 mV
            v_reset: -70.0,   // -70 mV
            r_m: 10.0,        // 10 MΩ
            t_ref: 2.0,       // 2 ms
        }
    }
}

impl LIFParams {
    /// Update state using forward Euler
    pub fn step(&self, state: &mut NeuronState, dt: f32) {
        // Check refractory period
        if state.t_since_spike < self.t_ref {
            state.v = self.v_reset;
            return;
        }

        // Leaky integration: τ_m * dV/dt = -(V - V_rest) + R_m * I
        let i_total = state.total_current();
        let dv = (-(state.v - self.v_rest) + self.r_m * i_total) / self.tau_m;
        state.v += dv * dt;

        // Spike detection
        if state.v >= self.v_thresh {
            state.v = self.v_reset;
            state.spiked = true;
        }
    }
}

/// Izhikevich neuron model (Level B)
#[derive(Debug, Clone, Copy)]
pub struct IzhikevichParams {
    /// Time scale of recovery variable
    pub a: f32,
    /// Sensitivity of recovery to subthreshold oscillations
    pub b: f32,
    /// After-spike reset of membrane potential
    pub c: f32,
    /// After-spike reset of recovery variable
    pub d: f32,
    /// Spike threshold
    pub v_peak: f32,
}

impl Default for IzhikevichParams {
    fn default() -> Self {
        // Regular spiking neuron parameters
        Self {
            a: 0.02,
            b: 0.2,
            c: -65.0,
            d: 8.0,
            v_peak: 30.0,
        }
    }
}

impl IzhikevichParams {
    /// Fast spiking neuron
    pub fn fast_spiking() -> Self {
        Self {
            a: 0.1,
            b: 0.2,
            c: -65.0,
            d: 2.0,
            v_peak: 30.0,
        }
    }

    /// Chattering neuron
    pub fn chattering() -> Self {
        Self {
            a: 0.02,
            b: 0.2,
            c: -50.0,
            d: 2.0,
            v_peak: 30.0,
        }
    }

    /// Low-threshold spiking
    pub fn lts() -> Self {
        Self {
            a: 0.02,
            b: 0.25,
            c: -65.0,
            d: 2.0,
            v_peak: 30.0,
        }
    }

    /// Update state using forward Euler
    pub fn step(&self, state: &mut NeuronState, dt: f32) {
        let i_total = state.total_current();

        // Izhikevich equations:
        // dv/dt = 0.04*v² + 5*v + 140 - u + I
        // du/dt = a*(b*v - u)
        let dv = 0.04 * state.v * state.v + 5.0 * state.v + 140.0 - state.u + i_total;
        let du = self.a * (self.b * state.v - state.u);

        state.v += dv * dt;
        state.u += du * dt;

        // Spike detection and reset
        if state.v >= self.v_peak {
            state.v = self.c;
            state.u += self.d;
            state.spiked = true;
        }
    }
}

/// Hodgkin-Huxley neuron model (Level D)
#[derive(Debug, Clone, Copy)]
pub struct HHParams {
    /// Membrane capacitance (μF/cm²)
    pub c_m: f32,
    /// Sodium conductance (mS/cm²)
    pub g_na: f32,
    /// Potassium conductance (mS/cm²)
    pub g_k: f32,
    /// Leak conductance (mS/cm²)
    pub g_l: f32,
    /// Sodium reversal potential (mV)
    pub e_na: f32,
    /// Potassium reversal potential (mV)
    pub e_k: f32,
    /// Leak reversal potential (mV)
    pub e_l: f32,
    /// Spike threshold for detection
    pub v_thresh: f32,
}

impl Default for HHParams {
    fn default() -> Self {
        Self {
            c_m: 1.0,        // 1 μF/cm²
            g_na: 120.0,     // 120 mS/cm²
            g_k: 36.0,       // 36 mS/cm²
            g_l: 0.3,        // 0.3 mS/cm²
            e_na: 50.0,      // 50 mV
            e_k: -77.0,      // -77 mV
            e_l: -54.387,    // -54.387 mV
            v_thresh: 0.0,   // 0 mV
        }
    }
}

impl HHParams {
    /// Update state using forward Euler
    pub fn step(&self, state: &mut NeuronState, dt: f32) {
        let v = state.v;
        let m = state.m;
        let h = state.h;
        let n = state.n;

        // Rate functions
        let alpha_m = 0.1 * (v + 40.0) / (1.0 - (-0.1 * (v + 40.0)).exp());
        let beta_m = 4.0 * (-0.0556 * (v + 65.0)).exp();
        let alpha_h = 0.07 * (-0.05 * (v + 65.0)).exp();
        let beta_h = 1.0 / (1.0 + (-0.1 * (v + 35.0)).exp());
        let alpha_n = 0.01 * (v + 55.0) / (1.0 - (-0.1 * (v + 55.0)).exp());
        let beta_n = 0.125 * (-0.0125 * (v + 65.0)).exp();

        // Gating variable derivatives
        let dm = alpha_m * (1.0 - m) - beta_m * m;
        let dh = alpha_h * (1.0 - h) - beta_h * h;
        let dn = alpha_n * (1.0 - n) - beta_n * n;

        // Ionic currents
        let i_na = self.g_na * m * m * m * h * (v - self.e_na);
        let i_k = self.g_k * n * n * n * n * (v - self.e_k);
        let i_l = self.g_l * (v - self.e_l);

        // Membrane potential derivative
        let i_total = state.total_current();
        let dv = (i_total - i_na - i_k - i_l) / self.c_m;

        // Update state
        state.v += dv * dt;
        state.m += dm * dt;
        state.h += dh * dt;
        state.n += dn * dt;

        // Clamp gating variables to [0, 1]
        state.m = state.m.clamp(0.0, 1.0);
        state.h = state.h.clamp(0.0, 1.0);
        state.n = state.n.clamp(0.0, 1.0);

        // Spike detection (crossing threshold from below)
        let was_below = v < self.v_thresh;
        let is_above = state.v >= self.v_thresh;
        if was_below && is_above {
            state.spiked = true;
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_neuron_creation() {
        let neuron = Neuron::new(0, "AVAL", NeuronClass::Interneuron)
            .with_position([100.0, 50.0, 200.0])
            .with_neurotransmitter(Neurotransmitter::Acetylcholine);

        assert_eq!(neuron.name, "AVAL");
        assert_eq!(neuron.class, NeuronClass::Interneuron);
        assert!(neuron.is_left());
    }

    #[test]
    fn test_lif_spiking() {
        let params = LIFParams::default();
        let mut state = NeuronState::resting(params.v_rest);

        // Apply strong input current
        state.set_external(5.0); // 5 nA

        // Should eventually spike
        for _ in 0..1000 {
            params.step(&mut state, 0.1);
            state.post_step(0.1);
            if state.t_since_spike < 1.0 {
                return; // Test passed
            }
        }

        panic!("LIF neuron should have spiked");
    }

    #[test]
    fn test_izhikevich_spiking() {
        let params = IzhikevichParams::default();
        let mut state = NeuronState::resting(-65.0);

        state.set_external(10.0);

        for _ in 0..1000 {
            params.step(&mut state, 0.1);
            state.post_step(0.1);
            if state.t_since_spike < 1.0 {
                return;
            }
        }

        panic!("Izhikevich neuron should have spiked");
    }
}
