//! Spiking Neural Network Simulation
//!
//! Simulates the C. elegans nervous system using biophysical neuron models.

use crate::connectome::Connectome;
use crate::models::{ModelLevel, ModelParams, NeuronModel};
use crate::neuron::{NeuronId, NeuronState};
use crate::synapse::{SynapticState, ChemicalSynapseParams, GapJunctionParams, SynapseType};
use crate::muscle_map::MuscleMap;

/// Spiking neural network simulator
pub struct SpikingNetwork {
    /// Model parameters
    params: ModelParams,
    /// Neuron states
    neuron_states: Vec<NeuronState>,
    /// Synaptic states
    synapse_states: Vec<SynapticState>,
    /// Connectome reference
    connectome: Connectome,
    /// Chemical synapse parameters
    chem_params: ChemicalSynapseParams,
    /// Gap junction parameters
    gap_params: GapJunctionParams,
    /// Muscle output mapping
    muscle_map: MuscleMap,
    /// Current simulation time (ms)
    time: f64,
    /// Spike history for STDP (neuron_id, spike_time)
    spike_history: Vec<(NeuronId, f64)>,
    /// Maximum spike history length
    max_spike_history: usize,
}

impl SpikingNetwork {
    /// Create a network from a connectome
    pub fn from_connectome(connectome: &Connectome, level: ModelLevel) -> Self {
        let params = level.default_params();
        let n_neurons = connectome.num_neurons();
        let n_synapses = connectome.num_synapses();

        // Initialize neuron states at resting potential
        let neuron_states: Vec<_> = (0..n_neurons)
            .map(|_| params.neuron_model.initial_state())
            .collect();

        // Initialize synaptic states
        let synapse_states: Vec<_> = (0..n_synapses)
            .map(|_| SynapticState::new())
            .collect();

        // Build muscle map
        let muscle_map = MuscleMap::from_connectome(connectome);

        Self {
            params,
            neuron_states,
            synapse_states,
            connectome: connectome.clone(),
            chem_params: ChemicalSynapseParams::default(),
            gap_params: GapJunctionParams::default(),
            muscle_map,
            time: 0.0,
            spike_history: Vec::new(),
            max_spike_history: 10000,
        }
    }

    /// Create for C. elegans with specified model level
    pub fn celegans(level: ModelLevel) -> Self {
        let connectome = Connectome::celegans();
        Self::from_connectome(&connectome, level)
    }

    /// Perform one simulation step
    pub fn step(&mut self, dt: f32) {
        let n_neurons = self.neuron_states.len();
        let n_synapses = self.synapse_states.len();

        // 1. Update synaptic conductances (decay)
        for i in 0..n_synapses {
            self.chem_params.step(&mut self.synapse_states[i], dt);
        }

        // 2. Compute synaptic currents
        self.compute_synaptic_currents(dt);

        // 3. Update neuron states
        let mut new_spikes: Vec<NeuronId> = Vec::new();

        for i in 0..n_neurons {
            let state = &mut self.neuron_states[i];

            // Update using appropriate model
            self.params.neuron_model.step(state, dt);

            // Record spikes
            if state.spiked {
                new_spikes.push(i as NeuronId);
            }

            // Post-step cleanup
            state.post_step(dt);
        }

        // 4. Process new spikes
        for &neuron_id in &new_spikes {
            self.on_spike(neuron_id);
        }

        // 5. Update time
        self.time += dt as f64;
    }

    /// Run simulation for specified duration
    pub fn run(&mut self, duration_ms: f64) {
        let dt = self.params.level.recommended_dt();
        let n_steps = (duration_ms / dt as f64).ceil() as usize;

        for _ in 0..n_steps {
            self.step(dt);
        }
    }

    /// Compute synaptic currents from all synapses
    fn compute_synaptic_currents(&mut self, _dt: f32) {
        let synapses = self.connectome.synapses();

        for (i, synapse) in synapses.iter().enumerate() {
            let pre_v = self.neuron_states[synapse.pre as usize].v;
            let post_v = self.neuron_states[synapse.post as usize].v;

            let current = match synapse.synapse_type {
                SynapseType::Chemical => {
                    // I = g * (V - E_rev)
                    let g = self.synapse_states[i].g * synapse.weight * self.params.weight_scale;
                    g * (post_v - synapse.e_rev)
                }
                SynapseType::GapJunction => {
                    // I = g * (V_pre - V_post)
                    let g = synapse.weight * self.params.gap_junction_g;
                    self.gap_params.current(pre_v, post_v) * g
                }
            };

            // Add to postsynaptic neuron
            self.neuron_states[synapse.post as usize].add_input(-current);

            // For gap junctions, also add to presynaptic (bidirectional)
            if synapse.synapse_type == SynapseType::GapJunction {
                self.neuron_states[synapse.pre as usize].add_input(current);
            }
        }
    }

    /// Handle spike event
    fn on_spike(&mut self, neuron_id: NeuronId) {
        // Record spike
        if self.spike_history.len() >= self.max_spike_history {
            self.spike_history.remove(0);
        }
        self.spike_history.push((neuron_id, self.time));

        // Activate outgoing synapses
        let synapses = self.connectome.synapses();
        for (i, synapse) in synapses.iter().enumerate() {
            if synapse.pre == neuron_id && synapse.synapse_type == SynapseType::Chemical {
                self.chem_params.on_spike(&mut self.synapse_states[i]);
            }
        }
    }

    /// Set external input to a neuron
    pub fn set_input(&mut self, neuron_id: NeuronId, current: f32) {
        if (neuron_id as usize) < self.neuron_states.len() {
            self.neuron_states[neuron_id as usize].set_external(current);
        }
    }

    /// Set external input by neuron name
    pub fn set_input_by_name(&mut self, name: &str, current: f32) {
        if let Some(id) = self.connectome.get_id(name) {
            self.set_input(id, current);
        }
    }

    /// Clear all external inputs
    pub fn clear_inputs(&mut self) {
        for state in &mut self.neuron_states {
            state.set_external(0.0);
        }
    }

    /// Get membrane potential of a neuron
    pub fn get_voltage(&self, neuron_id: NeuronId) -> f32 {
        self.neuron_states.get(neuron_id as usize)
            .map(|s| s.v)
            .unwrap_or(0.0)
    }

    /// Get membrane potential by neuron name
    pub fn get_voltage_by_name(&self, name: &str) -> Option<f32> {
        self.connectome.get_id(name).map(|id| self.get_voltage(id))
    }

    /// Get all membrane potentials
    pub fn get_voltages(&self) -> Vec<f32> {
        self.neuron_states.iter().map(|s| s.v).collect()
    }

    /// Check if neuron spiked this step
    pub fn did_spike(&self, neuron_id: NeuronId) -> bool {
        self.neuron_states.get(neuron_id as usize)
            .map(|s| s.spiked)
            .unwrap_or(false)
    }

    /// Get neurons that spiked this step
    pub fn get_spikes(&self) -> Vec<NeuronId> {
        self.neuron_states.iter().enumerate()
            .filter(|(_, s)| s.spiked)
            .map(|(i, _)| i as NeuronId)
            .collect()
    }

    /// Get muscle output (activation pattern)
    pub fn get_muscle_output(&self) -> [f32; 96] {
        self.muscle_map.compute_activation(&self.neuron_states)
    }

    /// Get spike history
    pub fn spike_history(&self) -> &[(NeuronId, f64)] {
        &self.spike_history
    }

    /// Get current time
    pub fn time(&self) -> f64 {
        self.time
    }

    /// Reset simulation
    pub fn reset(&mut self) {
        for state in &mut self.neuron_states {
            *state = self.params.neuron_model.initial_state();
        }
        for state in &mut self.synapse_states {
            *state = SynapticState::new();
        }
        self.spike_history.clear();
        self.time = 0.0;
    }

    /// Get network state for serialization
    pub fn state(&self) -> NetworkState {
        NetworkState {
            neuron_states: self.neuron_states.clone(),
            synapse_states: self.synapse_states.clone(),
            time: self.time,
        }
    }

    /// Restore network state
    pub fn set_state(&mut self, state: &NetworkState) {
        self.neuron_states = state.neuron_states.clone();
        self.synapse_states = state.synapse_states.clone();
        self.time = state.time;
    }

    /// Get number of neurons
    pub fn num_neurons(&self) -> usize {
        self.neuron_states.len()
    }

    /// Get number of synapses
    pub fn num_synapses(&self) -> usize {
        self.synapse_states.len()
    }

    /// Get connectome reference
    pub fn connectome(&self) -> &Connectome {
        &self.connectome
    }
}

/// Serializable network state
#[derive(Debug, Clone)]
pub struct NetworkState {
    /// Neuron states
    pub neuron_states: Vec<NeuronState>,
    /// Synaptic states
    pub synapse_states: Vec<SynapticState>,
    /// Simulation time
    pub time: f64,
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_network_creation() {
        let network = SpikingNetwork::celegans(ModelLevel::B);
        assert!(network.num_neurons() > 0);
        assert!(network.num_synapses() > 0);
    }

    #[test]
    fn test_network_simulation() {
        let mut network = SpikingNetwork::celegans(ModelLevel::A);

        // Apply input to sensory neurons
        network.set_input_by_name("ASHL", 10.0);

        // Run for 100ms
        network.run(100.0);

        assert!(network.time() > 0.0);

        // Should have some spikes
        assert!(!network.spike_history().is_empty());
    }

    #[test]
    fn test_muscle_output() {
        let mut network = SpikingNetwork::celegans(ModelLevel::B);

        // Stimulate command interneurons
        network.set_input_by_name("AVAL", 15.0);
        network.set_input_by_name("AVAR", 15.0);

        network.run(50.0);

        let muscles = network.get_muscle_output();

        // Some muscles should be activated
        let total: f32 = muscles.iter().sum();
        assert!(total > 0.0, "Should have muscle activation");
    }

    #[test]
    fn test_reset() {
        let mut network = SpikingNetwork::celegans(ModelLevel::A);

        network.set_input_by_name("AVAL", 10.0);
        network.run(50.0);

        let time_before = network.time();
        network.reset();

        assert_eq!(network.time(), 0.0);
        assert!(network.spike_history().is_empty());
    }
}
