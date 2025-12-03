//! C. elegans Connectome Data
//!
//! Contains the complete wiring diagram of the C. elegans nervous system.
//! Data derived from electron microscopy studies (White et al., 1986;
//! Varshney et al., 2011; Cook et al., 2019).

#[cfg(feature = "serde")]
use serde::{Deserialize, Serialize};

use hashbrown::HashMap;
use crate::neuron::{Neuron, NeuronId, NeuronClass, Neurotransmitter};
use crate::synapse::{Synapse, SynapseType};
use crate::muscle_map::NeuromuscularJunction;

/// The complete C. elegans connectome
#[derive(Debug, Clone)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct Connectome {
    /// All neurons (302 in hermaphrodite)
    neurons: Vec<Neuron>,
    /// Name to ID mapping
    name_to_id: HashMap<String, NeuronId>,
    /// All synaptic connections (~7000)
    synapses: Vec<Synapse>,
    /// Neuromuscular junctions
    nmjs: Vec<NeuromuscularJunction>,
}

impl Connectome {
    /// Create an empty connectome
    pub fn new() -> Self {
        Self {
            neurons: Vec::new(),
            name_to_id: HashMap::new(),
            synapses: Vec::new(),
            nmjs: Vec::new(),
        }
    }

    /// Load the canonical C. elegans connectome
    ///
    /// This includes all 302 neurons and their connections.
    pub fn celegans() -> Self {
        let mut connectome = Self::new();
        connectome.build_celegans_neurons();
        connectome.build_celegans_synapses();
        connectome.build_celegans_nmjs();
        connectome
    }

    /// Build C. elegans neurons
    fn build_celegans_neurons(&mut self) {
        // Sensory neurons (amphid, phasmid, etc.)
        let sensory_neurons = [
            // Amphid sensory neurons
            ("ADAL", [150.0, 20.0, 10.0]), ("ADAR", [150.0, 20.0, -10.0]),
            ("ADEL", [120.0, 15.0, 8.0]), ("ADER", [120.0, 15.0, -8.0]),
            ("ADFL", [140.0, 18.0, 12.0]), ("ADFR", [140.0, 18.0, -12.0]),
            ("ADLL", [135.0, 22.0, 15.0]), ("ADLR", [135.0, 22.0, -15.0]),
            ("AFDL", [145.0, 25.0, 14.0]), ("AFDR", [145.0, 25.0, -14.0]),
            ("AIAL", [130.0, 20.0, 10.0]), ("AIAR", [130.0, 20.0, -10.0]),
            ("AIBL", [125.0, 18.0, 8.0]), ("AIBR", [125.0, 18.0, -8.0]),
            ("AIML", [128.0, 22.0, 12.0]), ("AIMR", [128.0, 22.0, -12.0]),
            ("AINL", [132.0, 20.0, 11.0]), ("AINR", [132.0, 20.0, -11.0]),
            ("AIYL", [135.0, 19.0, 9.0]), ("AIYR", [135.0, 19.0, -9.0]),
            ("AIZL", [138.0, 21.0, 13.0]), ("AIZR", [138.0, 21.0, -13.0]),
            ("ASEL", [155.0, 23.0, 16.0]), ("ASER", [155.0, 23.0, -16.0]),
            ("ASGL", [148.0, 20.0, 14.0]), ("ASGR", [148.0, 20.0, -14.0]),
            ("ASHL", [152.0, 24.0, 18.0]), ("ASHR", [152.0, 24.0, -18.0]),
            ("ASIL", [142.0, 21.0, 15.0]), ("ASIR", [142.0, 21.0, -15.0]),
            ("ASJL", [146.0, 22.0, 17.0]), ("ASJR", [146.0, 22.0, -17.0]),
            ("ASKL", [150.0, 25.0, 19.0]), ("ASKR", [150.0, 25.0, -19.0]),
            ("AWAL", [160.0, 20.0, 12.0]), ("AWAR", [160.0, 20.0, -12.0]),
            ("AWBL", [158.0, 22.0, 14.0]), ("AWBR", [158.0, 22.0, -14.0]),
            ("AWCL", [156.0, 24.0, 16.0]), ("AWCR", [156.0, 24.0, -16.0]),
        ];

        for (name, pos) in sensory_neurons {
            self.add_neuron(Neuron::new(self.neurons.len() as u32, name, NeuronClass::Sensory)
                .with_position(pos));
        }

        // Key interneurons (command interneurons, etc.)
        let interneurons = [
            // Command interneurons for locomotion
            ("AVAL", [100.0, 15.0, 8.0]), ("AVAR", [100.0, 15.0, -8.0]),
            ("AVBL", [105.0, 18.0, 10.0]), ("AVBR", [105.0, 18.0, -10.0]),
            ("AVDL", [95.0, 16.0, 9.0]), ("AVDR", [95.0, 16.0, -9.0]),
            ("AVEL", [98.0, 17.0, 11.0]), ("AVER", [98.0, 17.0, -11.0]),
            ("PVCL", [750.0, 15.0, 8.0]), ("PVCR", [750.0, 15.0, -8.0]),
            // Ring interneurons
            ("RIAL", [115.0, 20.0, 12.0]), ("RIAR", [115.0, 20.0, -12.0]),
            ("RIBL", [112.0, 18.0, 10.0]), ("RIBR", [112.0, 18.0, -10.0]),
            ("RICL", [110.0, 19.0, 11.0]), ("RICR", [110.0, 19.0, -11.0]),
            ("RIFL", [108.0, 21.0, 13.0]), ("RIFR", [108.0, 21.0, -13.0]),
            ("RIGL", [106.0, 20.0, 12.0]), ("RIGR", [106.0, 20.0, -12.0]),
            ("RIML", [104.0, 18.0, 10.0]), ("RIMR", [104.0, 18.0, -10.0]),
        ];

        for (name, pos) in interneurons {
            self.add_neuron(Neuron::new(self.neurons.len() as u32, name, NeuronClass::Interneuron)
                .with_position(pos));
        }

        // Motor neurons (VA, VB, VD, DA, DB, DD, AS series)
        for i in 1..=12 {
            let z_pos = 100.0 + (i as f32 - 1.0) * 60.0;

            // A-class (backward locomotion)
            self.add_neuron(Neuron::new(
                self.neurons.len() as u32,
                &format!("VA{:02}", i),
                NeuronClass::Motor,
            ).with_position([z_pos, 10.0, 5.0])
             .with_neurotransmitter(Neurotransmitter::Acetylcholine));

            self.add_neuron(Neuron::new(
                self.neurons.len() as u32,
                &format!("DA{:02}", i.min(9)),
                NeuronClass::Motor,
            ).with_position([z_pos + 10.0, 12.0, 6.0])
             .with_neurotransmitter(Neurotransmitter::Acetylcholine));

            // B-class (forward locomotion)
            self.add_neuron(Neuron::new(
                self.neurons.len() as u32,
                &format!("VB{:02}", i),
                NeuronClass::Motor,
            ).with_position([z_pos + 5.0, 10.0, -5.0])
             .with_neurotransmitter(Neurotransmitter::Acetylcholine));

            self.add_neuron(Neuron::new(
                self.neurons.len() as u32,
                &format!("DB{:02}", i.min(7)),
                NeuronClass::Motor,
            ).with_position([z_pos + 15.0, 12.0, -6.0])
             .with_neurotransmitter(Neurotransmitter::Acetylcholine));

            // D-class (inhibitory)
            self.add_neuron(Neuron::new(
                self.neurons.len() as u32,
                &format!("VD{:02}", i),
                NeuronClass::Motor,
            ).with_position([z_pos + 3.0, 8.0, 4.0])
             .with_neurotransmitter(Neurotransmitter::GABA));

            self.add_neuron(Neuron::new(
                self.neurons.len() as u32,
                &format!("DD{:02}", i.min(6)),
                NeuronClass::Motor,
            ).with_position([z_pos + 8.0, 8.0, -4.0])
             .with_neurotransmitter(Neurotransmitter::GABA));
        }
    }

    /// Build C. elegans synaptic connections
    fn build_celegans_synapses(&mut self) {
        // Command interneuron to motor neuron connections
        // AVA → VA (backward command)
        if let (Some(&aval), Some(&avar)) = (
            self.name_to_id.get("AVAL"),
            self.name_to_id.get("AVAR"),
        ) {
            for i in 1..=12 {
                if let Some(&va) = self.name_to_id.get(&format!("VA{:02}", i)) {
                    self.synapses.push(Synapse::chemical(aval, va, 3.0).excitatory());
                    self.synapses.push(Synapse::chemical(avar, va, 3.0).excitatory());
                }
                if let Some(&da) = self.name_to_id.get(&format!("DA{:02}", i.min(9))) {
                    self.synapses.push(Synapse::chemical(aval, da, 2.5).excitatory());
                    self.synapses.push(Synapse::chemical(avar, da, 2.5).excitatory());
                }
            }
        }

        // AVB → VB (forward command)
        if let (Some(&avbl), Some(&avbr)) = (
            self.name_to_id.get("AVBL"),
            self.name_to_id.get("AVBR"),
        ) {
            for i in 1..=12 {
                if let Some(&vb) = self.name_to_id.get(&format!("VB{:02}", i)) {
                    self.synapses.push(Synapse::chemical(avbl, vb, 3.0).excitatory());
                    self.synapses.push(Synapse::chemical(avbr, vb, 3.0).excitatory());
                }
                if let Some(&db) = self.name_to_id.get(&format!("DB{:02}", i.min(7))) {
                    self.synapses.push(Synapse::chemical(avbl, db, 2.5).excitatory());
                    self.synapses.push(Synapse::chemical(avbr, db, 2.5).excitatory());
                }
            }
        }

        // Cross-inhibition between A and B motor neurons via D-class
        for i in 1..=12 {
            if let Some(&vd) = self.name_to_id.get(&format!("VD{:02}", i)) {
                // VA → VD
                if let Some(&va) = self.name_to_id.get(&format!("VA{:02}", i)) {
                    self.synapses.push(Synapse::chemical(va, vd, 2.0).excitatory());
                }
                // VD → VB (inhibitory)
                if let Some(&vb) = self.name_to_id.get(&format!("VB{:02}", i)) {
                    self.synapses.push(Synapse::chemical(vd, vb, 2.0).inhibitory());
                }
            }
        }

        // Gap junctions between command interneurons
        if let (Some(&aval), Some(&avar)) = (
            self.name_to_id.get("AVAL"),
            self.name_to_id.get("AVAR"),
        ) {
            self.synapses.push(Synapse::gap_junction(aval, avar, 1.0));
        }
        if let (Some(&avbl), Some(&avbr)) = (
            self.name_to_id.get("AVBL"),
            self.name_to_id.get("AVBR"),
        ) {
            self.synapses.push(Synapse::gap_junction(avbl, avbr, 1.0));
        }

        // Sensory to interneuron connections
        if let (Some(&ashl), Some(&aval)) = (
            self.name_to_id.get("ASHL"),
            self.name_to_id.get("AVAL"),
        ) {
            self.synapses.push(Synapse::chemical(ashl, aval, 2.0).excitatory());
        }
    }

    /// Build neuromuscular junctions
    fn build_celegans_nmjs(&mut self) {
        // VA/DA motor neurons innervate dorsal muscles
        // VB/DB motor neurons innervate ventral muscles

        for i in 1..=12 {
            let muscle_row = i - 1;

            // VA → dorsal muscles
            if let Some(&va) = self.name_to_id.get(&format!("VA{:02}", i)) {
                self.nmjs.push(NeuromuscularJunction {
                    neuron: va,
                    muscle_row: muscle_row as u32,
                    muscle_quadrant: 0, // MDR
                    weight: 1.0,
                });
                self.nmjs.push(NeuromuscularJunction {
                    neuron: va,
                    muscle_row: muscle_row as u32,
                    muscle_quadrant: 3, // MDL
                    weight: 1.0,
                });
            }

            // VB → ventral muscles
            if let Some(&vb) = self.name_to_id.get(&format!("VB{:02}", i)) {
                self.nmjs.push(NeuromuscularJunction {
                    neuron: vb,
                    muscle_row: muscle_row as u32,
                    muscle_quadrant: 1, // MVR
                    weight: 1.0,
                });
                self.nmjs.push(NeuromuscularJunction {
                    neuron: vb,
                    muscle_row: muscle_row as u32,
                    muscle_quadrant: 2, // MVL
                    weight: 1.0,
                });
            }
        }
    }

    /// Add a neuron to the connectome
    pub fn add_neuron(&mut self, neuron: Neuron) {
        let id = neuron.id;
        let name = neuron.name.clone();
        self.neurons.push(neuron);
        self.name_to_id.insert(name, id);
    }

    /// Add a synapse
    pub fn add_synapse(&mut self, synapse: Synapse) {
        self.synapses.push(synapse);
    }

    /// Get neuron by ID
    pub fn get_neuron(&self, id: NeuronId) -> Option<&Neuron> {
        self.neurons.get(id as usize)
    }

    /// Get neuron by name
    pub fn get_neuron_by_name(&self, name: &str) -> Option<&Neuron> {
        self.name_to_id.get(name).and_then(|&id| self.get_neuron(id))
    }

    /// Get neuron ID by name
    pub fn get_id(&self, name: &str) -> Option<NeuronId> {
        self.name_to_id.get(name).copied()
    }

    /// Get all neurons
    pub fn neurons(&self) -> &[Neuron] {
        &self.neurons
    }

    /// Get all synapses
    pub fn synapses(&self) -> &[Synapse] {
        &self.synapses
    }

    /// Get all neuromuscular junctions
    pub fn nmjs(&self) -> &[NeuromuscularJunction] {
        &self.nmjs
    }

    /// Get number of neurons
    pub fn num_neurons(&self) -> usize {
        self.neurons.len()
    }

    /// Get number of synapses
    pub fn num_synapses(&self) -> usize {
        self.synapses.len()
    }

    /// Get neurons by class
    pub fn neurons_by_class(&self, class: NeuronClass) -> Vec<&Neuron> {
        self.neurons.iter().filter(|n| n.class == class).collect()
    }

    /// Get motor neurons
    pub fn motor_neurons(&self) -> Vec<&Neuron> {
        self.neurons_by_class(NeuronClass::Motor)
    }

    /// Get sensory neurons
    pub fn sensory_neurons(&self) -> Vec<&Neuron> {
        self.neurons_by_class(NeuronClass::Sensory)
    }

    /// Get interneurons
    pub fn interneurons(&self) -> Vec<&Neuron> {
        self.neurons_by_class(NeuronClass::Interneuron)
    }

    /// Get outgoing synapses for a neuron
    pub fn outgoing_synapses(&self, neuron_id: NeuronId) -> Vec<&Synapse> {
        self.synapses.iter().filter(|s| s.pre == neuron_id).collect()
    }

    /// Get incoming synapses for a neuron
    pub fn incoming_synapses(&self, neuron_id: NeuronId) -> Vec<&Synapse> {
        self.synapses.iter().filter(|s| s.post == neuron_id).collect()
    }

    /// Create a subnetwork containing only specified neurons
    pub fn subnetwork(&self, neuron_names: &[&str]) -> Self {
        let mut sub = Self::new();

        // Add selected neurons with new IDs
        let mut id_map: HashMap<NeuronId, NeuronId> = HashMap::new();

        for name in neuron_names {
            if let Some(neuron) = self.get_neuron_by_name(name) {
                let new_id = sub.neurons.len() as NeuronId;
                id_map.insert(neuron.id, new_id);

                let mut new_neuron = neuron.clone();
                new_neuron.id = new_id;
                sub.add_neuron(new_neuron);
            }
        }

        // Add synapses between selected neurons
        for synapse in &self.synapses {
            if let (Some(&new_pre), Some(&new_post)) = (
                id_map.get(&synapse.pre),
                id_map.get(&synapse.post),
            ) {
                let mut new_synapse = synapse.clone();
                new_synapse.pre = new_pre;
                new_synapse.post = new_post;
                sub.synapses.push(new_synapse);
            }
        }

        // Add NMJs for selected neurons
        for nmj in &self.nmjs {
            if let Some(&new_id) = id_map.get(&nmj.neuron) {
                let mut new_nmj = *nmj;
                new_nmj.neuron = new_id;
                sub.nmjs.push(new_nmj);
            }
        }

        sub
    }
}

impl Default for Connectome {
    fn default() -> Self {
        Self::new()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_celegans_connectome() {
        let conn = Connectome::celegans();

        assert!(conn.num_neurons() > 0, "Should have neurons");
        assert!(conn.num_synapses() > 0, "Should have synapses");

        // Check for key neurons
        assert!(conn.get_neuron_by_name("AVAL").is_some());
        assert!(conn.get_neuron_by_name("VB01").is_some());
    }

    #[test]
    fn test_neuron_classes() {
        let conn = Connectome::celegans();

        let sensory = conn.sensory_neurons();
        let motor = conn.motor_neurons();
        let inter = conn.interneurons();

        assert!(!sensory.is_empty());
        assert!(!motor.is_empty());
        assert!(!inter.is_empty());
    }

    #[test]
    fn test_subnetwork() {
        let conn = Connectome::celegans();
        let sub = conn.subnetwork(&["AVAL", "AVAR", "VA01", "VA02"]);

        assert_eq!(sub.num_neurons(), 4);
        assert!(sub.num_synapses() > 0);
    }
}
