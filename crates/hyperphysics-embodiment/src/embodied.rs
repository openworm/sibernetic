//! Embodied Worm Simulation
//!
//! Main simulation coordinator that integrates neural network with body physics.

use crate::actuator::Actuator;
use crate::coupling::{CouplingConfig, CouplingMode, SegmentMapping};
use crate::proprioception::{
    compute_proprioceptive_state, ProprioceptiveOutput, ProprioceptiveState, ProprioceptiveSystem,
};
use crate::time_sync::{compute_schedule, IntegrationCycle, IntegrationPhase, TimeSync};
use crate::{EmbodimentError, Result};

use hyperphysics_connectome::{Connectome, ModelLevel, SpikingNetwork};
use hyperphysics_sph::{MuscleActivation, SphWorld};

/// Full simulation state snapshot
#[derive(Debug, Clone)]
pub struct SimulationState {
    /// Current time (ms)
    pub time: f64,

    /// Neural network state
    pub neural_voltages: Vec<f32>,

    /// Muscle activations
    pub muscle_activations: [f32; 96],

    /// Body center of mass
    pub center_of_mass: [f64; 3],

    /// Body velocity (center of mass)
    pub velocity: [f64; 3],

    /// Total kinetic energy
    pub kinetic_energy: f64,

    /// Proprioceptive state per segment
    pub proprioception: Vec<ProprioceptiveState>,

    /// Spike count in last window
    pub spike_count: usize,
}

impl Default for SimulationState {
    fn default() -> Self {
        Self {
            time: 0.0,
            neural_voltages: Vec::new(),
            muscle_activations: [0.0; 96],
            center_of_mass: [0.0; 3],
            velocity: [0.0; 3],
            kinetic_energy: 0.0,
            proprioception: Vec::new(),
            spike_count: 0,
        }
    }
}

/// Embodied C. elegans simulation
pub struct EmbodiedWorm {
    /// Neural network
    neural: SpikingNetwork,

    /// Body physics
    body: SphWorld,

    /// Coupling configuration
    coupling: CouplingConfig,

    /// Time synchronizer
    time_sync: TimeSync,

    /// Muscle actuator system
    actuator: Actuator,

    /// Proprioceptive system
    proprioception: ProprioceptiveSystem,

    /// Segment mapping
    segment_mapping: SegmentMapping,

    /// Previous segment centers (for velocity computation)
    prev_segment_centers: Vec<[f32; 3]>,

    /// Current simulation time (ms)
    time: f64,

    /// Step counter
    step_count: u64,

    /// Last muscle output from neural network
    last_muscle_output: [f32; 96],

    /// Cached proprioceptive outputs
    proprioceptive_outputs: Vec<ProprioceptiveOutput>,
}

impl EmbodiedWorm {
    /// Create new embodied worm simulation
    pub fn new(level: ModelLevel, coupling: CouplingConfig) -> Self {
        let neural = SpikingNetwork::celegans(level);
        let body = SphWorld::worm();
        let physics_dt = body.config().physics.time_step;
        let time_sync = TimeSync::for_model_level(level, physics_dt * 1000.0); // Convert to ms

        Self {
            neural,
            body,
            coupling,
            time_sync,
            actuator: Actuator::celegans(),
            proprioception: ProprioceptiveSystem::celegans(),
            segment_mapping: SegmentMapping::default(),
            prev_segment_centers: vec![[0.0; 3]; 24],
            time: 0.0,
            step_count: 0,
            last_muscle_output: [0.0; 96],
            proprioceptive_outputs: Vec::new(),
        }
    }

    /// Create with default settings
    pub fn default_level(level: ModelLevel) -> Self {
        Self::new(level, CouplingConfig::default())
    }

    /// Create fast simulation (simplified model)
    pub fn fast() -> Self {
        Self::new(ModelLevel::A, CouplingConfig::open_loop())
    }

    /// Create high-fidelity simulation
    pub fn high_fidelity() -> Self {
        Self::new(ModelLevel::C, CouplingConfig::high_fidelity())
    }

    /// Initialize body with worm geometry
    pub fn initialize_body(&mut self) -> Result<()> {
        // This would load a worm body model
        // For now, create a simple cylindrical approximation
        self.create_cylindrical_body(24, 8, 0.05)?;
        Ok(())
    }

    /// Create cylindrical worm body
    fn create_cylindrical_body(
        &mut self,
        num_segments: usize,
        particles_per_ring: usize,
        spacing: f32,
    ) -> Result<()> {
        use hyperphysics_sph::ParticleType;
        use std::f32::consts::PI;

        let radius = 0.01; // 10 μm radius (scaled)
        let length = spacing * num_segments as f32;

        let mut positions = Vec::new();

        // Create particle rings along body
        for seg in 0..num_segments {
            let x = seg as f32 * spacing;

            for ring_idx in 0..particles_per_ring {
                let angle = 2.0 * PI * ring_idx as f32 / particles_per_ring as f32;
                let y = radius * angle.cos();
                let z = radius * angle.sin();

                let idx = self
                    .body
                    .add_particle([x, y, z], [0.0, 0.0, 0.0], ParticleType::Elastic)?;

                positions.push([x, y, z]);

                // Connect to neighbors
                if ring_idx > 0 {
                    self.body.connect(idx - 1, idx, 100.0);
                }
                if seg > 0 {
                    let prev_ring_idx = (seg - 1) * particles_per_ring + ring_idx;
                    self.body.connect(prev_ring_idx, idx, 100.0);

                    // Muscle connections (assign to appropriate muscle)
                    let muscle_id = self.particle_to_muscle(seg, ring_idx, particles_per_ring);
                    self.body
                        .connect_muscle(prev_ring_idx, idx, 200.0, muscle_id);
                }
            }

            // Close ring
            let ring_start = seg * particles_per_ring;
            let ring_end = ring_start + particles_per_ring - 1;
            self.body.connect(ring_end, ring_start, 100.0);
        }

        // Build segment mapping
        self.segment_mapping = SegmentMapping::from_particles(&positions, num_segments);
        self.prev_segment_centers = self.segment_mapping.segment_centers.clone();

        Ok(())
    }

    /// Map particle to muscle index
    fn particle_to_muscle(&self, segment: usize, ring_idx: usize, ring_size: usize) -> i32 {
        // Divide ring into 4 quadrants
        let quadrant = (ring_idx * 4 / ring_size) % 4;
        (segment * 4 + quadrant) as i32
    }

    /// Perform one integration step
    pub fn step(&mut self) {
        let schedule = compute_schedule(&self.time_sync);
        let mut cycle = IntegrationCycle::new(schedule);

        while !cycle.is_complete() {
            match cycle.phase {
                IntegrationPhase::Neural => {
                    // Run neural simulation step
                    self.neural.step(self.time_sync.neural_dt);
                    self.time_sync.advance_neural();
                }

                IntegrationPhase::TransferToBody => {
                    // Get muscle activation from neural network
                    if self.coupling.mode.neural_drives_body() {
                        self.last_muscle_output = self.neural.get_muscle_output();
                        self.actuator
                            .update(&self.last_muscle_output, self.time_sync.physics_dt);
                        let activation = self.actuator.get_activation();
                        self.body.set_muscle_activations(&activation);
                    }
                }

                IntegrationPhase::Physics => {
                    // Run physics simulation step
                    self.body.step();
                    self.time_sync.advance_physics();
                }

                IntegrationPhase::TransferToNeural => {
                    // Compute proprioceptive feedback
                    if self.coupling.mode.body_drives_neural() && self.coupling.proprioception_enabled
                    {
                        self.update_proprioception();
                    }
                }

                IntegrationPhase::Complete => {}
            }

            cycle.advance();
        }

        self.time = self.time_sync.physics_time();
        self.step_count += 1;
    }

    /// Update proprioceptive feedback
    fn update_proprioception(&mut self) {
        // Update segment centers from body
        self.update_segment_centers();

        // Compute proprioceptive state
        let states = compute_proprioceptive_state(
            &self.segment_mapping.segment_centers,
            &self.prev_segment_centers,
            self.time_sync.physics_dt,
        );

        // Process through proprioceptive system
        self.proprioceptive_outputs = self
            .proprioception
            .process_all(&states, self.time_sync.physics_dt);

        // Apply proprioceptive input to neural network
        for (seg_idx, output) in self.proprioceptive_outputs.iter().enumerate() {
            // Apply stretch receptor outputs to sensory neurons
            // This is a simplified mapping - would need actual neuron IDs
            let base_current = output.magnitude() * self.coupling.proprioceptive_gain;

            // Apply to segment-specific proprioceptive neurons
            // In reality, this would target specific touch receptor neurons
            let neuron_id = (250 + seg_idx * 2) as u16;
            self.neural.set_input(neuron_id, base_current);
        }

        // Store current centers as previous for next step
        self.prev_segment_centers = self.segment_mapping.segment_centers.clone();
    }

    /// Update segment centers from particle positions
    fn update_segment_centers(&mut self) {
        for (seg_idx, particles) in self.segment_mapping.segment_particles.iter().enumerate() {
            if particles.is_empty() {
                continue;
            }

            let mut center = [0.0_f32; 3];
            for &p_idx in particles {
                let pos = self.body.get_position(p_idx);
                center[0] += pos[0];
                center[1] += pos[1];
                center[2] += pos[2];
            }

            let n = particles.len() as f32;
            self.segment_mapping.segment_centers[seg_idx] =
                [center[0] / n, center[1] / n, center[2] / n];
        }
    }

    /// Run simulation for given duration (ms)
    pub fn run(&mut self, duration_ms: f64) {
        let steps = (duration_ms / self.time_sync.physics_dt as f64).ceil() as u64;
        for _ in 0..steps {
            self.step();
        }
    }

    /// Apply external stimulus to sensory neurons
    pub fn stimulate(&mut self, neuron_name: &str, current: f32) {
        self.neural.set_input_by_name(neuron_name, current);
    }

    /// Apply touch stimulus (activates touch neurons based on position)
    pub fn touch(&mut self, position: [f32; 3], strength: f32) {
        // Find nearest segment
        let mut min_dist = f32::INFINITY;
        let mut nearest_seg = 0;

        for (seg_idx, center) in self.segment_mapping.segment_centers.iter().enumerate() {
            let dx = position[0] - center[0];
            let dy = position[1] - center[1];
            let dz = position[2] - center[2];
            let dist = (dx * dx + dy * dy + dz * dz).sqrt();

            if dist < min_dist {
                min_dist = dist;
                nearest_seg = seg_idx;
            }
        }

        // Activate touch neurons for that segment
        // In C. elegans: ALM, AVM, PLM for gentle touch; ASH, PVD for harsh
        let touch_neurons = match nearest_seg {
            0..=5 => vec!["ALML", "ALMR", "AVM"], // Anterior
            6..=17 => vec!["ALML", "ALMR", "PLML", "PLMR"], // Middle
            _ => vec!["PLML", "PLMR", "PVM"],     // Posterior
        };

        for name in touch_neurons {
            self.neural.set_input_by_name(name, strength);
        }
    }

    /// Get current simulation state
    pub fn state(&self) -> SimulationState {
        let com = self.body.center_of_mass();

        SimulationState {
            time: self.time,
            neural_voltages: self.neural.get_voltages(),
            muscle_activations: self.last_muscle_output,
            center_of_mass: com,
            velocity: [0.0; 3], // Would need to track previous COM
            kinetic_energy: self.body.kinetic_energy(),
            proprioception: vec![ProprioceptiveState::default(); 24], // Simplified
            spike_count: self.neural.spike_history().len(),
        }
    }

    /// Reset simulation
    pub fn reset(&mut self) {
        self.neural.reset();
        self.body.reset();
        self.actuator.reset();
        self.proprioception.reset();
        self.time_sync.reset();
        self.time = 0.0;
        self.step_count = 0;
        self.last_muscle_output = [0.0; 96];
    }

    // ========== Getters ==========

    /// Get current time (ms)
    pub fn time(&self) -> f64 {
        self.time
    }

    /// Get step count
    pub fn step_count(&self) -> u64 {
        self.step_count
    }

    /// Get neural network reference
    pub fn neural(&self) -> &SpikingNetwork {
        &self.neural
    }

    /// Get mutable neural network reference
    pub fn neural_mut(&mut self) -> &mut SpikingNetwork {
        &mut self.neural
    }

    /// Get body physics reference
    pub fn body(&self) -> &SphWorld {
        &self.body
    }

    /// Get mutable body physics reference
    pub fn body_mut(&mut self) -> &mut SphWorld {
        &mut self.body
    }

    /// Get muscle activations
    pub fn muscle_activations(&self) -> &[f32; 96] {
        &self.last_muscle_output
    }

    /// Get coupling configuration
    pub fn coupling(&self) -> &CouplingConfig {
        &self.coupling
    }

    /// Get energy expenditure
    pub fn energy_expenditure(&self) -> f32 {
        self.actuator.energy_expenditure()
    }
}

/// Builder for EmbodiedWorm
pub struct EmbodiedWormBuilder {
    level: ModelLevel,
    coupling: CouplingConfig,
    initialize_body: bool,
}

impl Default for EmbodiedWormBuilder {
    fn default() -> Self {
        Self {
            level: ModelLevel::B,
            coupling: CouplingConfig::default(),
            initialize_body: true,
        }
    }
}

impl EmbodiedWormBuilder {
    /// Create new builder
    pub fn new() -> Self {
        Self::default()
    }

    /// Set model level
    pub fn level(mut self, level: ModelLevel) -> Self {
        self.level = level;
        self
    }

    /// Set coupling configuration
    pub fn coupling(mut self, coupling: CouplingConfig) -> Self {
        self.coupling = coupling;
        self
    }

    /// Disable automatic body initialization
    pub fn no_body(mut self) -> Self {
        self.initialize_body = false;
        self
    }

    /// Build the simulation
    pub fn build(self) -> Result<EmbodiedWorm> {
        let mut worm = EmbodiedWorm::new(self.level, self.coupling);

        if self.initialize_body {
            worm.initialize_body()?;
        }

        Ok(worm)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_embodied_creation() {
        let worm = EmbodiedWorm::fast();
        assert_eq!(worm.time(), 0.0);
    }

    #[test]
    fn test_builder() {
        let worm = EmbodiedWormBuilder::new()
            .level(ModelLevel::A)
            .coupling(CouplingConfig::open_loop())
            .no_body()
            .build()
            .unwrap();

        assert_eq!(worm.coupling().proprioception_enabled, false);
    }

    #[test]
    fn test_simulation_step() {
        let mut worm = EmbodiedWorm::fast();

        // Step without body (just neural)
        worm.step();

        assert!(worm.time() > 0.0);
        assert_eq!(worm.step_count(), 1);
    }
}
