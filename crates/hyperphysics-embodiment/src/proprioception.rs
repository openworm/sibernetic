//! Proprioceptive Feedback
//!
//! Models sensory feedback from body mechanics to neural network.
//! C. elegans uses stretch receptors in body wall muscles and
//! mechanosensory neurons for proprioception.

#[cfg(feature = "serde")]
use serde::{Deserialize, Serialize};

/// Proprioceptive state for one body segment
#[derive(Debug, Clone, Copy, Default)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct ProprioceptiveState {
    /// Dorsal-ventral bend angle (radians)
    pub bend_angle: f32,

    /// Left-right bend angle (radians)
    pub lateral_angle: f32,

    /// Rate of bend change (rad/s)
    pub bend_velocity: f32,

    /// Muscle stretch (normalized, -1 to 1)
    pub stretch: [f32; 4], // DR, VR, VL, DL

    /// Local curvature
    pub curvature: f32,

    /// Segment velocity (body frame)
    pub velocity: [f32; 3],
}

impl ProprioceptiveState {
    /// Compute total bend magnitude
    pub fn bend_magnitude(&self) -> f32 {
        (self.bend_angle * self.bend_angle + self.lateral_angle * self.lateral_angle).sqrt()
    }

    /// Is this segment bending dorsally?
    pub fn is_dorsal_bend(&self) -> bool {
        self.bend_angle > 0.0
    }

    /// Is this segment bending ventrally?
    pub fn is_ventral_bend(&self) -> bool {
        self.bend_angle < 0.0
    }
}

/// Stretch receptor model
#[derive(Debug, Clone)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct StretchReceptor {
    /// Sensitivity to stretch
    pub sensitivity: f32,

    /// Threshold for activation
    pub threshold: f32,

    /// Saturation level
    pub saturation: f32,

    /// Adaptation time constant (ms)
    pub tau_adapt: f32,

    /// Current adaptation state
    adaptation: f32,
}

impl Default for StretchReceptor {
    fn default() -> Self {
        Self {
            sensitivity: 10.0,
            threshold: 0.05,
            saturation: 50.0,
            tau_adapt: 100.0,
            adaptation: 0.0,
        }
    }
}

impl StretchReceptor {
    /// Compute receptor output for given stretch
    pub fn respond(&mut self, stretch: f32, dt: f32) -> f32 {
        // Rectified stretch
        let s = (stretch.abs() - self.threshold).max(0.0);

        // Apply adaptation
        let adapted = s - self.adaptation;
        self.adaptation += (s - self.adaptation) * dt / self.tau_adapt;

        // Apply nonlinearity (saturating)
        let response = self.sensitivity * adapted;
        response.clamp(-self.saturation, self.saturation)
    }

    /// Reset adaptation
    pub fn reset(&mut self) {
        self.adaptation = 0.0;
    }
}

/// Proprioceptor for one body segment
/// Converts mechanical state to neural input
#[derive(Debug, Clone)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct Proprioceptor {
    /// Stretch receptors for each quadrant
    stretch_receptors: [StretchReceptor; 4],

    /// Curvature sensitivity
    curvature_sensitivity: f32,

    /// Velocity sensitivity
    velocity_sensitivity: f32,

    /// Target sensory neuron IDs
    pub target_neurons: Vec<u16>,

    /// Previous bend for velocity computation
    prev_bend: f32,
}

impl Default for Proprioceptor {
    fn default() -> Self {
        Self {
            stretch_receptors: [
                StretchReceptor::default(),
                StretchReceptor::default(),
                StretchReceptor::default(),
                StretchReceptor::default(),
            ],
            curvature_sensitivity: 5.0,
            velocity_sensitivity: 2.0,
            target_neurons: Vec::new(),
            prev_bend: 0.0,
        }
    }
}

impl Proprioceptor {
    /// Create proprioceptor targeting specific neurons
    pub fn with_targets(neurons: Vec<u16>) -> Self {
        Self {
            target_neurons: neurons,
            ..Default::default()
        }
    }

    /// Process mechanical state and return neural inputs
    pub fn process(&mut self, state: &ProprioceptiveState, dt: f32) -> ProprioceptiveOutput {
        // Stretch receptor responses
        let stretch_responses: [f32; 4] = [
            self.stretch_receptors[0].respond(state.stretch[0], dt),
            self.stretch_receptors[1].respond(state.stretch[1], dt),
            self.stretch_receptors[2].respond(state.stretch[2], dt),
            self.stretch_receptors[3].respond(state.stretch[3], dt),
        ];

        // Curvature response
        let curvature_response = state.curvature * self.curvature_sensitivity;

        // Bend velocity (differentiate)
        let bend_velocity = (state.bend_angle - self.prev_bend) / dt;
        self.prev_bend = state.bend_angle;
        let velocity_response = bend_velocity * self.velocity_sensitivity;

        ProprioceptiveOutput {
            stretch_responses,
            curvature_response,
            velocity_response,
            bend_angle: state.bend_angle,
        }
    }

    /// Reset state
    pub fn reset(&mut self) {
        for sr in &mut self.stretch_receptors {
            sr.reset();
        }
        self.prev_bend = 0.0;
    }
}

/// Output from proprioceptor processing
#[derive(Debug, Clone, Copy, Default)]
pub struct ProprioceptiveOutput {
    /// Stretch receptor outputs (one per quadrant)
    pub stretch_responses: [f32; 4],

    /// Curvature response
    pub curvature_response: f32,

    /// Velocity response
    pub velocity_response: f32,

    /// Current bend angle
    pub bend_angle: f32,
}

impl ProprioceptiveOutput {
    /// Total proprioceptive signal magnitude
    pub fn magnitude(&self) -> f32 {
        let stretch_sum: f32 = self.stretch_responses.iter().map(|x| x * x).sum();
        (stretch_sum
            + self.curvature_response * self.curvature_response
            + self.velocity_response * self.velocity_response)
            .sqrt()
    }
}

/// Full body proprioceptive system
#[derive(Debug, Clone)]
pub struct ProprioceptiveSystem {
    /// Proprioceptors for each segment
    segments: Vec<Proprioceptor>,

    /// Enable delay line for proprioceptive feedback
    delay_enabled: bool,

    /// Delay in steps
    delay_steps: usize,

    /// History buffer for delay
    history: Vec<Vec<ProprioceptiveState>>,
}

impl ProprioceptiveSystem {
    /// Create system with specified number of segments
    pub fn new(num_segments: usize) -> Self {
        Self {
            segments: vec![Proprioceptor::default(); num_segments],
            delay_enabled: false,
            delay_steps: 0,
            history: Vec::new(),
        }
    }

    /// Create C. elegans proprioceptive system (24 segments)
    pub fn celegans() -> Self {
        let mut system = Self::new(24);

        // Configure segment proprioceptors with appropriate target neurons
        // In C. elegans, PLM, PVD, and DVA neurons are involved in proprioception
        // This is a simplified mapping
        for (i, prop) in system.segments.iter_mut().enumerate() {
            // Map to appropriate sensory neurons based on segment
            let segment_neuron_offset = (i as u16) * 2;
            prop.target_neurons = vec![
                250 + segment_neuron_offset,     // Hypothetical proprioceptive neurons
                250 + segment_neuron_offset + 1,
            ];
        }

        system
    }

    /// Enable delay with given parameters
    pub fn set_delay(&mut self, delay_ms: f32, dt: f32) {
        self.delay_steps = (delay_ms / dt) as usize;
        self.delay_enabled = self.delay_steps > 0;
        self.history = vec![vec![ProprioceptiveState::default(); self.delay_steps]; self.segments.len()];
    }

    /// Process all segments
    pub fn process_all(
        &mut self,
        states: &[ProprioceptiveState],
        dt: f32,
    ) -> Vec<ProprioceptiveOutput> {
        states
            .iter()
            .enumerate()
            .map(|(i, state)| {
                let effective_state = if self.delay_enabled && i < self.history.len() {
                    // Use delayed state
                    let delayed = self.history[i].first().copied().unwrap_or_default();

                    // Update history
                    if !self.history[i].is_empty() {
                        self.history[i].remove(0);
                        self.history[i].push(*state);
                    }

                    delayed
                } else {
                    *state
                };

                if i < self.segments.len() {
                    self.segments[i].process(&effective_state, dt)
                } else {
                    ProprioceptiveOutput::default()
                }
            })
            .collect()
    }

    /// Get number of segments
    pub fn num_segments(&self) -> usize {
        self.segments.len()
    }

    /// Reset all proprioceptors
    pub fn reset(&mut self) {
        for seg in &mut self.segments {
            seg.reset();
        }
        for hist in &mut self.history {
            hist.fill(ProprioceptiveState::default());
        }
    }
}

/// Compute proprioceptive state from particle positions
pub fn compute_proprioceptive_state(
    segment_centers: &[[f32; 3]],
    prev_centers: &[[f32; 3]],
    dt: f32,
) -> Vec<ProprioceptiveState> {
    let n = segment_centers.len();
    let mut states = vec![ProprioceptiveState::default(); n];

    if n < 3 {
        return states;
    }

    for i in 1..n - 1 {
        let prev = segment_centers[i - 1];
        let curr = segment_centers[i];
        let next = segment_centers[i + 1];

        // Compute local tangent vectors
        let t1 = [
            curr[0] - prev[0],
            curr[1] - prev[1],
            curr[2] - prev[2],
        ];
        let t2 = [
            next[0] - curr[0],
            next[1] - curr[1],
            next[2] - curr[2],
        ];

        // Normalize
        let len1 = (t1[0] * t1[0] + t1[1] * t1[1] + t1[2] * t1[2]).sqrt();
        let len2 = (t2[0] * t2[0] + t2[1] * t2[1] + t2[2] * t2[2]).sqrt();

        if len1 > 1e-6 && len2 > 1e-6 {
            let t1_norm = [t1[0] / len1, t1[1] / len1, t1[2] / len1];
            let t2_norm = [t2[0] / len2, t2[1] / len2, t2[2] / len2];

            // Bend angle (from dot product)
            let dot = t1_norm[0] * t2_norm[0] + t1_norm[1] * t2_norm[1] + t1_norm[2] * t2_norm[2];
            let angle = dot.clamp(-1.0, 1.0).acos();

            // Determine direction (dorsal/ventral from y component)
            let cross_y = t1_norm[2] * t2_norm[0] - t1_norm[0] * t2_norm[2];
            states[i].bend_angle = if cross_y > 0.0 { angle } else { -angle };

            // Lateral angle from z component
            let cross_z = t1_norm[0] * t2_norm[1] - t1_norm[1] * t2_norm[0];
            states[i].lateral_angle = cross_z.atan2(dot);

            // Curvature
            states[i].curvature = 2.0 * angle / (len1 + len2);
        }

        // Velocity from position change
        if i < prev_centers.len() {
            states[i].velocity = [
                (curr[0] - prev_centers[i][0]) / dt,
                (curr[1] - prev_centers[i][1]) / dt,
                (curr[2] - prev_centers[i][2]) / dt,
            ];
        }

        // Stretch estimation from bend angle
        // Positive bend stretches dorsal muscles, compresses ventral
        let stretch = states[i].bend_angle * 0.5;
        states[i].stretch = [stretch, -stretch, -stretch, stretch];
    }

    states
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_stretch_receptor() {
        let mut receptor = StretchReceptor::default();

        // No response below threshold
        let r1 = receptor.respond(0.01, 0.1);
        assert!(r1.abs() < 0.01);

        // Response above threshold
        receptor.reset();
        let r2 = receptor.respond(0.2, 0.1);
        assert!(r2 > 0.0);
    }

    #[test]
    fn test_proprioceptive_state() {
        let state = ProprioceptiveState {
            bend_angle: 0.5,
            lateral_angle: 0.2,
            ..Default::default()
        };

        assert!(state.is_dorsal_bend());
        assert!(!state.is_ventral_bend());
        assert!(state.bend_magnitude() > 0.5);
    }

    #[test]
    fn test_compute_state() {
        // Straight line
        let centers: Vec<[f32; 3]> = (0..10)
            .map(|i| [i as f32 * 0.1, 0.0, 0.0])
            .collect();
        let prev = centers.clone();

        let states = compute_proprioceptive_state(&centers, &prev, 0.001);

        // Should have near-zero bend for straight line
        for state in &states[1..states.len() - 1] {
            assert!(state.bend_angle.abs() < 0.1);
        }
    }
}
