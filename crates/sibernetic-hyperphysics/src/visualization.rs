//! Visualization Support
//!
//! Data structures and utilities for visualization.

#[cfg(feature = "serde")]
use serde::{Deserialize, Serialize};

/// Visualization state
#[derive(Debug, Clone)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct VisualizationState {
    /// Current simulation time
    pub time: f64,

    /// Particle positions (flat array: x,y,z,w for each particle)
    pub positions: Vec<f32>,

    /// Particle velocities (flat array)
    pub velocities: Vec<f32>,

    /// Particle types
    pub types: Vec<u8>,

    /// Muscle activations (96 values)
    pub muscle_activations: [f32; 96],

    /// Elastic connections (pairs of particle indices)
    pub connections: Vec<(u32, u32)>,

    /// Membrane triangles (triplets of particle indices)
    pub membranes: Vec<(u32, u32, u32)>,

    /// Neural spike events (neuron IDs)
    pub spikes: Vec<u16>,

    /// Center of mass
    pub center_of_mass: [f64; 3],

    /// Total kinetic energy
    pub kinetic_energy: f64,
}

impl Default for VisualizationState {
    fn default() -> Self {
        Self {
            time: 0.0,
            positions: Vec::new(),
            velocities: Vec::new(),
            types: Vec::new(),
            muscle_activations: [0.0; 96],
            connections: Vec::new(),
            membranes: Vec::new(),
            spikes: Vec::new(),
            center_of_mass: [0.0; 3],
            kinetic_energy: 0.0,
        }
    }
}

impl VisualizationState {
    /// Create from simulation
    pub fn from_simulation(sim: &crate::WormSimulation) -> Self {
        let body = sim.embodiment().body();

        Self {
            time: sim.time(),
            positions: body.positions().to_vec(),
            velocities: body.velocities().to_vec(),
            types: body.particles().types.iter().map(|t| *t as u8).collect(),
            muscle_activations: *sim.get_muscle_activations(),
            connections: Vec::new(), // Would extract from elastic network
            membranes: Vec::new(),   // Would extract from membrane
            spikes: sim.get_spikes(),
            center_of_mass: sim.get_center_of_mass(),
            kinetic_energy: sim.get_kinetic_energy(),
        }
    }

    /// Get number of particles
    pub fn num_particles(&self) -> usize {
        self.positions.len() / 4
    }

    /// Get position of particle
    pub fn get_position(&self, idx: usize) -> Option<[f32; 3]> {
        let base = idx * 4;
        if base + 2 < self.positions.len() {
            Some([
                self.positions[base],
                self.positions[base + 1],
                self.positions[base + 2],
            ])
        } else {
            None
        }
    }
}

/// Render data for a single frame
#[derive(Debug, Clone)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct RenderData {
    /// Frame number
    pub frame: u64,

    /// Time
    pub time: f64,

    /// Vertex positions
    pub vertices: Vec<f32>,

    /// Vertex colors (RGBA)
    pub colors: Vec<f32>,

    /// Indices for triangles
    pub indices: Vec<u32>,

    /// Camera position
    pub camera_pos: [f32; 3],

    /// Camera target
    pub camera_target: [f32; 3],
}

impl Default for RenderData {
    fn default() -> Self {
        Self {
            frame: 0,
            time: 0.0,
            vertices: Vec::new(),
            colors: Vec::new(),
            indices: Vec::new(),
            camera_pos: [0.0, 5.0, 10.0],
            camera_target: [0.0, 0.0, 0.0],
        }
    }
}

/// Convert visualization state to render data
pub fn to_render_data(state: &VisualizationState, frame: u64) -> RenderData {
    let num_particles = state.num_particles();

    let mut vertices = Vec::with_capacity(num_particles * 3);
    let mut colors = Vec::with_capacity(num_particles * 4);

    for i in 0..num_particles {
        if let Some(pos) = state.get_position(i) {
            vertices.extend_from_slice(&pos);

            // Color based on particle type
            let color = match state.types.get(i).copied().unwrap_or(0) {
                0 => [0.3, 0.5, 0.8, 1.0], // Liquid - blue
                1 => [0.2, 0.8, 0.3, 1.0], // Elastic - green
                2 => [0.6, 0.6, 0.6, 1.0], // Boundary - gray
                _ => [1.0, 1.0, 1.0, 1.0], // Default - white
            };
            colors.extend_from_slice(&color);
        }
    }

    RenderData {
        frame,
        time: state.time,
        vertices,
        colors,
        indices: Vec::new(), // Point cloud, no indices needed
        camera_pos: [
            state.center_of_mass[0] as f32,
            state.center_of_mass[1] as f32 + 2.0,
            state.center_of_mass[2] as f32 + 5.0,
        ],
        camera_target: [
            state.center_of_mass[0] as f32,
            state.center_of_mass[1] as f32,
            state.center_of_mass[2] as f32,
        ],
    }
}

/// Muscle visualization
#[derive(Debug, Clone)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct MuscleVisualization {
    /// 24 rows x 4 quadrants
    pub grid: [[f32; 4]; 24],

    /// Quadrant names
    pub quadrant_names: [&'static str; 4],

    /// Color map (activation -> RGB)
    pub color_map: ColorMap,
}

impl Default for MuscleVisualization {
    fn default() -> Self {
        Self {
            grid: [[0.0; 4]; 24],
            quadrant_names: ["MDR", "MVR", "MVL", "MDL"],
            color_map: ColorMap::HeatMap,
        }
    }
}

impl MuscleVisualization {
    /// Update from muscle activations
    pub fn update(&mut self, activations: &[f32; 96]) {
        for row in 0..24 {
            for quad in 0..4 {
                self.grid[row][quad] = activations[row * 4 + quad];
            }
        }
    }

    /// Get color for activation value
    pub fn activation_color(&self, activation: f32) -> [f32; 3] {
        self.color_map.map(activation)
    }
}

/// Color mapping for visualization
#[derive(Debug, Clone, Copy)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub enum ColorMap {
    /// Blue -> Red
    HeatMap,
    /// Black -> White
    Grayscale,
    /// Green -> Yellow -> Red
    Traffic,
    /// Blue -> Green -> Yellow -> Red
    Viridis,
}

impl ColorMap {
    /// Map value (0-1) to RGB color
    pub fn map(&self, value: f32) -> [f32; 3] {
        let v = value.clamp(0.0, 1.0);

        match self {
            ColorMap::HeatMap => {
                // Blue -> Cyan -> Green -> Yellow -> Red
                if v < 0.25 {
                    [0.0, v * 4.0, 1.0]
                } else if v < 0.5 {
                    [0.0, 1.0, 1.0 - (v - 0.25) * 4.0]
                } else if v < 0.75 {
                    [(v - 0.5) * 4.0, 1.0, 0.0]
                } else {
                    [1.0, 1.0 - (v - 0.75) * 4.0, 0.0]
                }
            }
            ColorMap::Grayscale => [v, v, v],
            ColorMap::Traffic => {
                if v < 0.5 {
                    [v * 2.0, 1.0, 0.0]
                } else {
                    [1.0, 1.0 - (v - 0.5) * 2.0, 0.0]
                }
            }
            ColorMap::Viridis => {
                // Approximation of viridis colormap
                let r = 0.267 + v * (0.993 - 0.267);
                let g = 0.004 + v * 0.5 * (1.0 - v);
                let b = 0.329 * (1.0 - v);
                [r, g.max(0.0), b]
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_visualization_state() {
        let state = VisualizationState::default();
        assert_eq!(state.num_particles(), 0);
    }

    #[test]
    fn test_color_map() {
        let heat = ColorMap::HeatMap;

        let cold = heat.map(0.0);
        let hot = heat.map(1.0);

        // Cold should be blue-ish
        assert!(cold[2] > cold[0]);

        // Hot should be red-ish
        assert!(hot[0] > hot[2]);
    }

    #[test]
    fn test_muscle_visualization() {
        let mut viz = MuscleVisualization::default();
        let mut activations = [0.0_f32; 96];
        activations[0] = 0.5;
        activations[95] = 1.0;

        viz.update(&activations);

        assert!((viz.grid[0][0] - 0.5).abs() < 1e-6);
        assert!((viz.grid[23][3] - 1.0).abs() < 1e-6);
    }
}
