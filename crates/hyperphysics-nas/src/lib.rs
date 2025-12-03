//! HyperPhysics NAS - Neural Architecture Search
//!
//! Evolutionary algorithms for discovering optimal neural network architectures.
//!
//! # Overview
//!
//! This crate provides tools for evolving neural network topologies,
//! inspired by NEAT (NeuroEvolution of Augmenting Topologies) and
//! modern differentiable architecture search methods.
//!
//! # Key Features
//!
//! - **Genome Encoding**: Efficient representation of network topology
//! - **Mutation Operators**: Add/remove neurons, connections, modify weights
//! - **Crossover**: Combine successful architectures
//! - **Speciation**: Protect topological innovations
//! - **Fitness Evaluation**: Multi-objective optimization support
//!
//! # Integration with Consciousness Metrics
//!
//! The NAS system can optimize for consciousness-related metrics:
//! - Integrated Information (Φ)
//! - Causal Density
//! - Metastability
//! - Information Transfer

mod genome;
mod mutation;
mod crossover;
mod speciation;
mod fitness;
mod population;
mod evolution;

pub use genome::{Genome, ConnectionGene, NodeGene, NodeType, InnovationTracker};
pub use mutation::{MutationParams, MutationOperator};
pub use crossover::{CrossoverParams, CrossoverOperator};
pub use speciation::{Species, SpeciationParams, SpeciesManager};
pub use fitness::{FitnessFunction, MultiFitness, FitnessMetric};
pub use population::{Population, Individual, PopulationParams};
pub use evolution::{EvolutionEngine, EvolutionParams, GenerationStats};

use thiserror::Error;

/// NAS errors
#[derive(Debug, Error)]
pub enum NasError {
    #[error("Invalid genome: {0}")]
    InvalidGenome(String),

    #[error("Population extinction: all individuals have zero fitness")]
    PopulationExtinction,

    #[error("Innovation tracker overflow")]
    InnovationOverflow,

    #[error("Incompatible genomes for crossover")]
    IncompatibleGenomes,

    #[error("Fitness evaluation failed: {0}")]
    FitnessError(String),
}

pub type Result<T> = std::result::Result<T, NasError>;

/// Configuration presets for different evolution scenarios
pub mod presets {
    use super::*;

    /// Small network evolution (quick exploration)
    pub fn small_network() -> EvolutionParams {
        EvolutionParams {
            population_size: 50,
            max_generations: 100,
            mutation: MutationParams::exploratory(),
            crossover: CrossoverParams::default(),
            speciation: SpeciationParams::default(),
            elitism: 2,
            ..Default::default()
        }
    }

    /// Large network evolution (thorough search)
    pub fn large_network() -> EvolutionParams {
        EvolutionParams {
            population_size: 200,
            max_generations: 500,
            mutation: MutationParams::default(),
            crossover: CrossoverParams::default(),
            speciation: SpeciationParams::protective(),
            elitism: 5,
            ..Default::default()
        }
    }

    /// Consciousness optimization
    pub fn consciousness() -> EvolutionParams {
        EvolutionParams {
            population_size: 100,
            max_generations: 1000,
            mutation: MutationParams::balanced(),
            crossover: CrossoverParams::default(),
            speciation: SpeciationParams::protective(),
            elitism: 3,
            multi_objective: true,
            ..Default::default()
        }
    }
}
