//! Evolution Engine
//!
//! High-level evolution controller for neural architecture search.

#[cfg(feature = "serde")]
use serde::{Deserialize, Serialize};

use rand::SeedableRng;
use rand_distr::rand::rngs::StdRng;

use crate::crossover::CrossoverParams;
use crate::fitness::FitnessFunction;
use crate::genome::Genome;
use crate::mutation::MutationParams;
use crate::population::{Individual, Population, PopulationParams};
use crate::speciation::SpeciationParams;
use crate::Result;

/// Evolution parameters
#[derive(Debug, Clone)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct EvolutionParams {
    /// Population size
    pub population_size: usize,

    /// Maximum generations
    pub max_generations: usize,

    /// Number of inputs
    pub num_inputs: usize,

    /// Number of outputs
    pub num_outputs: usize,

    /// Elitism count (best individuals preserved)
    pub elitism: usize,

    /// Target fitness (stop when reached)
    pub target_fitness: Option<f32>,

    /// Mutation parameters
    pub mutation: MutationParams,

    /// Crossover parameters
    pub crossover: CrossoverParams,

    /// Speciation parameters
    pub speciation: SpeciationParams,

    /// Enable multi-objective optimization
    pub multi_objective: bool,

    /// Random seed (None for random)
    pub seed: Option<u64>,

    /// Stagnation threshold (generations without improvement)
    pub stagnation_limit: usize,
}

impl Default for EvolutionParams {
    fn default() -> Self {
        Self {
            population_size: 100,
            max_generations: 500,
            num_inputs: 3,
            num_outputs: 1,
            elitism: 2,
            target_fitness: None,
            mutation: MutationParams::default(),
            crossover: CrossoverParams::default(),
            speciation: SpeciationParams::default(),
            multi_objective: false,
            seed: None,
            stagnation_limit: 50,
        }
    }
}

/// Generation statistics
#[derive(Debug, Clone, Default)]
pub struct GenerationStats {
    /// Generation number
    pub generation: usize,

    /// Best fitness
    pub best_fitness: f32,

    /// Average fitness
    pub avg_fitness: f32,

    /// Worst fitness
    pub worst_fitness: f32,

    /// Number of species
    pub num_species: usize,

    /// Average genome complexity (nodes + connections)
    pub avg_complexity: f32,

    /// Number of unique innovations
    pub num_innovations: u32,
}

/// Callback for evolution events
pub trait EvolutionCallback: Send + Sync {
    /// Called after each generation
    fn on_generation(&mut self, stats: &GenerationStats, population: &Population);

    /// Called when target is reached
    fn on_target_reached(&mut self, _genome: &Genome, _generation: usize) {}

    /// Called on stagnation
    fn on_stagnation(&mut self, _generation: usize) {}

    /// Called when evolution completes
    fn on_complete(&mut self, _champion: Option<&Individual>) {}
}

/// No-op callback
pub struct NoopCallback;
impl EvolutionCallback for NoopCallback {
    fn on_generation(&mut self, _stats: &GenerationStats, _population: &Population) {}
}

/// Evolution engine
pub struct EvolutionEngine {
    /// Parameters
    params: EvolutionParams,

    /// Population
    population: Population,

    /// Random number generator
    rng: StdRng,

    /// Generation statistics history
    history: Vec<GenerationStats>,

    /// Best fitness seen
    best_fitness: f32,

    /// Generations since improvement
    stagnation_counter: usize,
}

impl EvolutionEngine {
    /// Create new evolution engine
    pub fn new(params: EvolutionParams) -> Self {
        let seed = params.seed.unwrap_or_else(rand::random);
        let rng = StdRng::seed_from_u64(seed);

        let pop_params = PopulationParams {
            size: params.population_size,
            num_inputs: params.num_inputs,
            num_outputs: params.num_outputs,
            ..Default::default()
        };

        let population = Population::new(
            pop_params,
            params.mutation.clone(),
            params.crossover.clone(),
            params.speciation.clone(),
        );

        Self {
            params,
            population,
            rng,
            history: Vec::new(),
            best_fitness: f32::NEG_INFINITY,
            stagnation_counter: 0,
        }
    }

    /// Run evolution with fitness function
    pub fn run<F: FitnessFunction, C: EvolutionCallback>(
        &mut self,
        fitness_fn: &F,
        callback: &mut C,
    ) -> Result<Option<Genome>> {
        for _ in 0..self.params.max_generations {
            // Evaluate population
            self.population.evaluate(fitness_fn);

            // Compute statistics
            let stats = self.compute_stats();

            // Check for improvement
            if stats.best_fitness > self.best_fitness {
                self.best_fitness = stats.best_fitness;
                self.stagnation_counter = 0;
            } else {
                self.stagnation_counter += 1;
            }

            // Store history
            self.history.push(stats.clone());

            // Callback
            callback.on_generation(&stats, &self.population);

            // Check target
            if let Some(target) = self.params.target_fitness {
                if stats.best_fitness >= target {
                    if let Some(best) = self.population.best() {
                        callback.on_target_reached(&best.genome, stats.generation);
                        callback.on_complete(Some(best));
                        return Ok(Some(best.genome.clone()));
                    }
                }
            }

            // Check stagnation
            if self.stagnation_counter >= self.params.stagnation_limit {
                callback.on_stagnation(stats.generation);
                // Could add population reset or other recovery strategies here
            }

            // Evolve
            self.population.evolve(&mut self.rng, self.params.elitism)?;
        }

        // Evolution complete
        let champion = self.population.champion().cloned();
        callback.on_complete(champion.as_ref());

        Ok(champion.map(|c| c.genome))
    }

    /// Run for a single generation
    pub fn step<F: FitnessFunction>(&mut self, fitness_fn: &F) -> Result<GenerationStats> {
        self.population.evaluate(fitness_fn);
        let stats = self.compute_stats();
        self.history.push(stats.clone());
        self.population.evolve(&mut self.rng, self.params.elitism)?;
        Ok(stats)
    }

    /// Compute generation statistics
    fn compute_stats(&self) -> GenerationStats {
        let individuals = self.population.individuals();

        if individuals.is_empty() {
            return GenerationStats {
                generation: self.population.generation(),
                ..Default::default()
            };
        }

        let fitnesses: Vec<f32> = individuals.iter().map(|i| i.genome.fitness).collect();
        let best = fitnesses.iter().copied().fold(f32::NEG_INFINITY, f32::max);
        let worst = fitnesses.iter().copied().fold(f32::INFINITY, f32::min);
        let avg = fitnesses.iter().sum::<f32>() / fitnesses.len() as f32;

        let complexities: Vec<f32> = individuals
            .iter()
            .map(|i| (i.genome.nodes.len() + i.genome.connections.len()) as f32)
            .collect();
        let avg_complexity = complexities.iter().sum::<f32>() / complexities.len() as f32;

        GenerationStats {
            generation: self.population.generation(),
            best_fitness: best,
            avg_fitness: avg,
            worst_fitness: worst,
            num_species: self.population.num_species(),
            avg_complexity,
            num_innovations: self.population.innovation_tracker().innovation_count(),
        }
    }

    // ========== Getters ==========

    /// Get current population
    pub fn population(&self) -> &Population {
        &self.population
    }

    /// Get mutable population
    pub fn population_mut(&mut self) -> &mut Population {
        &mut self.population
    }

    /// Get evolution history
    pub fn history(&self) -> &[GenerationStats] {
        &self.history
    }

    /// Get current generation
    pub fn generation(&self) -> usize {
        self.population.generation()
    }

    /// Get best fitness seen
    pub fn best_fitness(&self) -> f32 {
        self.best_fitness
    }

    /// Get champion genome
    pub fn champion(&self) -> Option<&Individual> {
        self.population.champion()
    }
}

/// Quick evolution helper
pub fn evolve<F: FitnessFunction>(
    num_inputs: usize,
    num_outputs: usize,
    fitness_fn: &F,
    max_generations: usize,
) -> Result<Option<Genome>> {
    let params = EvolutionParams {
        num_inputs,
        num_outputs,
        max_generations,
        ..Default::default()
    };

    let mut engine = EvolutionEngine::new(params);
    let mut callback = NoopCallback;

    engine.run(fitness_fn, &mut callback)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::fitness::SimpleFitness;

    #[test]
    fn test_evolution_engine() {
        let params = EvolutionParams {
            population_size: 20,
            max_generations: 10,
            num_inputs: 2,
            num_outputs: 1,
            ..Default::default()
        };

        let mut engine = EvolutionEngine::new(params);
        let fitness = SimpleFitness::new(|g: &Genome| g.connections.len() as f32);
        let mut callback = NoopCallback;

        let result = engine.run(&fitness, &mut callback);
        assert!(result.is_ok());

        // Should have history
        assert_eq!(engine.history().len(), 10);
    }

    #[test]
    fn test_quick_evolve() {
        let fitness = SimpleFitness::new(|_: &Genome| 1.0);
        let result = evolve(2, 1, &fitness, 5);
        assert!(result.is_ok());
    }
}
