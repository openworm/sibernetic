//! Population Management
//!
//! Maintains and evolves a population of genomes.

#[cfg(feature = "serde")]
use serde::{Deserialize, Serialize};

use rand::Rng;

use crate::crossover::{CrossoverOperator, CrossoverParams};
use crate::fitness::{FitnessFunction, MultiFitness};
use crate::genome::{Genome, InnovationTracker};
use crate::mutation::{MutationOperator, MutationParams};
use crate::speciation::{SpeciationParams, SpeciesManager};
use crate::Result;

/// Population parameters
#[derive(Debug, Clone)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct PopulationParams {
    /// Population size
    pub size: usize,

    /// Number of inputs
    pub num_inputs: usize,

    /// Number of outputs
    pub num_outputs: usize,

    /// Probability of crossover (vs asexual reproduction)
    pub crossover_prob: f32,

    /// Only allow crossover within species
    pub interspecies_mating: bool,

    /// Probability of interspecies mating (if allowed)
    pub interspecies_prob: f32,
}

impl Default for PopulationParams {
    fn default() -> Self {
        Self {
            size: 100,
            num_inputs: 3,
            num_outputs: 1,
            crossover_prob: 0.75,
            interspecies_mating: true,
            interspecies_prob: 0.001,
        }
    }
}

/// An individual in the population
#[derive(Debug, Clone)]
pub struct Individual {
    /// Genome
    pub genome: Genome,

    /// Multi-objective fitness
    pub fitness: MultiFitness,

    /// Generation born
    pub generation: usize,
}

impl Individual {
    pub fn new(genome: Genome, generation: usize) -> Self {
        Self {
            genome,
            fitness: MultiFitness::new(),
            generation,
        }
    }
}

/// The evolving population
pub struct Population {
    /// Population parameters
    params: PopulationParams,

    /// Current individuals
    individuals: Vec<Individual>,

    /// Innovation tracker
    innovation_tracker: InnovationTracker,

    /// Species manager
    species_manager: SpeciesManager,

    /// Mutation operator
    mutation: MutationOperator,

    /// Crossover operator
    crossover: CrossoverOperator,

    /// Current generation
    generation: usize,

    /// Best individual ever
    champion: Option<Individual>,
}

impl Population {
    /// Create new population
    pub fn new(
        params: PopulationParams,
        mutation_params: MutationParams,
        crossover_params: CrossoverParams,
        speciation_params: SpeciationParams,
    ) -> Self {
        let num_initial_nodes = params.num_inputs + 1 + params.num_outputs;
        let mut innovation_tracker = InnovationTracker::with_nodes(num_initial_nodes as u32);

        // Create initial population
        let individuals: Vec<Individual> = (0..params.size)
            .map(|_| {
                let genome = Genome::minimal(params.num_inputs, params.num_outputs, &mut innovation_tracker);
                Individual::new(genome, 0)
            })
            .collect();

        Self {
            params,
            individuals,
            innovation_tracker,
            species_manager: SpeciesManager::new(speciation_params),
            mutation: MutationOperator::new(mutation_params),
            crossover: CrossoverOperator::new(crossover_params),
            generation: 0,
            champion: None,
        }
    }

    /// Evaluate all individuals using fitness function
    pub fn evaluate<F: FitnessFunction>(&mut self, fitness_fn: &F) {
        for individual in &mut self.individuals {
            individual.fitness = fitness_fn.evaluate_multi(&individual.genome);
            individual.genome.fitness = individual.fitness.primary();
        }

        // Update champion
        if let Some(best) = self.individuals.iter().max_by(|a, b| {
            a.genome.fitness.partial_cmp(&b.genome.fitness).unwrap_or(std::cmp::Ordering::Equal)
        }) {
            if self.champion.is_none() || best.genome.fitness > self.champion.as_ref().unwrap().genome.fitness {
                self.champion = Some(best.clone());
            }
        }
    }

    /// Perform speciation
    pub fn speciate(&mut self) {
        let mut genomes: Vec<Genome> = self.individuals.iter().map(|i| i.genome.clone()).collect();
        self.species_manager.speciate(&mut genomes);
        self.species_manager.update_fitness(&mut genomes);

        // Update individuals with speciation results
        for (individual, genome) in self.individuals.iter_mut().zip(genomes.into_iter()) {
            individual.genome = genome;
        }
    }

    /// Create next generation
    pub fn evolve<R: Rng>(&mut self, rng: &mut R, elitism: usize) -> Result<()> {
        self.speciate();

        // Allocate offspring to species
        let target_size = self.params.size.saturating_sub(elitism);
        self.species_manager.allocate_offspring(target_size);

        // Collect elite individuals
        let mut elites: Vec<Individual> = Vec::new();
        if elitism > 0 {
            let mut sorted: Vec<_> = self.individuals.iter().cloned().collect();
            sorted.sort_by(|a, b| {
                b.genome.fitness.partial_cmp(&a.genome.fitness).unwrap_or(std::cmp::Ordering::Equal)
            });
            elites = sorted.into_iter().take(elitism).collect();
        }

        // Create offspring
        let mut offspring: Vec<Individual> = Vec::new();

        for species in self.species_manager.species() {
            if species.offspring_count == 0 {
                continue;
            }

            // Get species members
            let members: Vec<&Individual> = species
                .members
                .iter()
                .map(|&idx| &self.individuals[idx])
                .collect();

            if members.is_empty() {
                continue;
            }

            for _ in 0..species.offspring_count {
                let child_genome = if rng.gen::<f32>() < self.params.crossover_prob && members.len() > 1 {
                    // Sexual reproduction
                    let parent1_idx = self.select_parent(&members, rng);
                    let parent2_idx = self.select_parent(&members, rng);

                    let parent1 = &members[parent1_idx].genome;
                    let parent2 = &members[parent2_idx].genome;

                    // Order by fitness
                    let (better, worse) = if parent1.fitness >= parent2.fitness {
                        (parent1, parent2)
                    } else {
                        (parent2, parent1)
                    };

                    let mut child = self.crossover.crossover(better, worse, rng);
                    self.mutation.mutate(&mut child, &mut self.innovation_tracker, rng);
                    child
                } else {
                    // Asexual reproduction
                    let parent_idx = self.select_parent(&members, rng);
                    let mut child = members[parent_idx].genome.clone();
                    self.mutation.mutate(&mut child, &mut self.innovation_tracker, rng);
                    child
                };

                offspring.push(Individual::new(child_genome, self.generation + 1));
            }
        }

        // Combine elites and offspring
        self.individuals = elites;
        self.individuals.extend(offspring);

        // Pad if necessary
        while self.individuals.len() < self.params.size {
            let genome = Genome::minimal(
                self.params.num_inputs,
                self.params.num_outputs,
                &mut self.innovation_tracker,
            );
            self.individuals.push(Individual::new(genome, self.generation + 1));
        }

        // Truncate if necessary
        self.individuals.truncate(self.params.size);

        self.generation += 1;

        // Remove stagnant species
        self.species_manager.remove_stagnant();

        Ok(())
    }

    /// Tournament selection
    fn select_parent<'a, R: Rng>(&self, candidates: &[&'a Individual], rng: &mut R) -> usize {
        if candidates.len() < 2 {
            return 0;
        }

        let tournament_size = 3.min(candidates.len());
        let mut best_idx = rng.gen_range(0..candidates.len());
        let mut best_fitness = candidates[best_idx].genome.fitness;

        for _ in 1..tournament_size {
            let idx = rng.gen_range(0..candidates.len());
            if candidates[idx].genome.fitness > best_fitness {
                best_idx = idx;
                best_fitness = candidates[idx].genome.fitness;
            }
        }

        best_idx
    }

    // ========== Getters ==========

    /// Get current generation
    pub fn generation(&self) -> usize {
        self.generation
    }

    /// Get population size
    pub fn size(&self) -> usize {
        self.individuals.len()
    }

    /// Get best individual in current population
    pub fn best(&self) -> Option<&Individual> {
        self.individuals.iter().max_by(|a, b| {
            a.genome.fitness.partial_cmp(&b.genome.fitness).unwrap_or(std::cmp::Ordering::Equal)
        })
    }

    /// Get all-time champion
    pub fn champion(&self) -> Option<&Individual> {
        self.champion.as_ref()
    }

    /// Get average fitness
    pub fn avg_fitness(&self) -> f32 {
        if self.individuals.is_empty() {
            return 0.0;
        }
        let sum: f32 = self.individuals.iter().map(|i| i.genome.fitness).sum();
        sum / self.individuals.len() as f32
    }

    /// Get number of species
    pub fn num_species(&self) -> usize {
        self.species_manager.num_species()
    }

    /// Get all individuals
    pub fn individuals(&self) -> &[Individual] {
        &self.individuals
    }

    /// Get innovation tracker
    pub fn innovation_tracker(&self) -> &InnovationTracker {
        &self.innovation_tracker
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::fitness::SimpleFitness;

    #[test]
    fn test_population_creation() {
        let pop = Population::new(
            PopulationParams {
                size: 10,
                num_inputs: 2,
                num_outputs: 1,
                ..Default::default()
            },
            MutationParams::default(),
            CrossoverParams::default(),
            SpeciationParams::default(),
        );

        assert_eq!(pop.size(), 10);
        assert_eq!(pop.generation(), 0);
    }

    #[test]
    fn test_evolution() {
        let mut rng = rand::thread_rng();

        let mut pop = Population::new(
            PopulationParams {
                size: 20,
                num_inputs: 2,
                num_outputs: 1,
                ..Default::default()
            },
            MutationParams::exploratory(),
            CrossoverParams::default(),
            SpeciationParams::default(),
        );

        // Simple fitness: more connections = better
        let fitness = SimpleFitness::new(|g: &Genome| g.connections.len() as f32);

        for _ in 0..5 {
            pop.evaluate(&fitness);
            pop.evolve(&mut rng, 2).unwrap();
        }

        assert_eq!(pop.generation(), 5);
    }
}
