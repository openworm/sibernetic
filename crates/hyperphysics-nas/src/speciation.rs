//! Speciation
//!
//! Groups genomes into species to protect topological innovations.

#[cfg(feature = "serde")]
use serde::{Deserialize, Serialize};

use crate::genome::Genome;

/// Speciation parameters
#[derive(Debug, Clone)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct SpeciationParams {
    /// Compatibility threshold
    pub threshold: f32,

    /// Coefficient for excess genes
    pub c1: f32,

    /// Coefficient for disjoint genes
    pub c2: f32,

    /// Coefficient for weight differences
    pub c3: f32,

    /// Target number of species
    pub target_species: usize,

    /// Threshold adjustment rate
    pub threshold_adjustment: f32,

    /// Minimum species size before elimination
    pub min_species_size: usize,

    /// Generations without improvement before stagnation
    pub stagnation_threshold: usize,
}

impl Default for SpeciationParams {
    fn default() -> Self {
        Self {
            threshold: 3.0,
            c1: 1.0,
            c2: 1.0,
            c3: 0.4,
            target_species: 10,
            threshold_adjustment: 0.1,
            min_species_size: 2,
            stagnation_threshold: 15,
        }
    }
}

impl SpeciationParams {
    /// More protective of species
    pub fn protective() -> Self {
        Self {
            threshold: 2.0,
            target_species: 15,
            stagnation_threshold: 20,
            ..Default::default()
        }
    }

    /// More aggressive species merging
    pub fn aggressive() -> Self {
        Self {
            threshold: 4.0,
            target_species: 5,
            stagnation_threshold: 10,
            ..Default::default()
        }
    }
}

/// A species of similar genomes
#[derive(Debug, Clone)]
pub struct Species {
    /// Species ID
    pub id: usize,

    /// Representative genome (for compatibility testing)
    pub representative: Genome,

    /// Member genome indices
    pub members: Vec<usize>,

    /// Best fitness in this species
    pub best_fitness: f32,

    /// Average fitness
    pub avg_fitness: f32,

    /// Generations since fitness improved
    pub stagnation: usize,

    /// Offspring allocation
    pub offspring_count: usize,
}

impl Species {
    pub fn new(id: usize, representative: Genome) -> Self {
        Self {
            id,
            representative,
            members: Vec::new(),
            best_fitness: f32::NEG_INFINITY,
            avg_fitness: 0.0,
            stagnation: 0,
            offspring_count: 0,
        }
    }

    /// Check if a genome is compatible with this species
    pub fn is_compatible(&self, genome: &Genome, params: &SpeciationParams) -> bool {
        let distance = self.representative.distance(genome, params.c1, params.c2, params.c3);
        distance < params.threshold
    }

    /// Update species statistics
    pub fn update_stats(&mut self, fitnesses: &[f32]) {
        if fitnesses.is_empty() {
            return;
        }

        let max_fitness = fitnesses.iter().copied().fold(f32::NEG_INFINITY, f32::max);
        let sum: f32 = fitnesses.iter().sum();
        self.avg_fitness = sum / fitnesses.len() as f32;

        if max_fitness > self.best_fitness {
            self.best_fitness = max_fitness;
            self.stagnation = 0;
        } else {
            self.stagnation += 1;
        }
    }

    /// Check if species is stagnant
    pub fn is_stagnant(&self, threshold: usize) -> bool {
        self.stagnation >= threshold
    }
}

/// Manages species throughout evolution
pub struct SpeciesManager {
    /// All species
    species: Vec<Species>,

    /// Parameters
    params: SpeciationParams,

    /// Next species ID
    next_id: usize,

    /// Current compatibility threshold (adaptive)
    current_threshold: f32,
}

impl SpeciesManager {
    pub fn new(params: SpeciationParams) -> Self {
        let current_threshold = params.threshold;
        Self {
            species: Vec::new(),
            params,
            next_id: 0,
            current_threshold,
        }
    }

    /// Speciate a population
    pub fn speciate(&mut self, population: &mut [Genome]) {
        // Clear existing members
        for species in &mut self.species {
            species.members.clear();
        }

        // Assign each genome to a species
        for (idx, genome) in population.iter_mut().enumerate() {
            let mut found = false;

            for species in &mut self.species {
                if species.is_compatible(genome, &self.params) {
                    species.members.push(idx);
                    genome.species_id = species.id;
                    found = true;
                    break;
                }
            }

            if !found {
                // Create new species
                let mut new_species = Species::new(self.next_id, genome.clone());
                new_species.members.push(idx);
                genome.species_id = self.next_id;
                self.species.push(new_species);
                self.next_id += 1;
            }
        }

        // Remove empty species
        self.species.retain(|s| !s.members.is_empty());

        // Update representatives for next generation
        for species in &mut self.species {
            if !species.members.is_empty() {
                let rep_idx = species.members[0];
                species.representative = population[rep_idx].clone();
            }
        }

        // Adjust threshold to target number of species
        self.adjust_threshold();
    }

    /// Adjust compatibility threshold
    fn adjust_threshold(&mut self) {
        let num_species = self.species.len();
        let target = self.params.target_species;

        if num_species < target {
            self.current_threshold -= self.params.threshold_adjustment;
        } else if num_species > target {
            self.current_threshold += self.params.threshold_adjustment;
        }

        self.current_threshold = self.current_threshold.max(0.5);
    }

    /// Update species statistics and compute adjusted fitness
    pub fn update_fitness(&mut self, population: &mut [Genome]) {
        for species in &mut self.species {
            let fitnesses: Vec<f32> = species
                .members
                .iter()
                .map(|&idx| population[idx].fitness)
                .collect();

            species.update_stats(&fitnesses);

            // Compute adjusted fitness (fitness sharing)
            let species_size = species.members.len() as f32;
            for &idx in &species.members {
                population[idx].adjusted_fitness = population[idx].fitness / species_size;
            }
        }
    }

    /// Allocate offspring to species
    pub fn allocate_offspring(&mut self, total_offspring: usize) {
        // Calculate total adjusted fitness
        let total_adj_fitness: f32 = self
            .species
            .iter()
            .map(|s| s.avg_fitness.max(0.0))
            .sum();

        if total_adj_fitness <= 0.0 {
            // Equal distribution
            let per_species = total_offspring / self.species.len().max(1);
            for species in &mut self.species {
                species.offspring_count = per_species;
            }
            return;
        }

        // Proportional allocation
        let mut allocated = 0;
        for species in &mut self.species {
            let proportion = species.avg_fitness.max(0.0) / total_adj_fitness;
            species.offspring_count = (proportion * total_offspring as f32).floor() as usize;
            allocated += species.offspring_count;
        }

        // Distribute remaining
        let remaining = total_offspring.saturating_sub(allocated);
        for (i, species) in self.species.iter_mut().enumerate() {
            if i < remaining {
                species.offspring_count += 1;
            }
        }
    }

    /// Remove stagnant species
    pub fn remove_stagnant(&mut self) {
        self.species.retain(|s| {
            !s.is_stagnant(self.params.stagnation_threshold)
                || s.members.len() >= self.params.min_species_size
        });
    }

    /// Get all species
    pub fn species(&self) -> &[Species] {
        &self.species
    }

    /// Get number of species
    pub fn num_species(&self) -> usize {
        self.species.len()
    }

    /// Get current threshold
    pub fn threshold(&self) -> f32 {
        self.current_threshold
    }

    /// Get species by ID
    pub fn get_species(&self, id: usize) -> Option<&Species> {
        self.species.iter().find(|s| s.id == id)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::genome::InnovationTracker;

    #[test]
    fn test_speciation() {
        let mut tracker = InnovationTracker::new();
        let mut population: Vec<Genome> = (0..10)
            .map(|_| Genome::minimal(2, 1, &mut tracker))
            .collect();

        let mut manager = SpeciesManager::new(SpeciationParams::default());
        manager.speciate(&mut population);

        // All identical genomes should be in one species
        assert_eq!(manager.num_species(), 1);
    }

    #[test]
    fn test_offspring_allocation() {
        let mut tracker = InnovationTracker::new();
        let mut population: Vec<Genome> = (0..10)
            .map(|i| {
                let mut g = Genome::minimal(2, 1, &mut tracker);
                g.fitness = i as f32;
                g
            })
            .collect();

        let mut manager = SpeciesManager::new(SpeciationParams::default());
        manager.speciate(&mut population);
        manager.update_fitness(&mut population);
        manager.allocate_offspring(10);

        let total: usize = manager.species().iter().map(|s| s.offspring_count).sum();
        assert_eq!(total, 10);
    }
}
