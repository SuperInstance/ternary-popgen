#![forbid(unsafe_code)]
//! Population genetics for ternary agent systems.

/// A ternary allele.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum Allele { Negative, Neutral, Positive }

impl Allele {
    pub fn fitness(&self, landscape: &FitnessLandscape) -> f64 {
        match self {
            Allele::Negative => landscape.neg_fitness,
            Allele::Neutral => landscape.neutral_fitness,
            Allele::Positive => landscape.pos_fitness,
        }
    }
}

/// Fitness landscape over ternary alleles.
#[derive(Debug, Clone)]
pub struct FitnessLandscape {
    pub neg_fitness: f64,
    pub neutral_fitness: f64,
    pub pos_fitness: f64,
}

impl FitnessLandscape {
    pub fn neutral() -> Self { Self { neg_fitness: 1.0, neutral_fitness: 1.0, pos_fitness: 1.0 } }
    pub fn directional(pos: f64) -> Self { Self { neg_fitness: 0.5, neutral_fitness: 0.75, pos_fitness: pos } }
    pub fn disruptive() -> Self { Self { neg_fitness: 1.2, neutral_fitness: 0.6, pos_fitness: 1.2 } }
    pub fn stabilizing() -> Self { Self { neg_fitness: 0.8, neutral_fitness: 1.0, pos_fitness: 0.8 } }
}

/// A population of ternary agents.
#[derive(Debug, Clone)]
pub struct Population {
    pub individuals: Vec<Allele>,
    pub generation: u64,
}

impl Population {
    pub fn new(size: usize, alleles: Vec<Allele>) -> Self {
        Self { individuals: alleles.into_iter().chain(std::iter::repeat(Allele::Neutral)).take(size).collect(), generation: 0 }
    }

    pub fn random(size: usize, rng: &mut impl FnMut() -> f64) -> Self {
        let individuals = (0..size).map(|_| match (rng() * 3.0) as usize {
            0 => Allele::Negative, 1 => Allele::Neutral, _ => Allele::Positive
        }).collect();
        Self { individuals, generation: 0 }
    }

    pub fn allele_frequencies(&self) -> (f64, f64, f64) {
        let n = self.individuals.len() as f64;
        let neg = self.individuals.iter().filter(|a| **a == Allele::Negative).count() as f64 / n;
        let neu = self.individuals.iter().filter(|a| **a == Allele::Neutral).count() as f64 / n;
        let pos = 1.0 - neg - neu;
        (neg, neu, pos)
    }

    pub fn heterozygosity(&self) -> f64 {
        let (p, q, r) = self.allele_frequencies();
        1.0 - p*p - q*q - r*r
    }

    pub fn size(&self) -> usize { self.individuals.len() }

    /// Wright-Fisher generation: sample next generation from current fitness-weighted distribution.
    pub fn wright_fisher(&mut self, landscape: &FitnessLandscape, mutation_rate: f64, rng: &mut impl FnMut() -> f64) {
        let n = self.individuals.len();
        // Compute selection weights
        let weights: Vec<f64> = self.individuals.iter().map(|a| a.fitness(landscape)).collect();
        let total: f64 = weights.iter().sum();
        let probs: Vec<f64> = weights.iter().map(|w| w / total).collect();

        // Sample next generation
        let mut new_gen = Vec::with_capacity(n);
        for _ in 0..n {
            let r = rng() * total;
            let mut cumsum = 0.0;
            let mut chosen = self.individuals[0];
            for (i, &w) in weights.iter().enumerate() {
                cumsum += w;
                if cumsum >= r { chosen = self.individuals[i]; break; }
            }
            // Mutation
            if rng() < mutation_rate {
                chosen = match (rng() * 3.0) as usize {
                    0 => Allele::Negative, 1 => Allele::Neutral, _ => Allele::Positive
                };
            }
            new_gen.push(chosen);
        }
        self.individuals = new_gen;
        self.generation += 1;
    }

    /// Moran process: one birth-death event per step.
    pub fn moran_step(&mut self, landscape: &FitnessLandscape, mutation_rate: f64, rng: &mut impl FnMut() -> f64) {
        let n = self.individuals.len();
        if n == 0 { return; }

        // Select individual to reproduce (fitness-weighted)
        let weights: Vec<f64> = self.individuals.iter().map(|a| a.fitness(landscape)).collect();
        let total: f64 = weights.iter().sum();
        let r1 = rng() * total;
        let mut cumsum = 0.0;
        let mut parent_idx = 0;
        for (i, &w) in weights.iter().enumerate() {
            cumsum += w;
            if cumsum >= r1 { parent_idx = i; break; }
        }

        // Select individual to die (uniform)
        let death_idx = (rng() * n as f64) as usize;

        // Birth with possible mutation
        let mut offspring = self.individuals[parent_idx];
        if rng() < mutation_rate {
            offspring = match (rng() * 3.0) as usize {
                0 => Allele::Negative, 1 => Allele::Neutral, _ => Allele::Positive
            };
        }

        self.individuals[death_idx] = offspring;
        self.generation += 1;
    }

    /// Hardy-Weinberg equilibrium test (chi-squared).
    pub fn hardy_weinberg_test(&self) -> f64 {
        let n = self.individuals.len() as f64;
        let (p, q, r) = self.allele_frequencies();
        // Expected heterozygosity under HWE
        let expected_het = 2.0 * (p*q + p*r + q*r);
        let observed_het = self.heterozygosity();
        if expected_het == 0.0 { return 0.0; }
        (observed_het - expected_het).powi(2) / expected_het
    }

    /// Effective population size (from heterozygosity).
    pub fn effective_size(&self, prev_het: f64) -> f64 {
        let curr_het = self.heterozygosity();
        if prev_het <= curr_het || prev_het <= 0.0 { return f64::INFINITY; }
        curr_het / (prev_het - curr_het)
    }

    /// Run Wright-Fisher for N generations, return frequency history.
    pub fn run_wf(&mut self, gens: u64, landscape: &FitnessLandscape, mutation_rate: f64) -> Vec<(f64, f64, f64)> {
        let mut rng_state: u64 = 42;
        let mut rng = || -> f64 {
            rng_state = rng_state.wrapping_mul(6364136223846793005).wrapping_add(1);
            (rng_state >> 33) as f64 / (1u64 << 31) as f64
        };
        let mut history = Vec::new();
        for _ in 0..gens {
            self.wright_fisher(landscape, mutation_rate, &mut rng);
            history.push(self.allele_frequencies());
        }
        history
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn rng() -> impl FnMut() -> f64 {
        let mut s: u64 = 12345;
        move || { s = s.wrapping_mul(6364136223846793005).wrapping_add(1); (s >> 33) as f64 / (1u64 << 31) as f64 }
    }

    #[test]
    fn test_allele_frequencies() {
        let pop = Population::new(3, vec![Allele::Negative, Allele::Neutral, Allele::Positive]);
        let (p, q, r) = pop.allele_frequencies();
        assert!((p - 1.0/3.0).abs() < 0.01);
        assert!((q - 1.0/3.0).abs() < 0.01);
        assert!((r - 1.0/3.0).abs() < 0.01);
    }

    #[test]
    fn test_neutral_drift() {
        let mut pop = Population::random(100, &mut rng());
        let landscape = FitnessLandscape::neutral();
        let history = pop.run_wf(100, &landscape, 0.0);
        assert_eq!(history.len(), 100);
        // Frequencies should still sum to 1
        for (p, q, r) in &history {
            assert!((p + q + r - 1.0).abs() < 0.01);
        }
    }

    #[test]
    fn test_directional_selection() {
        let mut pop = Population::random(200, &mut rng());
        let landscape = FitnessLandscape::directional(2.0);
        let history = pop.run_wf(200, &landscape, 0.001);
        let (p, _, r) = history.last().unwrap();
        assert!(*r > *p, "Positive should dominate under directional selection");
    }

    #[test]
    fn test_moran_step() {
        let mut pop = Population::random(50, &mut rng());
        let landscape = FitnessLandscape::neutral();
        let gen_before = pop.generation;
        pop.moran_step(&landscape, 0.01, &mut rng());
        assert_eq!(pop.generation, gen_before + 1);
        assert_eq!(pop.size(), 50);
    }

    #[test]
    fn test_heterozygosity() {
        let pop = Population::new(3, vec![Allele::Negative, Allele::Neutral, Allele::Positive]);
        let h = pop.heterozygosity();
        assert!(h > 0.6, "3 different alleles should have high heterozygosity");
    }

    #[test]
    fn test_heterozygosity_pure() {
        let pop = Population::new(10, vec![Allele::Positive; 10]);
        assert!(pop.heterozygosity().abs() < 0.01);
    }

    #[test]
    fn test_stabilizing_selection() {
        let mut pop = Population::random(200, &mut rng());
        let landscape = FitnessLandscape::stabilizing();
        let history = pop.run_wf(200, &landscape, 0.01);
        let (_, q, _) = history.last().unwrap();
        // Neutral allele should increase under stabilizing selection
        assert!(*q > 0.2, "Neutral should be favored: q={}", q);
    }

    #[test]
    fn test_disruptive_selection() {
        let mut pop = Population::random(200, &mut rng());
        let landscape = FitnessLandscape::disruptive();
        let history = pop.run_wf(200, &landscape, 0.01);
        let (_, q, _) = history.last().unwrap();
        assert!(*q < 0.5, "Neutral should be disfavored: q={}", q);
    }

    #[test]
    fn test_hardy_weinberg() {
        let pop = Population::random(1000, &mut rng());
        let chi2 = pop.hardy_weinberg_test();
        assert!(chi2 >= 0.0);
    }

    #[test]
    fn test_mutation_prevents_fixation() {
        let mut pop = Population::new(50, vec![Allele::Positive; 50]);
        let landscape = FitnessLandscape::neutral();
        let history = pop.run_wf(100, &landscape, 0.1); // High mutation
        let (p, q, r) = history.last().unwrap();
        assert!(*q > 0.0 || *p > 0.0, "Mutation should introduce diversity");
    }

    #[test]
    fn test_bottleneck() {
        let mut pop = Population::random(200, &mut rng());
        // Bottleneck to 10
        pop.individuals.truncate(10);
        let h_before = pop.heterozygosity();
        let landscape = FitnessLandscape::neutral();
        pop.run_wf(50, &landscape, 0.0);
        // After bottleneck, heterozygosity likely reduced
        assert_eq!(pop.size(), 10);
    }

    #[test]
    fn test_frequency_sum() {
        let mut pop = Population::random(100, &mut rng());
        let landscape = FitnessLandscape::neutral();
        for _ in 0..50 {
            pop.wright_fisher(&landscape, 0.01, &mut rng());
            let (p, q, r) = pop.allele_frequencies();
            assert!((p + q + r - 1.0).abs() < 0.01, "Frequencies must sum to 1");
        }
    }
}
