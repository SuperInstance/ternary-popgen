# ternary-popgen

**Population genetics for ternary agent systems. Fitness, selection, mutation, migration.**

Population genetics studies how allele frequencies change over time under evolutionary forces: natural selection (some alleles are fitter), genetic drift (random sampling), mutation (spontaneous change), and migration (gene flow between populations). This crate implements all four forces for ternary populations where every individual carries one of three alleles: `Negative (-1)`, `Neutral (0)`, or `Positive (+1)`.

## What's Inside

- **`Allele`** — the three ternary alleles: `Negative`, `Neutral`, `Positive`
- **`FitnessLandscape`** — fitness values for each allele. Presets: `neutral()`, `directional()`, `disruptive()`, `stabilizing()`
- **`Population`** — a collection of individuals. Initialize random or structured
- **`Selection`** — fitness-proportional selection. Fitter alleles reproduce more
- **`Mutation`** — spontaneous allele changes at a given rate per generation
- **`Migration`** — gene flow between populations. Exchange individuals at a given rate
- **`HardyWeinberg`** — test whether a population is at Hardy-Weinberg equilibrium
- **`allele_frequencies(pop)`** — current frequency of each allele
- **`heterozygosity(pop)`** — genetic diversity measure

## Quick Example

```rust
use ternary_popgen::*;

// Create a population of 100 agents
let mut pop = Population::random(100, &mut || 0.42);

// Define fitness landscape: positive allele is fittest
let landscape = FitnessLandscape::directional(1.5);

// Evolve for 100 generations
for _ in 0..100 {
    pop = pop.select(&landscape, &mut || 0.42);
    pop = pop.mutate(0.01, &mut || 0.42);
}

// Check allele frequencies
let freqs = allele_frequencies(&pop);
println!("Positive allele frequency: {:.2}", freqs[2]);
// Should be high — selection favored it

// Hardy-Weinberg equilibrium test
let hw = HardyWeinberg::test(&pop);
println!("At equilibrium: {}", hw.is_equilibrium(0.05));
```

## The Deeper Truth

**The neutral allele (0) is the pivot.** In a directional fitness landscape, the neutral allele is intermediate — less fit than positive, more fit than negative. It acts as a reservoir of genetic variation: under stabilizing selection, it's the *fittest* allele (everything converges to neutral). Under disruptive selection, it's the *least fit* (both extremes are favored). The neutral allele is to population genetics what the 0 state is to ternary dynamics in general: the pivot that everything else turns around.

The Hardy-Weinberg test is particularly informative for ternary populations. With three alleles, the equilibrium genotype frequencies follow a clear pattern: P(AA) = p², P(AB) = 2pq, for all three allele pairs. Deviations from HW equilibrium reveal which evolutionary forces are active.

**Use cases:**
- **Evolutionary computation** — fitness landscapes and selection for genetic algorithms
- **Population modeling** — simulate genetic change over time
- **Conservation genetics** — assess genetic diversity and extinction risk
- **Multi-agent systems** — agents evolve their ternary strategies over generations
- **Education** — the simplest non-binary population genetics framework

## See Also

- **ternary-drift** — genetic drift (the random force, without selection)
- **ternary-experiment** — parameter sweeps over fitness landscapes
- **ternary-ga** — genetic algorithms (optimization using these forces)
- **ternary-grace** — the finding that forgiveness (grace) kills, trust heals

## Install

```bash
cargo add ternary-popgen
```

## License

MIT
