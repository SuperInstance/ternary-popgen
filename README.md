# ternary-popgen

Population genetics for ternary {-1, 0, +1} agents. Wright-Fisher model, Moran process, fitness landscapes (directional, stabilizing, disruptive), heterozygosity, Hardy-Weinberg equilibrium testing, and effective population size estimation.

## Why It Matters

Population genetics describes how allele frequencies change over time under the forces of selection, mutation, and drift. By mapping ternary agent states to alleles — Negative (-1), Neutral (0), Positive (+1) — we can apply the full machinery of theoretical population genetics to ternary agent systems:

- **Wright-Fisher model**: discrete-generation drift and selection
- **Moran process**: continuous-time birth-death dynamics
- **Fitness landscapes**: directional, stabilizing, disruptive, and neutral selection
- **Hardy-Weinberg equilibrium**: test for evolutionary forces
- **Effective population size**: quantify drift strength from heterozygosity loss

This lets us predict: Will the +1 allele go to fixation? How fast? Under what mutation rate is diversity maintained?

## How It Works

### Allele Model

Three alleles map directly to ternary values:

| Allele | Value | Biological Analog |
|---|---|---|
| `Negative` | -1 | Deleterious / loss-of-function |
| `Neutral` | 0 | Wild-type / silent |
| `Positive` | +1 | Advantageous / gain-of-function |

### Fitness Landscape

Each allele has a fitness value $w_a$ determining its reproductive success:

| Landscape | $w_{-1}$ | $w_0$ | $w_{+1}$ | Effect |
|---|---|---|---|---|
| Neutral | 1.0 | 1.0 | 1.0 | Drift only |
| Directional | 0.5 | 0.75 | $w^+$ | Selection favors +1 |
| Stabilizing | 0.8 | 1.0 | 0.8 | Selection favors 0 |
| Disruptive | 1.2 | 0.6 | 1.2 | Selection favors extremes |

### Wright-Fisher Model

Discrete generations. The entire population is replaced each generation by sampling from the fitness-weighted distribution:

$$P(a_i \to a_j) = \frac{w_{a_j}}{\sum_k w_{a_k}}$$

With mutation rate $\mu$, each offspring has probability $\mu$ of randomizing to any of the three alleles.

**Properties:**
- Genetic drift: $P(\text{fixation of neutral allele}) = f_0$ (initial frequency)
- Expected heterozygosity loss per generation: $\Delta H = -H / (2N_e)$

**Complexity:** O($N$) per generation for population size $N$.

### Moran Process

One individual reproduces (fitness-weighted) and one dies (uniform random) per step. The population size $N$ is constant.

**Time to fixation:** O($N^2$) steps for neutral alleles, O($N$) for strongly selected alleles.

**Complexity:** O($N$) per step (fitness-weighted sampling).

### Allele Frequencies and Heterozygosity

Frequencies: $p = f_{-1}$, $q = f_0$, $r = f_{+1}$, where $p + q + r = 1$.

Gene diversity (expected heterozygosity):

$$H = 1 - p^2 - q^2 - r^2 = 2(pq + pr + qr)$$

- $H = 0$: monomorphic (fixed population)
- $H = 2/3$: uniform distribution (maximum diversity for 3 alleles)

### Hardy-Weinberg Test

Under Hardy-Weinberg equilibrium (random mating, no selection), expected heterozygosity is:

$$H_{\text{exp}} = 2(pq + pr + qr)$$

Chi-squared statistic:

$$\chi^2 = \frac{(H_{\text{obs}} - H_{\text{exp}})^2}{H_{\text{exp}}}$$

A significant $\chi^2$ indicates selection, non-random mating, or population structure.

### Effective Population Size

From heterozygosity decline:

$$N_e = \frac{H_t}{H_0 - H_t}$$

where $H_0$ is the initial heterozygosity and $H_t$ is the current value.

## Quick Start

```rust
use ternary_popgen::*;

let mut rng = || -> f64 {
    // deterministic RNG closure
    # let mut s = 42u64;
    // s = s.wrapping_mul(6364136223846793005).wrapping_add(1);
    // (s >> 33) as f64 / (1u64 << 31) as f64
    0.5
};

// Create population
let mut pop = Population::new(100,
    vec![Allele::Positive, Allele::Neutral, Allele::Negative]);
let (p, q, r) = pop.allele_frequencies();

// Directional selection: Positive should fix
let landscape = FitnessLandscape::directional(2.0);
let history = pop.run_wf(200, &landscape, 0.001);
let (p_final, _, r_final) = *history.last().unwrap();
// r_final (Positive) > p_final (Negative)

// Stabilizing selection favors Neutral
let mut pop2 = Population::random(200, &mut rng);
let landscape2 = FitnessLandscape::stabilizing();
pop2.run_wf(100, &landscape2, 0.01);
let (_, q2, _) = pop2.allele_frequencies();
// q2 (Neutral frequency) should be high

// Moran process
let mut pop3 = Population::random(50, &mut rng);
pop3.moran_step(&FitnessLandscape::neutral(), 0.01, &mut rng);

// Hardy-Weinberg test
let chi2 = pop.hardy_weinberg_test();
println!("HWE chi² = {:.4}", chi2);
```

## API

| Type / Method | Description |
|---|---|
| `Allele::Negative / Neutral / Positive` | The three ternary alleles |
| `FitnessLandscape::neutral/directional/disruptive/stabilizing()` | Predefined selection regimes |
| `Population::new(size, alleles)` | Create initial population |
| `Population::random(size, rng)` | Random initial distribution |
| `.allele_frequencies() → (f64, f64, f64)` | Frequencies (p, q, r) |
| `.heterozygosity() → f64` | Gene diversity $H$ |
| `.wright_fisher(landscape, mutation_rate, rng)` | One WF generation |
| `.moran_step(landscape, mutation_rate, rng)` | One Moran event |
| `.hardy_weinberg_test() → f64` | χ² statistic |
| `.effective_size(prev_het) → f64` | $N_e$ from heterozygosity |
| `.run_wf(gens, landscape, rate) → Vec<(f64,f64,f64)>` | N-generation history |

## Architecture Notes

Population genetics is the natural domain of the **γ + η = C** conservation identity. The Positive allele (+1) represents constructive mass γ, the Negative allele (-1) represents inhibitory mass η, and the Neutral allele (0) is the ancestral/substrate state. The conserved total $C = N$ (population size) is fixed by the Moran constraint (birth = death) and by the Wright-Fisher resampling.

The key insight: allele frequencies $p + q + r = 1$ is exactly $\gamma/C + \eta/C + \text{neutral}/C = 1$. Selection changes the *equilibrium* distribution between γ and η, but mutation and drift control the *rate* of conversion. The fixation probability of a beneficial +1 allele is approximately:

$$P_{\text{fix}} \approx \frac{2s}{1 - e^{-2N_es}}$$

where $s = w_{+1}/w_0 - 1$ is the selection coefficient — directly measuring how much γ-mass the constructive allele accumulates per generation.

## References

- Hartl, D. L. & Clark, A. G. (2007). *Principles of Population Genetics.* 4th ed. Sinauer.
- Ewens, W. J. (2004). *Mathematical Population Genetics.* 2nd ed. Springer.
- Wright, S. (1931). *Evolution in Mendelian Populations.* Genetics, 16(2). ($N_e$ theory)
- Moran, P. A. P. (1958). *Random Processes in Genetics.* Proc. Cambridge Phil. Soc. (Moran process)

## License

MIT
