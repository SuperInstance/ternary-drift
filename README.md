# Ternary Drift

Genetic drift simulation for **ternary populations** — modeling the stochastic forces that change allele frequencies {-1, 0, +1} in finite populations through the Wright-Fisher model. Computes heterozygosity decay, fixation probabilities, effective population size, and diversity loss.

## Why It Matters

In infinite populations, allele frequencies are deterministic — governed by selection pressure alone. But real populations are finite. **Genetic drift** — the random fluctuation of allele frequencies due to sampling — is the dominant evolutionary force in small populations. It causes:

- **Loss of diversity**: Alleles randomly fix or go extinct
- **Differentiation**: Isolated populations diverge genetically
- **Inbreeding**: Effective population size shrinks over generations

For ternary populations {-1, 0, +1}, drift is especially consequential because the three-allele system has more fixation states than binary (3 possible fixations vs. 2). The rate of heterozygosity decay — the "molecular clock" of drift — follows:

$$H(t) = H_0 \left(1 - \frac{1}{N_e}\right)^t$$

where $N_e$ is the effective population size. This crate implements the full simulation loop and the analytical formulas.

## How It Works

### Wright-Fisher Model

Each generation, every individual samples a parent uniformly at random from the previous generation:

$$X_i^{(t+1)} = X_J^{(t)}, \quad J \sim \text{Uniform}\{1, \ldots, N\}$$

This is a Markov chain on $\{-1, 0, +1\}^N$ with absorbing states at fixation (all individuals have the same allele). The expected heterozygosity after one generation:

$$\mathbb{E}[H_{t+1}] = H_t \left(1 - \frac{1}{2N}\right)$$

For diploids: $\mathbb{E}[H_{t+1}] = H_t \left(1 - \frac{1}{2N_e}\right)$. We use the haploid approximation (factor of $N$, not $2N$).

### Heterozygosity

The expected heterozygosity (probability two random individuals differ):

$$H = 1 - \sum_{a \in \{-1,0,+1\}} p_a^2$$

Maximum $H = 1 - 3 \cdot (1/3)^2 = 2/3$ when uniform. $H = 0$ when fixed.

### Fixation Probability

Under neutral drift, the probability that allele $a$ with frequency $p_a$ eventually fixes:

$$P_{\text{fix}}(a) = p_a = \frac{k_a}{N}$$

This is the fundamental theorem of neutral evolution — fixation probability equals initial frequency.

### Effective Population Size

From entropy decay observations, we can estimate $N_e$:

$$N_e \approx \frac{t}{2 \ln(H_0 / H_t)}$$

where $H_0$ and $H_t$ are initial and final heterozygosity over $t$ generations.

### Entropy

Shannon entropy of the allele frequency distribution:

$$S = -\sum_a p_a \ln p_a$$

Maximum $S = \ln 3 \approx 1.099$ nats when uniform.

### Complexity

| Operation | Time |
|-----------|------|
| `drift_step(pop, rng)` | O(N) |
| `run_drift(pop, G, rng)` | O(G · N) |
| `heterozygosity(pop)` | O(N) |
| `entropy(pop)` | O(N) |
| `fixation_probability(k, N)` | O(1) |
| `effective_pop_size(H₀, Hₜ, t)` | O(1) |

Where N = population size, G = generations.

## Quick Start

```rust
use ternary_drift::{run_drift, heterozygosity, fixation_probability, effective_pop_size};

// Initial population: mixed ternary alleles
let mut pop = vec![0i8, 1, -1, 0, 1, -1, 0, 1, -1, 0];

// Deterministic RNG for reproducibility
let mut s: u64 = 42;
let mut rng = || {
    s = s.wrapping_mul(6364136223846793005).wrapping_add(1);
    (s >> 33) as f64 / (1u64 << 31) as f64
};

let initial_h = heterozygosity(&pop);
let result = run_drift(&mut pop, 100, &mut rng);

println!("Final distribution: {:?}", result.final_fracs);
println!("Diversity loss: {:.1}%", result.diversity_loss * 100.0);
println!("Fixation generation: {:?}", result.fixation_gen);

// Fixation probability (neutral theory)
let p_fix = fixation_probability(3, 10); // 3 copies out of 10
assert!((p_fix - 0.3).abs() < 1e-10);

// Effective population size from entropy decay
let ne = effective_pop_size(1.0, 0.5, 100.0);
```

## API

### Simulation

| Function | Description |
|----------|-------------|
| `drift_step(pop, rng)` | One Wright-Fisher generation |
| `run_drift(pop, generations, rng) → DriftResult` | Full simulation |

### Statistics

| Function | Description |
|----------|-------------|
| `heterozygosity(pop) → f64` | $1 - \sum p_a^2$ |
| `entropy(pop) → f64` | Shannon entropy (natural log) |
| `fixation_probability(count, size) → f64` | Neutral fixation probability |
| `effective_pop_size(H₀, Hₜ, t) → f64` | Estimated $N_e$ from entropy decay |

### DriftResult

| Field | Description |
|-------|-------------|
| `final_fracs: [f64; 3]` | Final allele frequencies (neg, zero, pos) |
| `fixation_gen: Option<usize>` | Generation where fixation occurred |
| `diversity_loss: f64` | Fraction of initial heterozygosity lost |

## Architecture Notes

The drift model demonstrates the **γ + η = C** conservation principle in population genetics:

- **γ (structure)**: the population size $N$ — the fixed carrying capacity that constrains drift dynamics
- **η (dynamics)**: stochastic sampling — the random parent selection that perturbs allele frequencies each generation
- **C (conservation)**: allele count conservation — $\sum_i X_i$ is a martingale under neutral drift (no selection), so the expected total allele count is conserved

The key insight: drift destroys diversity (η increases as populations diverge), but the expected frequency of each allele remains constant (C holds in expectation). The variance of allele frequencies grows as $\sim t/N$, quantifying when drift overwhelms deterministic forces.

## References

- Wright, S. (1931). *Evolution in Mendelian Populations*. Genetics. — Foundational paper on drift.
- Fisher, R.A. (1930). *The Genetical Theory of Natural Selection*. — Connection between drift and selection.
| Kimura, M. (1962). *On the Probability of Fixation of Mutant Genes in a Population*. Genetics.
| Crow, J.F. & Kimura, M. (1970). *An Introduction to Population Genetics Theory*. Harper & Row.
| Ewens, W.J. (2004). *Mathematical Population Genetics* (2nd ed.). Springer.

## License: MIT
