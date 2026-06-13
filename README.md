# Ternary Drift — Genetic Drift Simulation for Ternary Populations

**Ternary Drift** simulates Wright-Fisher genetic drift on populations where each individual carries one of three alleles {-1, 0, +1}. It provides single-generation drift steps, multi-generation simulation with fixation detection, heterozygosity tracking, effective population size estimation, and Shannon entropy analysis — all using the ternary allele model where the neutral allele (0) represents the critical "spindle" state.

## Why It Matters

Genetic drift — random fluctuations in allele frequency — is the dominant evolutionary force in small populations. In ternary agent systems, drift explains how fleet diversity collapses over time even without selection pressure: neutral agents (state 0) can fixate through pure chance, causing the system to lose the {-1, +1} diversity that enables adaptation. Understanding drift rates is critical for setting fleet sizes: too few agents, and diversity is lost in O(N) generations; too many, and coordination becomes impossible. This crate provides the quantitative tools to measure, predict, and mitigate drift-induced diversity loss.

## How It Works

### Wright-Fisher Model

Each generation, every individual in the population independently samples a parent uniformly at random from the previous generation. This is a Moran/Wright-Fisher process:

```
new_pop[i] = old_pop[random(0..N)]
```

for each i in 0..N. This is O(N) per generation. The key property: allele frequencies follow a random walk bounded by [0, N], and absorption (fixation) is guaranteed eventually.

### Heterozygosity

Heterozygosity measures genetic diversity:

```
H = 1 - Σ(pᵢ²)   where pᵢ = frequency of allele i
```

H = 0 when one allele is fixed; H = 2/3 when all three alleles are equally frequent (maximum diversity). Heterozygosity decreases monotonically under pure drift, with expected decrease:

```
E[H(t+1)] = (1 - 1/N) · H(t)
```

### Effective Population Size

From the entropy decay rate:

```
Ne ≈ t / (2 · |ln(H_final / H_initial)|)
```

When Ne ≪ N (census size), it indicates that population structure (variance in reproductive success) reduces the effective gene pool.

### Fixation

Fixation occurs when all individuals carry the same allele. The probability of ultimate fixation of any allele equals its initial frequency. Expected fixation time is O(N) generations for a neutral allele.

## Quick Start

```rust
use ternary_drift::{drift_step, run_drift, heterozygosity};

// Start with equal frequencies of {-1, 0, +1}
let mut pop: Vec<i8> = (0..300).map(|i| match i % 3 { 0 => -1, 1 => 0, _ => 1 }).collect();

let h0 = heterozygosity(&pop);
let result = run_drift(&mut pop, 100, &mut rand::random);
let h_final = heterozygosity(&pop);

println!("Diversity loss: {:.1}%", (1.0 - h_final / h0) * 100.0);
```

```bash
cargo add ternary-drift
```

## API

| Type / Function | Description |
|---|---|
| `drift_step(&mut [i8], rng)` | One Wright-Fisher generation (O(N)) |
| `run_drift(&mut [i8], gens, rng)` | Multi-gen simulation → `DriftResult` |
| `heterozygosity(&[i8]) → f64` | Diversity measure (0 = fixed, ⅔ = uniform) |
| `effective_pop_size(h₀, h_f, gens)` | Ne estimation from entropy decay |
| `DriftResult` | `{ final_fracs, fixation_gen, diversity_loss }` |

## Architecture Notes

Drift models how **SuperInstance** fleet diversity degrades over time. Without active intervention, the η (entropy) term in γ + η = C decreases as drift fixates agents into a single state. The fleet must inject controlled diversity (γ pulses) to maintain the conservation balance. See [Architecture](https://github.com/SuperInstance/SuperInstance/blob/main/ARCHITECTURE.md).

## References

- Wright, Sewall. "Evolution in Mendelian Populations," *Genetics*, 16(2), 1931 — Wright-Fisher model.
- Crow, James F. & Kimura, Motoo. *An Introduction to Population Genetics Theory*, Harper & Row, 1970.
- Ewens, Warren J. *Mathematical Population Genetics*, 2nd ed., Springer, 2004.

## License

MIT
