# ternary-drift

**Genetic drift for ternary populations. The math of how diversity dies.**

Genetic drift is the slow, random loss of diversity in a finite population. It's not selection — there's no fitness difference. It's just chance. In a population of 10, one allele might vanish in a few generations. In a population of 10,000, it takes thousands. The math is the same either way: random sampling, generation after generation, until one variant *fixes* and the others are gone.

This crate implements Wright-Fisher drift for ternary populations where each individual carries one of three alleles: `{-1, 0, +1}`. You provide a random number generator, we provide the dynamics.

## What's Inside

- **`drift_step(pop, rng)`** — one generation of Wright-Fisher sampling
- **`run_drift(pop, generations, rng)`** — full simulation, returns `DriftResult`
- **`DriftResult`** — final allele fractions, fixation generation, diversity loss
- **`heterozygosity(pop)`** — Shannon diversity: `H = 1 - Σ(pᵢ²)`. Zero when fixed, maximal when uniform
- **`effective_pop_size()`** — estimate Ne from entropy decay rate
- **`fixation_probability(pop, allele)`** — probability that a specific allele takes over
- **`time_to_fixation(pop, rng)`** — expected generations until one allele dominates

## Quick Example

```rust
use ternary_drift::*;

// Population of 100 ternary agents: equal thirds
let mut pop: Vec<i8> = vec![-1; 33].into_iter()
    .chain(vec![0; 34])
    .chain(vec![1; 33])
    .collect();

let mut rng = || { /* your RNG */ 0.5 };

// Run 500 generations of drift
let result = run_drift(&mut pop, 500, &mut rng);

println!("Diversity loss: {:.1}%", result.diversity_loss * 100.0);
println!("Fixed at generation: {:?}", result.fixation_gen);
println!("Final fractions: -1={:.2}, 0={:.2}, +1={:.2}",
    result.final_fracs[0], result.final_fracs[1], result.final_fracs[2]);

// Heterozygosity: how diverse is the population?
let h = heterozygosity(&pop); // 0.0 = monoculture, ~0.667 = maximum diversity
```

## The Insight

**Drift is inevitable in finite populations.** You can't stop it — you can only slow it down with larger populations. In ternary systems, the 0 state has a peculiar property: once the population drifts to 0-fixation, it's a *monoculture* — no variation, no adaptability, no resilience. This is why tunneling (random reactivation) matters in ternary agent systems: it's the only force that counteracts drift.

**Use cases:**
- **Population genetics** — model allele frequency dynamics in finite populations
- **Evolutionary computation** — understand diversity loss in genetic algorithms
- **Multi-agent systems** — analyze how agent diversity degrades over time
- **Conservation biology** — estimate effective population size from genetic data
- **Cultural evolution** — how ideas/habits fix or drift in social populations

## See Also

- **ternary-experiment** — parameter sweeps for drift experiments at scale
- **ternary-percolation** — spatial drift on grids
- **ternary-minority** — the minority rule that fights fixation
- **ternary-life** — lifecycle dynamics beyond simple drift

## Install

```bash
cargo add ternary-drift
```

## License

MIT
