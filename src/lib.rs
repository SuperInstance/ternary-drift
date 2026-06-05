#![forbid(unsafe_code)]

/// Result of a genetic drift simulation.
#[derive(Debug, Clone)]
pub struct DriftResult {
    pub final_fracs: [f64; 3],
    pub fixation_gen: Option<usize>,
    pub diversity_loss: f64,
}

/// Perform one generation of Wright-Fisher drift on a ternary population.
/// Each individual randomly selects a parent from the previous generation.
pub fn drift_step(pop: &mut [i8], rng: &mut impl FnMut() -> f64) {
    let old = pop.to_vec();
    let n = old.len();
    for i in 0..n {
        let parent_idx = (rng() * n as f64) as usize;
        pop[i] = old[parent_idx.min(n - 1)];
    }
}

/// Run drift for a given number of generations.
pub fn run_drift(
    pop: &mut [i8],
    generations: usize,
    rng: &mut impl FnMut() -> f64,
) -> DriftResult {
    let initial_h = heterozygosity(pop);

    let mut fixation_gen = None;
    for gen in 0..generations {
        drift_step(pop, rng);
        let h = heterozygosity(pop);
        if h == 0.0 && fixation_gen.is_none() {
            fixation_gen = Some(gen + 1);
        }
    }

    let final_h = heterozygosity(pop);
    let diversity_loss = if initial_h > 0.0 {
        1.0 - final_h / initial_h
    } else {
        0.0
    };

    let mut counts = [0usize; 3];
    for &v in pop.iter() {
        let idx = (v + 1).max(0).min(2) as usize;
        counts[idx] += 1;
    }
    let n_pop = pop.len() as f64;
    let final_fracs = [
        counts[0] as f64 / n_pop,
        counts[1] as f64 / n_pop,
        counts[2] as f64 / n_pop,
    ];

    DriftResult {
        final_fracs,
        fixation_gen,
        diversity_loss,
    }
}

/// Compute effective population size from entropy change over generations.
pub fn effective_pop_size(initial_entropy: f64, final_entropy: f64, generations: f64) -> f64 {
    if final_entropy <= 0.0 || initial_entropy <= 0.0 || generations <= 0.0 {
        return f64::INFINITY;
    }
    let ratio = final_entropy / initial_entropy;
    if ratio >= 1.0 {
        return f64::INFINITY;
    }
    // Ne ≈ generations / (2 * ln(initial/final))
    generations / (2.0 * ratio.ln().abs())
}

/// Shannon-like heterozygosity for ternary alleles.
/// H = 1 - sum(p_i^2), which is 0 when fixed, maximal when uniform.
pub fn heterozygosity(pop: &[i8]) -> f64 {
    if pop.is_empty() {
        return 0.0;
    }
    let n = pop.len() as f64;
    let mut counts = [0usize; 3];
    for &v in pop.iter() {
        let idx = (v + 1).max(0).min(2) as usize;
        counts[idx] += 1;
    }
    let homo = counts.iter().map(|&c| (c as f64 / n).powi(2)).sum::<f64>();
    1.0 - homo
}

/// Fixation probability for an allele with `allele_count` copies in `pop_size`.
/// For neutral drift: k/N. For ternary, same principle applies.
pub fn fixation_probability(allele_count: usize, pop_size: usize) -> f64 {
    if pop_size == 0 {
        return 0.0;
    }
    allele_count as f64 / pop_size as f64
}

/// Helper: compute entropy of a population.
pub fn entropy(pop: &[i8]) -> f64 {
    if pop.is_empty() {
        return 0.0;
    }
    let n = pop.len() as f64;
    let mut counts = [0usize; 3];
    for &v in pop.iter() {
        let idx = (v + 1).max(0).min(2) as usize;
        counts[idx] += 1;
    }
    let mut h = 0.0;
    for &c in &counts {
        if c > 0 {
            let p = c as f64 / n;
            h -= p * p.ln();
        }
    }
    h
}

#[cfg(test)]
mod tests {
    use super::*;

    fn deterministic_rng(seed: u64) -> impl FnMut() -> f64 {
        let mut s = seed;
        move || {
            s = s.wrapping_mul(6364136223846793005).wrapping_add(1);
            (s >> 33) as f64 / (1u64 << 31) as f64
        }
    }

    #[test]
    fn test_drift_step_preserves_length() {
        let mut pop = vec![0i8, 1, -1, 0, 1];
        let mut rng = deterministic_rng(42);
        drift_step(&mut pop, &mut rng);
        assert_eq!(pop.len(), 5);
    }

    #[test]
    fn test_drift_step_values_valid() {
        let mut pop = vec![0i8, 1, -1, 0, 1];
        let mut rng = deterministic_rng(42);
        drift_step(&mut pop, &mut rng);
        for &v in &pop {
            assert!(v >= -1 && v <= 1);
        }
    }

    #[test]
    fn test_drift_step_deterministic() {
        let mut pop1 = vec![0i8, 1, -1, 0, 1, -1];
        let mut pop2 = vec![0i8, 1, -1, 0, 1, -1];
        let mut rng1 = deterministic_rng(123);
        let mut rng2 = deterministic_rng(123);
        drift_step(&mut pop1, &mut rng1);
        drift_step(&mut pop2, &mut rng2);
        assert_eq!(pop1, pop2);
    }

    #[test]
    fn test_heterozygosity_uniform() {
        let pop = vec![0i8, 1, -1];
        let h = heterozygosity(&pop);
        assert!((h - (1.0 - 3.0_f64 * (1.0_f64 / 3.0_f64).powi(2))).abs() < 1e-10);
    }

    #[test]
    fn test_heterozygosity_fixed() {
        let pop = vec![1i8, 1, 1, 1];
        assert_eq!(heterozygosity(&pop), 0.0);
    }

    #[test]
    fn test_heterozygosity_empty() {
        assert_eq!(heterozygosity(&[]), 0.0);
    }

    #[test]
    fn test_heterozygosity_range() {
        let pop = vec![0i8, 0, 1, -1];
        let h = heterozygosity(&pop);
        assert!(h >= 0.0 && h <= 1.0);
    }

    #[test]
    fn test_fixation_probability_basic() {
        assert!((fixation_probability(1, 4) - 0.25).abs() < 1e-10);
    }

    #[test]
    fn test_fixation_probability_zero() {
        assert_eq!(fixation_probability(0, 10), 0.0);
    }

    #[test]
    fn test_fixation_probability_full() {
        assert!((fixation_probability(10, 10) - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_fixation_probability_empty_pop() {
        assert_eq!(fixation_probability(5, 0), 0.0);
    }

    #[test]
    fn test_run_drift_fixation() {
        // All same allele — already fixed, heterozygosity stays 0 throughout
        let mut pop = vec![1i8; 10];
        let mut rng = deterministic_rng(42);
        let result = run_drift(&mut pop, 100, &mut rng);
        // Since initial_h == 0, diversity_loss == 0
        assert_eq!(result.diversity_loss, 0.0);
        // fixation_gen will be Some(1) since h==0 on gen 0
        assert!(result.fixation_gen.is_some());
    }

    #[test]
    fn test_run_drift_diversity_loss_range() {
        let mut pop = vec![0i8, 1, -1, 0, 1, -1, 0, 1];
        let mut rng = deterministic_rng(42);
        let result = run_drift(&mut pop, 50, &mut rng);
        assert!(result.diversity_loss >= 0.0 && result.diversity_loss <= 1.0);
    }

    #[test]
    fn test_run_drift_final_fracs_sum_to_one() {
        let mut pop = vec![0i8, 1, -1, 0, 1, -1];
        let mut rng = deterministic_rng(99);
        let result = run_drift(&mut pop, 20, &mut rng);
        let sum: f64 = result.final_fracs.iter().sum();
        assert!((sum - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_effective_pop_size_basic() {
        let ne = effective_pop_size(1.0, 0.5, 100.0);
        assert!(ne > 0.0 && ne.is_finite());
    }

    #[test]
    fn test_effective_pop_size_no_change() {
        let ne = effective_pop_size(1.0, 1.0, 100.0);
        assert!(ne.is_infinite());
    }

    #[test]
    fn test_entropy_uniform() {
        let pop = vec![0i8, 1, -1];
        let e = entropy(&pop);
        // ln(3) ≈ 1.0986
        assert!((e - 3.0_f64.ln()).abs() < 1e-10);
    }
}
