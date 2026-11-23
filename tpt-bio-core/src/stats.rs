//! Sequence statistics and discrete distributions.
//!
//! Discrete distributions (Poisson, Binomial) delegate to `tpt-math` and need the
//! `std` feature. Invalid parameters yield `NaN`, matching the `tpt-math` convention.

use crate::nucleotide::Base;

/// Counts of A, C, G and T (case-insensitive). Other bytes are ignored.
#[must_use]
pub fn base_counts(seq: &[u8]) -> [u64; 4] {
    let mut counts = [0u64; 4];
    for &b in seq {
        if let Some(base) = Base::from_ascii(b) {
            counts[usize::from(base.to_2bit())] += 1;
        }
    }
    counts
}

/// Fraction of unambiguous bases that are G or C.
///
/// Returns `None` when the sequence has no A, C, G or T.
#[must_use]
pub fn gc_content(seq: &[u8]) -> Option<f64> {
    let [a, c, g, t] = base_counts(seq);
    let total = a + c + g + t;
    if total == 0 { None } else { Some((c + g) as f64 / total as f64) }
}

/// Shannon entropy in bits of a discrete distribution given by counts.
///
/// Returns `None` when all counts are zero. Zero counts contribute nothing.
#[must_use]
pub fn shannon_entropy_bits(counts: &[u64]) -> Option<f64> {
    let total: u64 = counts.iter().sum();
    if total == 0 {
        return None;
    }
    let n = total as f64;
    let h = counts
        .iter()
        .filter(|&&c| c > 0)
        .map(|&c| {
            let p = c as f64 / n;
            -p * libm::log2(p)
        })
        .sum::<f64>();
    Some(h)
}

/// Shannon entropy in bits of the base composition of a sequence.
///
/// Returns `None` when the sequence has no A, C, G or T.
#[must_use]
pub fn sequence_entropy_bits(seq: &[u8]) -> Option<f64> {
    shannon_entropy_bits(&base_counts(seq))
}

/// Poisson probability mass `P(X = k)` for rate `lambda`.
#[cfg(feature = "std")]
#[must_use]
pub fn poisson_pmf(lambda: f64, k: u64) -> f64 {
    tpt_math_prob_dist::density::poisson_pmf(lambda, k)
}

/// Poisson cumulative probability `P(X <= k)` for rate `lambda`.
#[cfg(feature = "std")]
#[must_use]
pub fn poisson_cdf(lambda: f64, k: u64) -> f64 {
    tpt_math_prob_dist::density::poisson_cdf(lambda, k)
}

/// Poisson quantile: the smallest `k` with `P(X <= k) >= p`.
#[cfg(feature = "std")]
#[must_use]
pub fn poisson_quantile(lambda: f64, p: f64) -> f64 {
    tpt_math_prob_dist::density::poisson_quantile(lambda, p)
}

/// Binomial probability mass `P(X = k)` for `n` trials with success probability `p`.
#[cfg(feature = "std")]
#[must_use]
pub fn binomial_pmf(n: u64, p: f64, k: u64) -> f64 {
    tpt_math_prob_dist::density::binomial_pmf(n, p, k)
}

/// Binomial cumulative probability `P(X <= k)`.
#[cfg(feature = "std")]
#[must_use]
pub fn binomial_cdf(n: u64, p: f64, k: u64) -> f64 {
    tpt_math_prob_dist::density::binomial_cdf(n, p, k)
}

/// Binomial quantile: the smallest `k` in `0..=n` with `P(X <= k) >= p`.
///
/// Returns `None` if `p` is outside `[0, 1]` or is NaN.
#[cfg(feature = "std")]
#[must_use]
pub fn binomial_quantile(n: u64, p: f64, q: f64) -> Option<u64> {
    if !(0.0..=1.0).contains(&q) || !(0.0..=1.0).contains(&p) {
        return None;
    }
    // Linear scan of the CDF. Binomial quantiles are only needed for small n in
    // practice, so the simple search is clearer than a bisection.
    (0..=n).find(|&k| binomial_cdf(n, p, k) >= q).or(Some(n))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn gc_and_counts() {
        assert_eq!(base_counts(b"AcGTnN"), [1, 1, 1, 1]);
        assert_eq!(gc_content(b"GGCC"), Some(1.0));
        assert_eq!(gc_content(b"ATAT"), Some(0.0));
        assert_eq!(gc_content(b"ACGT"), Some(0.5));
        assert_eq!(gc_content(b"NNNN"), None);
    }

    #[test]
    fn entropy_known_values() {
        assert!((shannon_entropy_bits(&[1, 1, 1, 1]).unwrap() - 2.0).abs() < 1e-12);
        assert!(shannon_entropy_bits(&[5, 0, 0, 0]).unwrap().abs() < 1e-12);
        assert_eq!(shannon_entropy_bits(&[0, 0]), None);
    }

    #[cfg(feature = "std")]
    #[test]
    fn poisson_matches_closed_form() {
        // P(X = 0) = e^-2
        assert!((poisson_pmf(2.0, 0) - libm::exp(-2.0)).abs() < 1e-12);
        assert!((poisson_cdf(2.0, 0) - libm::exp(-2.0)).abs() < 1e-12);
    }

    #[cfg(feature = "std")]
    #[test]
    fn binomial_quantile_is_consistent_with_cdf() {
        let (n, p) = (20, 0.3);
        let k = binomial_quantile(n, p, 0.5).unwrap();
        assert!(binomial_cdf(n, p, k) >= 0.5);
        if k > 0 {
            assert!(binomial_cdf(n, p, k - 1) < 0.5);
        }
        assert_eq!(binomial_quantile(n, p, 1.5), None);
    }
}
