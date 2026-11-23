//! Sequence statistics and discrete distributions.
//!
//! Poisson is implemented on `libm` (log-gamma and the regularised incomplete gamma
//! function) and works in `no_std`. Binomial delegates to `tpt-math` and needs the
//! `std` feature. Invalid parameters yield `NaN`.

use core::fmt;

use crate::nucleotide::Base;

pub mod cluster;
pub mod linalg;
pub mod normalize;
pub mod pca;
pub mod phenotype;
pub mod qc;
pub(crate) mod rng;
pub mod tsne;
pub mod umap;

pub use linalg::Matrix;

/// Errors from the numeric routines in this module.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum StatsError {
    /// The input has no rows, no columns, or no values.
    Empty,
    /// Dimensions disagree or a shape requirement is not met.
    Shape(&'static str),
    /// An input value is NaN or infinite.
    NonFinite,
    /// A parameter is out of range, such as `k = 0` for clustering.
    InvalidArgument(&'static str),
    /// The statistic is undefined for this data, such as zero variance.
    Degenerate(&'static str),
}

impl fmt::Display for StatsError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Empty => f.write_str("input is empty"),
            Self::Shape(msg) => write!(f, "shape error: {msg}"),
            Self::NonFinite => f.write_str("input contains NaN or infinity"),
            Self::InvalidArgument(msg) => write!(f, "invalid argument: {msg}"),
            Self::Degenerate(msg) => write!(f, "statistic undefined: {msg}"),
        }
    }
}

impl core::error::Error for StatsError {}

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

/// Relative tolerance at which the incomplete gamma iterations stop.
const GAMMA_EPS: f64 = 1e-15;
/// Smallest magnitude used as a stand-in for zero in the continued fraction.
const GAMMA_FPMIN: f64 = 1e-300;
/// Iteration cap for the incomplete gamma series and continued fraction.
const GAMMA_MAX_ITER: usize = 200_000;

/// Regularised upper incomplete gamma `Q(a, x) = Γ(a, x) / Γ(a)` for `a > 0`, `x >= 0`.
///
/// Uses the power series for `x < a + 1` and Lentz's continued fraction otherwise,
/// as in Numerical Recipes, so that both tails keep full relative accuracy.
#[allow(clippy::many_single_char_names)]
fn gamma_q(a: f64, x: f64) -> f64 {
    if x <= 0.0 {
        return 1.0;
    }
    // Common prefactor exp(-x + a ln x - ln Γ(a)).
    let front = libm::exp(-x + a * libm::log(x) - libm::lgamma(a));
    if x < a + 1.0 {
        // P(a, x) by series, then Q = 1 - P.
        let mut ap = a;
        let mut term = 1.0 / a;
        let mut sum = term;
        for _ in 0..GAMMA_MAX_ITER {
            ap += 1.0;
            term *= x / ap;
            sum += term;
            if libm::fabs(term) < libm::fabs(sum) * GAMMA_EPS {
                break;
            }
        }
        1.0 - front * sum
    } else {
        // Q(a, x) by modified Lentz continued fraction.
        let mut b = x + 1.0 - a;
        let mut c = 1.0 / GAMMA_FPMIN;
        let mut d = 1.0 / b;
        let mut h = d;
        for i in 1..GAMMA_MAX_ITER {
            let i = i as f64;
            let an = -i * (i - a);
            b += 2.0;
            d = an * d + b;
            if libm::fabs(d) < GAMMA_FPMIN {
                d = GAMMA_FPMIN;
            }
            c = b + an / c;
            if libm::fabs(c) < GAMMA_FPMIN {
                c = GAMMA_FPMIN;
            }
            d = 1.0 / d;
            let delta = d * c;
            h *= delta;
            if libm::fabs(delta - 1.0) < GAMMA_EPS {
                break;
            }
        }
        front * h
    }
}

/// Poisson probability mass `P(X = k)` for rate `lambda`.
///
/// Computed in log space as `exp(k ln λ - λ - ln k!)`, with `ln k!` from `lgamma`.
/// Returns `NaN` if `lambda` is negative, NaN or infinite.
#[must_use]
pub fn poisson_pmf(lambda: f64, k: u64) -> f64 {
    if !lambda.is_finite() || lambda < 0.0 {
        return f64::NAN;
    }
    if lambda <= 0.0 {
        return if k == 0 { 1.0 } else { 0.0 };
    }
    let kf = k as f64;
    libm::exp(kf * libm::log(lambda) - lambda - libm::lgamma(kf + 1.0))
}

/// Poisson cumulative probability `P(X <= k)` for rate `lambda`.
///
/// Uses the identity `P(X <= k) = Q(k + 1, λ)` with the regularised upper incomplete
/// gamma function. Returns `NaN` if `lambda` is negative, NaN or infinite.
#[must_use]
pub fn poisson_cdf(lambda: f64, k: u64) -> f64 {
    if !lambda.is_finite() || lambda < 0.0 {
        return f64::NAN;
    }
    if lambda <= 0.0 {
        return 1.0;
    }
    gamma_q(k as f64 + 1.0, lambda)
}

/// Poisson quantile: the smallest `k` with `P(X <= k) >= p`.
///
/// Returns `NaN` if `lambda` is invalid or `p` is outside `[0, 1]`. `p = 0` gives `0`
/// and `p = 1` gives infinity. The result is returned as `f64` because the quantile
/// of a Poisson distribution is unbounded.
#[must_use]
pub fn poisson_quantile(lambda: f64, p: f64) -> f64 {
    if !lambda.is_finite() || lambda < 0.0 || !(0.0..=1.0).contains(&p) {
        return f64::NAN;
    }
    if p <= 0.0 {
        return 0.0;
    }
    if p >= 1.0 {
        return f64::INFINITY;
    }
    // Double the upper bound until it covers p, then bisect over integers.
    let mut hi: u64 = 1;
    while poisson_cdf(lambda, hi) < p && hi < (1 << 62) {
        hi *= 2;
    }
    let mut lo: u64 = 0;
    while lo < hi {
        let mid = lo + (hi - lo) / 2;
        if poisson_cdf(lambda, mid) >= p {
            hi = mid;
        } else {
            lo = mid + 1;
        }
    }
    lo as f64
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

    #[test]
    fn poisson_pmf_hand_computed() {
        // lambda = 2: P(0) = e^-2, P(1) = P(2) = 2e^-2, P(3) = (8/6)e^-2.
        let e2 = 0.135_335_283_236_612_7_f64;
        assert!((poisson_pmf(2.0, 0) - e2).abs() < 1e-12);
        assert!((poisson_pmf(2.0, 1) - 0.270_670_566_473_225_4).abs() < 1e-12);
        assert!((poisson_pmf(2.0, 2) - 0.270_670_566_473_225_4).abs() < 1e-12);
        assert!((poisson_pmf(2.0, 3) - 0.180_447_044_315_483_6).abs() < 1e-12);
        assert!((poisson_pmf(0.0, 0) - 1.0).abs() < 1e-15);
        assert!(poisson_pmf(0.0, 4).abs() < 1e-15);
        assert!(poisson_pmf(-1.0, 0).is_nan());
        assert!(poisson_pmf(f64::NAN, 0).is_nan());
    }

    #[test]
    fn poisson_cdf_hand_computed() {
        assert!((poisson_cdf(2.0, 0) - 0.135_335_283_236_612_7).abs() < 1e-12);
        assert!((poisson_cdf(2.0, 1) - 0.406_005_849_709_838_1).abs() < 1e-12);
        assert!((poisson_cdf(2.0, 2) - 0.676_676_416_183_063_5).abs() < 1e-12);
        assert!((poisson_cdf(2.0, 3) - 0.857_123_460_498_547_1).abs() < 1e-12);
        assert!((poisson_cdf(0.0, 0) - 1.0).abs() < 1e-15);
        assert!(poisson_cdf(-0.5, 1).is_nan());
    }

    #[test]
    fn poisson_cdf_equals_sum_of_pmf_terms() {
        for &lambda in &[0.3, 2.0, 5.0, 17.5, 60.0] {
            let mut running = 0.0;
            for k in 0..=120u64 {
                running += poisson_pmf(lambda, k);
                let cdf = poisson_cdf(lambda, k);
                assert!(
                    (cdf - running).abs() < 1e-11,
                    "lambda={lambda} k={k} cdf={cdf} sum={running}"
                );
            }
        }
    }

    #[test]
    fn poisson_quantile_round_trip() {
        for &lambda in &[0.5, 3.0, 10.0, 50.0, 400.0] {
            for &p in &[0.001, 0.1, 0.5, 0.9, 0.999] {
                let k = poisson_quantile(lambda, p) as u64;
                assert!(poisson_cdf(lambda, k) >= p, "lambda={lambda} p={p} k={k}");
                if k > 0 {
                    assert!(poisson_cdf(lambda, k - 1) < p, "lambda={lambda} p={p} k={k}");
                }
            }
        }
        assert!(poisson_quantile(2.0, 0.0).abs() < 1e-15);
        assert!(poisson_quantile(2.0, 1.0).is_infinite());
        assert!(poisson_quantile(2.0, 1.5).is_nan());
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
