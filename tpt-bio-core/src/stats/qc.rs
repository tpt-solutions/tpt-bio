//! Assay quality metrics computed from positive and negative control wells.
//!
//! Both metrics use the sample standard deviation (divide by `n - 1`).
//!
//! - Z-prime (Zhang, Chung and Oldenburg, 1999):
//!   `Z' = 1 - 3 (sd_p + sd_n) / |mean_p - mean_n|`. Values above 0.5 are usually
//!   read as an excellent assay.
//! - SSMD (strictly standardized mean difference):
//!   `(mean_p - mean_n) / sqrt(sd_p^2 + sd_n^2)`.

use super::StatsError;

/// Mean and sample standard deviation (`n - 1` denominator).
///
/// # Errors
/// Returns [`StatsError::Degenerate`] with fewer than two values, and
/// [`StatsError::NonFinite`] for NaN or infinite input.
pub fn mean_and_sd(values: &[f64]) -> Result<(f64, f64), StatsError> {
    if values.len() < 2 {
        return Err(StatsError::Degenerate("need at least two values per control"));
    }
    if values.iter().any(|v| !v.is_finite()) {
        return Err(StatsError::NonFinite);
    }
    let n = values.len() as f64;
    let mean = values.iter().sum::<f64>() / n;
    let ss: f64 = values.iter().map(|v| (v - mean) * (v - mean)).sum();
    Ok((mean, libm::sqrt(ss / (n - 1.0))))
}

/// Z-prime between positive and negative controls.
///
/// # Errors
/// Returns [`StatsError::Degenerate`] when the two means are equal, since the
/// separation term then divides by zero. Also propagates [`mean_and_sd`] errors.
pub fn z_prime(positive: &[f64], negative: &[f64]) -> Result<f64, StatsError> {
    let (mp, sp) = mean_and_sd(positive)?;
    let (mn, sn) = mean_and_sd(negative)?;
    let separation = (mp - mn).abs();
    if separation == 0.0 {
        return Err(StatsError::Degenerate("control means are equal"));
    }
    Ok(1.0 - 3.0 * (sp + sn) / separation)
}

/// Strictly standardized mean difference, `mean_p - mean_n` over the pooled spread.
///
/// # Errors
/// Returns [`StatsError::Degenerate`] when both controls have zero spread, and
/// propagates [`mean_and_sd`] errors.
pub fn ssmd(positive: &[f64], negative: &[f64]) -> Result<f64, StatsError> {
    let (mp, sp) = mean_and_sd(positive)?;
    let (mn, sn) = mean_and_sd(negative)?;
    let spread = libm::sqrt(sp * sp + sn * sn);
    if spread == 0.0 {
        return Err(StatsError::Degenerate("both controls have zero spread"));
    }
    Ok((mp - mn) / spread)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn sample_sd_by_hand() {
        // [2, 4, 4, 4, 5, 5, 7, 9]: mean 5, sum of squares 32, sample variance 32/7.
        let (m, s) = mean_and_sd(&[2.0, 4.0, 4.0, 4.0, 5.0, 5.0, 7.0, 9.0]).unwrap();
        assert!((m - 5.0).abs() < 1e-12);
        assert!((s - libm::sqrt(32.0 / 7.0)).abs() < 1e-12);
    }

    #[test]
    fn z_prime_by_hand() {
        // Positive: mean 12, sd 2. Negative: mean 2, sd 2.
        // Z' = 1 - 3 (2 + 2) / |12 - 2| = 1 - 12/10 = -0.2.
        let z = z_prime(&[10.0, 12.0, 14.0], &[0.0, 2.0, 4.0]).unwrap();
        assert!((z - (-0.2)).abs() < 1e-12, "{z}");
    }

    #[test]
    fn ssmd_by_hand() {
        // (12 - 2) / sqrt(2^2 + 2^2) = 10 / sqrt(8).
        let s = ssmd(&[10.0, 12.0, 14.0], &[0.0, 2.0, 4.0]).unwrap();
        assert!((s - 10.0 / libm::sqrt(8.0)).abs() < 1e-12, "{s}");
    }

    #[test]
    fn zero_spread_controls_give_perfect_z_prime() {
        // sd 0 for both: Z' = 1 - 0 / 10 = 1. SSMD is undefined.
        assert!((z_prime(&[10.0, 10.0], &[0.0, 0.0]).unwrap() - 1.0).abs() < 1e-12);
        assert!(matches!(ssmd(&[10.0, 10.0], &[0.0, 0.0]), Err(StatsError::Degenerate(_))));
    }

    #[test]
    fn invalid_controls_are_rejected() {
        assert!(matches!(z_prime(&[1.0], &[0.0, 1.0]), Err(StatsError::Degenerate(_))));
        // Both means are 2, so the separation is zero.
        assert!(matches!(z_prime(&[1.0, 3.0], &[0.0, 4.0]), Err(StatsError::Degenerate(_))));
        assert_eq!(ssmd(&[1.0, f64::NAN], &[0.0, 1.0]), Err(StatsError::NonFinite));
    }
}
