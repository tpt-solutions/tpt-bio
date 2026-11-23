//! Column-wise feature normalization, optionally within groups such as plates.
//!
//! Rows are samples (wells, cells, images) and columns are features. Each column
//! is normalized independently. With grouping, the statistics are computed
//! separately inside each group, so a plate effect is removed without mixing
//! plates together.
//!
//! A column with zero spread (standard deviation or MAD of zero) maps to all zeros
//! rather than dividing by zero.

use alloc::collections::BTreeMap;
use alloc::vec::Vec;

use super::{Matrix, StatsError};

/// The centering and scaling rule.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Method {
    /// `(x - mean) / sd`, with the population standard deviation (divide by n).
    ZScore,
    /// `(x - median) / (1.4826 * MAD)`. The factor 1.4826 makes the MAD consistent
    /// with the standard deviation for normally distributed data.
    Robust,
}

/// Scale factor that makes the median absolute deviation estimate sigma under normality.
pub const MAD_SCALE: f64 = 1.4826;

/// Normalizes every column of `data` over all rows.
///
/// # Errors
/// Returns an error only if the matrix is empty, which [`Matrix`] already prevents.
pub fn normalize_columns(data: &Matrix, method: Method) -> Result<Matrix, StatsError> {
    let rows: Vec<usize> = (0..data.rows()).collect();
    let mut out = data.clone();
    normalize_subset(data, &rows, method, &mut out);
    Ok(out)
}

/// Normalizes every column separately inside each group of rows.
///
/// `groups[i]` is the label of row `i`. Labels can be any ordered type, such as a
/// plate name or a plate number.
///
/// # Errors
/// Returns [`StatsError::Shape`] when `groups.len()` differs from the row count.
pub fn normalize_by_group<K: Ord>(
    data: &Matrix,
    groups: &[K],
    method: Method,
) -> Result<Matrix, StatsError> {
    if groups.len() != data.rows() {
        return Err(StatsError::Shape("one group label is required per row"));
    }
    let mut members: BTreeMap<&K, Vec<usize>> = BTreeMap::new();
    for (i, g) in groups.iter().enumerate() {
        members.entry(g).or_default().push(i);
    }
    let mut out = data.clone();
    for rows in members.values() {
        normalize_subset(data, rows, method, &mut out);
    }
    Ok(out)
}

/// Writes the normalized values of `rows` into `out`, using statistics from those rows only.
fn normalize_subset(data: &Matrix, rows: &[usize], method: Method, out: &mut Matrix) {
    for c in 0..data.cols() {
        let values: Vec<f64> = rows.iter().map(|&r| data.get(r, c)).collect();
        let (center, scale) = match method {
            Method::ZScore => {
                let n = values.len() as f64;
                let mean = values.iter().sum::<f64>() / n;
                let var = values.iter().map(|v| (v - mean) * (v - mean)).sum::<f64>() / n;
                (mean, libm::sqrt(var))
            }
            Method::Robust => {
                let med = median(&values);
                let deviations: Vec<f64> = values.iter().map(|v| (v - med).abs()).collect();
                (med, MAD_SCALE * median(&deviations))
            }
        };
        for (&r, &v) in rows.iter().zip(&values) {
            let z = if scale == 0.0 { 0.0 } else { (v - center) / scale };
            out.set(r, c, z);
        }
    }
}

/// Median of a non-empty slice. Averages the two middle values for even lengths.
///
/// # Panics
/// Panics if `values` is empty.
pub(crate) fn median(values: &[f64]) -> f64 {
    let mut sorted = values.to_vec();
    sorted.sort_by(f64::total_cmp);
    let n = sorted.len();
    if n % 2 == 1 { sorted[n / 2] } else { f64::midpoint(sorted[n / 2 - 1], sorted[n / 2]) }
}

#[cfg(test)]
mod tests {
    use super::*;
    use alloc::vec;

    fn approx(a: f64, b: f64) -> bool {
        (a - b).abs() < 1e-12
    }

    #[test]
    fn zscore_of_known_column() {
        // Mean 5, population variance (16+4+0+4+16)/5 = 8.
        let m = Matrix::new(5, 1, vec![1.0, 3.0, 5.0, 7.0, 9.0]).unwrap();
        let z = normalize_columns(&m, Method::ZScore).unwrap();
        let sd = libm::sqrt(8.0);
        assert!(approx(z.get(0, 0), -4.0 / sd));
        assert!(approx(z.get(2, 0), 0.0));
        assert!(approx(z.get(4, 0), 4.0 / sd));
        let mean: f64 = z.column(0).iter().sum::<f64>() / 5.0;
        assert!(approx(mean, 0.0));
    }

    #[test]
    fn robust_uses_median_and_mad() {
        // Median 3, deviations [2,1,0,1,97] -> MAD 1. The outlier does not move the centre.
        let m = Matrix::new(5, 1, vec![1.0, 2.0, 3.0, 4.0, 100.0]).unwrap();
        let z = normalize_columns(&m, Method::Robust).unwrap();
        assert!(approx(z.get(2, 0), 0.0));
        assert!(approx(z.get(0, 0), -2.0 / MAD_SCALE));
        assert!(approx(z.get(4, 0), 97.0 / MAD_SCALE));
    }

    #[test]
    fn constant_column_maps_to_zero() {
        let m = Matrix::new(3, 2, vec![1.0, 7.0, 1.0, 8.0, 1.0, 9.0]).unwrap();
        let z = normalize_columns(&m, Method::ZScore).unwrap();
        assert_eq!(z.column(0), vec![0.0, 0.0, 0.0]);
        assert!(approx(z.get(0, 1), -1.0 / libm::sqrt(2.0 / 3.0)));
    }

    #[test]
    fn per_plate_removes_plate_offset() {
        // Plate B is plate A shifted by +100. Per-plate z-scores must agree.
        let m = Matrix::new(4, 1, vec![1.0, 3.0, 101.0, 103.0]).unwrap();
        let plates = ["A", "A", "B", "B"];
        let z = normalize_by_group(&m, &plates, Method::ZScore).unwrap();
        assert!(approx(z.get(0, 0), z.get(2, 0)));
        assert!(approx(z.get(1, 0), z.get(3, 0)));
        assert!(approx(z.get(0, 0), -1.0));
    }

    #[test]
    fn group_label_count_must_match() {
        let m = Matrix::new(2, 1, vec![1.0, 2.0]).unwrap();
        assert!(matches!(
            normalize_by_group(&m, &[1u8], Method::ZScore),
            Err(StatsError::Shape(_))
        ));
    }
}
