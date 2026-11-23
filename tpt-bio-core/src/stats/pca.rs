//! Principal component analysis by eigendecomposition of the sample covariance.

use alloc::vec;
use alloc::vec::Vec;

use super::Matrix;
use super::StatsError;
use super::linalg::symmetric_eigen;

/// A fitted PCA model with the top `k` components.
#[derive(Debug, Clone, PartialEq)]
pub struct Pca {
    mean: Vec<f64>,
    /// `k x features`. Row `i` is the unit-length principal axis `i`.
    components: Matrix,
    explained_variance: Vec<f64>,
    explained_variance_ratio: Vec<f64>,
}

impl Pca {
    /// Fits PCA on `data` (rows are samples, columns are features), keeping `n_components`.
    ///
    /// The covariance uses the `n - 1` denominator. Each component is normalized to
    /// unit length, and its sign is chosen so that its largest-magnitude entry is
    /// positive, which makes the output deterministic.
    ///
    /// The explained variance ratio is the eigenvalue divided by the total variance
    /// (the trace of the covariance), so it is relative to all features, not only
    /// the retained components. If the total variance is zero, every ratio is zero.
    ///
    /// # Errors
    /// Returns [`StatsError::InvalidArgument`] when `n_components` is zero or exceeds
    /// the feature count, and [`StatsError::Shape`] when there are fewer than two rows.
    pub fn fit(data: &Matrix, n_components: usize) -> Result<Self, StatsError> {
        let (n, d) = (data.rows(), data.cols());
        if n < 2 {
            return Err(StatsError::Shape("PCA needs at least two rows"));
        }
        if n_components == 0 || n_components > d {
            return Err(StatsError::InvalidArgument(
                "n_components must be between 1 and the feature count",
            ));
        }
        let mean: Vec<f64> =
            (0..d).map(|c| data.column(c).iter().sum::<f64>() / n as f64).collect();

        let mut cov = Matrix::zeros(d, d);
        for r in 0..n {
            let row = data.row(r);
            for i in 0..d {
                let xi = row[i] - mean[i];
                for j in 0..d {
                    let v = cov.get(i, j) + xi * (row[j] - mean[j]);
                    cov.set(i, j, v);
                }
            }
        }
        let denom = (n - 1) as f64;
        for v in cov.as_mut_slice() {
            *v /= denom;
        }

        let total: f64 = (0..d).map(|i| cov.get(i, i)).sum();
        let eig = symmetric_eigen(&cov)?;

        let mut components = Matrix::zeros(n_components, d);
        for k in 0..n_components {
            let mut axis = eig.vectors.column(k);
            let pivot = axis
                .iter()
                .copied()
                .fold(0.0_f64, |best, x| if x.abs() > best.abs() { x } else { best });
            if pivot < 0.0 {
                for x in &mut axis {
                    *x = -*x;
                }
            }
            for (j, &x) in axis.iter().enumerate() {
                components.set(k, j, x);
            }
        }

        let explained_variance: Vec<f64> = eig.values[..n_components].to_vec();
        let explained_variance_ratio =
            explained_variance.iter().map(|&v| if total > 0.0 { v / total } else { 0.0 }).collect();
        Ok(Self { mean, components, explained_variance, explained_variance_ratio })
    }

    /// Projects `data` onto the fitted components. The result is `rows x k`.
    ///
    /// # Errors
    /// Returns [`StatsError::Shape`] when the feature count differs from the fit.
    pub fn transform(&self, data: &Matrix) -> Result<Matrix, StatsError> {
        if data.cols() != self.mean.len() {
            return Err(StatsError::Shape("feature count differs from the fitted model"));
        }
        let k = self.n_components();
        let mut scores = vec![0.0; data.rows() * k];
        for r in 0..data.rows() {
            let row = data.row(r);
            for c in 0..k {
                let axis = self.components.row(c);
                scores[r * k + c] =
                    row.iter().zip(&self.mean).zip(axis).map(|((x, m), a)| (x - m) * a).sum();
            }
        }
        Matrix::new(data.rows(), k, scores)
    }

    /// Number of retained components.
    #[must_use]
    pub fn n_components(&self) -> usize {
        self.components.rows()
    }

    /// The per-feature mean subtracted before projection.
    #[must_use]
    pub fn mean(&self) -> &[f64] {
        &self.mean
    }

    /// The retained principal axes as a `k x features` matrix (row per component).
    #[must_use]
    pub const fn components(&self) -> &Matrix {
        &self.components
    }

    /// Variance along each retained component, in descending order.
    #[must_use]
    pub fn explained_variance(&self) -> &[f64] {
        &self.explained_variance
    }

    /// Fraction of total variance along each retained component.
    #[must_use]
    pub fn explained_variance_ratio(&self) -> &[f64] {
        &self.explained_variance_ratio
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn recovers_known_principal_direction() {
        // Points spread along (1,1)/sqrt2 with a small orthogonal wobble along (1,-1)/sqrt2.
        let s = libm::sqrt(0.5);
        let mut rows = Vec::new();
        for t in -4..=4 {
            let t = f64::from(t);
            let e = if (t as i64) % 2 == 0 { 0.05 } else { -0.05 };
            rows.push(vec![(t + e) * s, (t - e) * s]);
        }
        let data = Matrix::from_rows(&rows).unwrap();
        let pca = Pca::fit(&data, 1).unwrap();
        let axis = pca.components().row(0);
        assert!((axis[0] - s).abs() < 1e-3, "axis {axis:?}");
        assert!((axis[1] - s).abs() < 1e-3, "axis {axis:?}");
        assert!(pca.explained_variance_ratio()[0] > 0.99);
    }

    #[test]
    fn scores_have_expected_variance_and_mean() {
        let data =
            Matrix::from_rows(&[vec![2.0, 0.0], vec![0.0, 1.0], vec![-2.0, 0.0], vec![0.0, -1.0]])
                .unwrap();
        let pca = Pca::fit(&data, 2).unwrap();
        // Covariance is diag(8/3, 2/3): trace 10/3, so the ratios are 0.8 and 0.2.
        assert!((pca.explained_variance()[0] - 8.0 / 3.0).abs() < 1e-9);
        assert!((pca.explained_variance_ratio()[0] - 0.8).abs() < 1e-9);
        assert!((pca.explained_variance_ratio()[1] - 0.2).abs() < 1e-9);
        let scores = pca.transform(&data).unwrap();
        assert_eq!(scores.rows(), 4);
        assert_eq!(scores.cols(), 2);
        // The first score of (2,0) is 2 along axis (1,0) (sign-normalized).
        assert!((scores.get(0, 0).abs() - 2.0).abs() < 1e-9);
    }

    #[test]
    fn rejects_invalid_arguments() {
        let data = Matrix::from_rows(&[vec![1.0, 2.0], vec![3.0, 5.0]]).unwrap();
        assert!(matches!(Pca::fit(&data, 0), Err(StatsError::InvalidArgument(_))));
        assert!(matches!(Pca::fit(&data, 3), Err(StatsError::InvalidArgument(_))));
        let one = Matrix::from_rows(&[vec![1.0, 2.0]]).unwrap();
        assert!(matches!(Pca::fit(&one, 1), Err(StatsError::Shape(_))));
    }
}
