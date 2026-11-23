//! Dense row-major matrix and a Jacobi eigensolver for symmetric matrices.

use alloc::vec;
use alloc::vec::Vec;

use super::StatsError;

/// Maximum number of Jacobi sweeps before giving up on convergence.
const MAX_SWEEPS: usize = 100;

/// A dense row-major matrix of finite `f64` values.
#[derive(Debug, Clone, PartialEq)]
///
/// # Examples
///
/// ```
/// use tpt_bio_core::stats::Matrix;
///
/// let m = Matrix::from_rows(&[vec![1.0, 2.0], vec![3.0, 4.0]]).unwrap();
/// assert_eq!((m.rows(), m.cols()), (2, 2));
/// assert_eq!(m.get(1, 0), 3.0);
/// ```
pub struct Matrix {
    rows: usize,
    cols: usize,
    data: Vec<f64>,
}

impl Matrix {
    /// Builds a matrix from row-major data.
    ///
    /// # Errors
    /// Returns [`StatsError::Empty`] for a zero dimension, [`StatsError::Shape`] when
    /// `data.len() != rows * cols`, and [`StatsError::NonFinite`] for NaN or infinity.
    pub fn new(rows: usize, cols: usize, data: Vec<f64>) -> Result<Self, StatsError> {
        if rows == 0 || cols == 0 {
            return Err(StatsError::Empty);
        }
        if rows.checked_mul(cols) != Some(data.len()) {
            return Err(StatsError::Shape("data length does not equal rows * cols"));
        }
        if data.iter().any(|v| !v.is_finite()) {
            return Err(StatsError::NonFinite);
        }
        Ok(Self { rows, cols, data })
    }

    /// Builds a matrix from a slice of rows. All rows must have the same length.
    ///
    /// # Errors
    /// Same as [`Matrix::new`], plus [`StatsError::Shape`] for ragged rows.
    pub fn from_rows(rows: &[Vec<f64>]) -> Result<Self, StatsError> {
        let cols = rows.first().map_or(0, Vec::len);
        if rows.iter().any(|r| r.len() != cols) {
            return Err(StatsError::Shape("rows have different lengths"));
        }
        Self::new(rows.len(), cols, rows.concat())
    }

    /// A matrix of zeros. Used internally where the shape is already known to be valid.
    pub(crate) fn zeros(rows: usize, cols: usize) -> Self {
        Self { rows, cols, data: vec![0.0; rows * cols] }
    }

    /// Number of rows.
    #[must_use]
    pub const fn rows(&self) -> usize {
        self.rows
    }

    /// Number of columns.
    #[must_use]
    pub const fn cols(&self) -> usize {
        self.cols
    }

    /// The element at `(row, col)`.
    ///
    /// # Panics
    /// Panics if either index is out of range.
    #[must_use]
    pub fn get(&self, row: usize, col: usize) -> f64 {
        assert!(row < self.rows && col < self.cols, "matrix index out of range");
        self.data[row * self.cols + col]
    }

    /// Row `r` as a slice.
    ///
    /// # Panics
    /// Panics if `r` is out of range.
    #[must_use]
    pub fn row(&self, r: usize) -> &[f64] {
        &self.data[r * self.cols..(r + 1) * self.cols]
    }

    /// Column `c` copied into a new vector.
    ///
    /// # Panics
    /// Panics if `c` is out of range.
    #[must_use]
    pub fn column(&self, c: usize) -> Vec<f64> {
        (0..self.rows).map(|r| self.data[r * self.cols + c]).collect()
    }

    /// The row-major backing data.
    #[must_use]
    pub fn as_slice(&self) -> &[f64] {
        &self.data
    }

    pub(crate) fn set(&mut self, row: usize, col: usize, value: f64) {
        self.data[row * self.cols + col] = value;
    }

    pub(crate) fn as_mut_slice(&mut self) -> &mut [f64] {
        &mut self.data
    }

    /// Consumes the matrix and returns the row-major data.
    #[must_use]
    pub fn into_data(self) -> Vec<f64> {
        self.data
    }
}

/// Eigendecomposition of a symmetric matrix, sorted by descending eigenvalue.
#[derive(Debug, Clone, PartialEq)]
pub struct Eigen {
    /// Eigenvalues in descending order.
    pub values: Vec<f64>,
    /// Eigenvectors stored as columns: column `j` belongs to `values[j]`.
    pub vectors: Matrix,
}

/// Eigendecomposition of a real symmetric matrix by cyclic Jacobi rotations.
///
/// Only the values of the input are used; symmetry is assumed, not checked.
///
/// # Errors
/// Returns [`StatsError::Shape`] if the matrix is not square.
#[allow(clippy::many_single_char_names)] // Jacobi indices p, q, k follow the textbook notation.
pub fn symmetric_eigen(a: &Matrix) -> Result<Eigen, StatsError> {
    let n = a.rows;
    if n != a.cols {
        return Err(StatsError::Shape("matrix must be square"));
    }
    let mut m = a.data.clone();
    let mut v = Matrix::zeros(n, n);
    for i in 0..n {
        v.set(i, i, 1.0);
    }

    for _ in 0..MAX_SWEEPS {
        let mut off = 0.0;
        let mut total = 0.0;
        for p in 0..n {
            for q in 0..n {
                let x = m[p * n + q];
                total += x * x;
                if p != q {
                    off += x * x;
                }
            }
        }
        if off <= 1e-30 * total || off == 0.0 {
            break;
        }
        for p in 0..n {
            for q in (p + 1)..n {
                let apq = m[p * n + q];
                if apq == 0.0 {
                    continue;
                }
                let app = m[p * n + p];
                let aqq = m[q * n + q];
                let theta = (aqq - app) / (2.0 * apq);
                let sign = if theta >= 0.0 { 1.0 } else { -1.0 };
                let t = sign / (theta.abs() + libm::sqrt(theta * theta + 1.0));
                let c = 1.0 / libm::sqrt(t * t + 1.0);
                let s = t * c;
                for k in 0..n {
                    if k != p && k != q {
                        let akp = m[k * n + p];
                        let akq = m[k * n + q];
                        m[k * n + p] = c * akp - s * akq;
                        m[p * n + k] = m[k * n + p];
                        m[k * n + q] = c * akq + s * akp;
                        m[q * n + k] = m[k * n + q];
                    }
                    let vkp = v.data[k * n + p];
                    let vkq = v.data[k * n + q];
                    v.data[k * n + p] = c * vkp - s * vkq;
                    v.data[k * n + q] = c * vkq + s * vkp;
                }
                m[p * n + p] = app - t * apq;
                m[q * n + q] = aqq + t * apq;
                m[p * n + q] = 0.0;
                m[q * n + p] = 0.0;
            }
        }
    }

    let mut order: Vec<usize> = (0..n).collect();
    order.sort_by(|&i, &j| m[j * n + j].total_cmp(&m[i * n + i]));
    let values = order.iter().map(|&i| m[i * n + i]).collect();
    let mut vectors = Matrix::zeros(n, n);
    for (new_col, &old_col) in order.iter().enumerate() {
        for r in 0..n {
            vectors.set(r, new_col, v.data[r * n + old_col]);
        }
    }
    Ok(Eigen { values, vectors })
}

#[cfg(test)]
mod tests {
    use super::*;

    fn close(a: f64, b: f64) -> bool {
        (a - b).abs() < 1e-9
    }

    #[test]
    fn rejects_bad_shapes() {
        assert_eq!(Matrix::new(0, 3, vec![]), Err(StatsError::Empty));
        assert!(matches!(Matrix::new(2, 2, vec![1.0; 3]), Err(StatsError::Shape(_))));
        assert_eq!(Matrix::new(1, 1, vec![f64::NAN]), Err(StatsError::NonFinite));
        assert!(matches!(
            Matrix::from_rows(&[vec![1.0, 2.0], vec![3.0]]),
            Err(StatsError::Shape(_))
        ));
        let rect = Matrix::new(1, 2, vec![1.0, 2.0]).unwrap();
        assert!(matches!(symmetric_eigen(&rect), Err(StatsError::Shape(_))));
    }

    #[test]
    fn two_by_two_known_eigenvalues() {
        // [[2, 1], [1, 2]] has eigenvalues 3 and 1 with eigenvectors (1,1)/sqrt2 and (1,-1)/sqrt2.
        let a = Matrix::new(2, 2, vec![2.0, 1.0, 1.0, 2.0]).unwrap();
        let e = symmetric_eigen(&a).unwrap();
        assert!(close(e.values[0], 3.0));
        assert!(close(e.values[1], 1.0));
        let r = libm::sqrt(0.5);
        assert!(close(e.vectors.get(0, 0).abs(), r));
        assert!(close(e.vectors.get(1, 0).abs(), r));
        assert!(close(e.vectors.get(0, 0), e.vectors.get(1, 0)));
    }

    #[test]
    fn three_by_three_reconstructs_matrix() {
        let data = vec![4.0, 1.0, -2.0, 1.0, 2.0, 0.0, -2.0, 0.0, 3.0];
        let a = Matrix::new(3, 3, data.clone()).unwrap();
        let e = symmetric_eigen(&a).unwrap();
        // A = V diag(values) V^T
        for i in 0..3 {
            for j in 0..3 {
                let s: f64 =
                    (0..3).map(|k| e.vectors.get(i, k) * e.values[k] * e.vectors.get(j, k)).sum();
                assert!(close(s, data[i * 3 + j]), "entry ({i},{j}) = {s}");
            }
        }
        assert!(e.values.windows(2).all(|w| w[0] >= w[1]));
    }
}
