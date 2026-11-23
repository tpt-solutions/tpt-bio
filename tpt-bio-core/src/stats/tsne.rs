//! t-SNE: gradient descent on the Kullback-Leibler objective, computed exactly.
//!
//! This is a simplified implementation, not a reference one. Every pairwise distance
//! is computed, so each iteration costs O(n²) time and memory. Each point's Gaussian
//! precision is found by bisection so that its conditional distribution has the
//! requested perplexity. The optimiser is gradient descent with momentum and
//! per-coordinate gains, early exaggeration for the first 250 iterations, and a
//! seeded random start. Output has not been compared with scikit-learn or any other
//! reference implementation, so the coordinates will not match theirs.

use alloc::vec;
use alloc::vec::Vec;

use super::rng::SplitMix64;
use super::{Matrix, StatsError};

/// Number of initial iterations run with early exaggeration.
const EXAGGERATION_ITERS: usize = 250;
/// Early exaggeration factor applied to the input affinities.
const EXAGGERATION: f64 = 12.0;
/// Half-width of the uniform initial layout.
const INIT_SCALE: f64 = 1e-4;
/// Bisection steps for the perplexity calibration.
const CALIBRATION_STEPS: usize = 100;

/// Settings for [`tsne`].
#[derive(Debug, Clone, PartialEq)]
pub struct TsneParams {
    /// Effective number of neighbours. Must be at least 1 and below `n - 1`.
    pub perplexity: f64,
    /// Number of gradient steps. Must be at least 1.
    pub iterations: usize,
    /// Step size for gradient descent. Must be positive and finite.
    pub learning_rate: f64,
    /// Seed for the random initial layout.
    pub seed: u64,
}

impl Default for TsneParams {
    fn default() -> Self {
        Self { perplexity: 30.0, iterations: 1000, learning_rate: 200.0, seed: 0 }
    }
}

/// Output of [`tsne`].
#[derive(Debug, Clone, PartialEq)]
pub struct TsneResult {
    /// The `n x 2` layout, one row per input row.
    pub embedding: Matrix,
    /// The final KL divergence between the input and output affinities, measured
    /// without early exaggeration.
    pub kl_divergence: f64,
}

/// Embeds the rows of `data` in two dimensions with t-SNE.
///
/// # Errors
/// Returns [`StatsError::InvalidArgument`] when there are fewer than three rows, when
/// the perplexity is outside `[1, n - 1)`, or when the iteration count or learning
/// rate is not valid.
#[allow(clippy::many_single_char_names)]
pub fn tsne(data: &Matrix, params: &TsneParams) -> Result<TsneResult, StatsError> {
    let n = data.rows();
    if n < 3 {
        return Err(StatsError::InvalidArgument("t-SNE needs at least 3 rows"));
    }
    if !(params.perplexity >= 1.0 && params.perplexity < (n - 1) as f64) {
        return Err(StatsError::InvalidArgument("perplexity must be in [1, n - 1)"));
    }
    if params.iterations == 0 {
        return Err(StatsError::InvalidArgument("iterations must be at least 1"));
    }
    if !(params.learning_rate > 0.0 && params.learning_rate.is_finite()) {
        return Err(StatsError::InvalidArgument("learning rate must be positive and finite"));
    }

    let p = joint_affinities(data, params.perplexity);
    let mut rng = SplitMix64::new(params.seed);
    let mut y: Vec<f64> = (0..2 * n).map(|_| (2.0 * rng.next_unit() - 1.0) * INIT_SCALE).collect();
    let mut velocity = vec![0.0; 2 * n];
    let mut gains: Vec<f64> = vec![1.0; 2 * n];
    let mut grad = vec![0.0; 2 * n];
    let mut w = vec![0.0; n * n];

    for iter in 0..params.iterations {
        let (exaggeration, momentum) =
            if iter < EXAGGERATION_ITERS { (EXAGGERATION, 0.5) } else { (1.0, 0.8) };
        let z = student_weights(&y, n, &mut w);

        grad.fill(0.0);
        for i in 0..n {
            for j in 0..n {
                if i == j {
                    continue;
                }
                let wij = w[i * n + j];
                let mult = 4.0 * (exaggeration * p[i * n + j] - wij / z) * wij;
                grad[2 * i] += mult * (y[2 * i] - y[2 * j]);
                grad[2 * i + 1] += mult * (y[2 * i + 1] - y[2 * j + 1]);
            }
        }

        for k in 0..2 * n {
            // Gains grow while the gradient keeps its sign and shrink when it flips.
            gains[k] = if (grad[k] > 0.0) == (velocity[k] > 0.0) {
                (gains[k] * 0.8).max(0.01)
            } else {
                gains[k] + 0.2
            };
            velocity[k] = momentum * velocity[k] - params.learning_rate * gains[k] * grad[k];
            y[k] += velocity[k];
        }
        center(&mut y);
    }

    let z = student_weights(&y, n, &mut w);
    let mut kl = 0.0;
    for i in 0..n {
        for j in 0..n {
            let pij = p[i * n + j];
            if i != j && pij > 0.0 {
                kl += pij * libm::log(pij / (w[i * n + j] / z));
            }
        }
    }

    let embedding = Matrix::new(n, 2, y).map_err(|_| StatsError::Shape("t-SNE layout"))?;
    Ok(TsneResult { embedding, kl_divergence: kl })
}

/// Symmetrised joint affinities `P_ij = (p_j|i + p_i|j) / (2n)`, as a flat `n x n` array.
fn joint_affinities(data: &Matrix, perplexity: f64) -> Vec<f64> {
    let n = data.rows();
    let mut cond = vec![0.0; n * n];
    let mut others: Vec<usize> = Vec::with_capacity(n - 1);
    let mut dists: Vec<f64> = Vec::with_capacity(n - 1);
    let mut row = vec![0.0; n - 1];
    for i in 0..n {
        others.clear();
        dists.clear();
        for j in (0..n).filter(|&j| j != i) {
            others.push(j);
            dists.push(sq_dist(data.row(i), data.row(j)));
        }
        conditional_row(&dists, perplexity, &mut row);
        for (&j, &pj) in others.iter().zip(&row) {
            cond[i * n + j] = pj;
        }
    }
    let scale = 1.0 / (2.0 * n as f64);
    let mut p = vec![0.0; n * n];
    for i in 0..n {
        for j in 0..n {
            p[i * n + j] = (cond[i * n + j] + cond[j * n + i]) * scale;
        }
    }
    p
}

/// Gaussian conditional probabilities over the squared distances `sq_dists`, with the
/// given perplexity, written to `out`. The precision is found by bisection on entropy.
fn conditional_row(sq_dists: &[f64], perplexity: f64, out: &mut [f64]) {
    let target = libm::log(perplexity);
    let dmin = sq_dists.iter().copied().fold(f64::INFINITY, f64::min);
    let mut beta: f64 = 1.0;
    let (mut lo, mut hi) = (0.0_f64, f64::INFINITY);
    for _ in 0..CALIBRATION_STEPS {
        let mut sum = 0.0;
        let mut weighted = 0.0;
        for (o, &d) in out.iter_mut().zip(sq_dists) {
            // Shifting by the minimum keeps the largest term at 1 and avoids underflow.
            let shifted = d - dmin;
            *o = libm::exp(-beta * shifted);
            sum += *o;
            weighted += *o * shifted;
        }
        // Entropy in nats of the normalised distribution.
        let entropy = libm::log(sum) + beta * weighted / sum;
        let diff = entropy - target;
        if libm::fabs(diff) < 1e-5 {
            break;
        }
        if diff > 0.0 {
            lo = beta;
            beta = if hi.is_infinite() { beta * 2.0 } else { f64::midpoint(beta, hi) };
        } else {
            hi = beta;
            beta = f64::midpoint(beta, lo);
        }
    }
    let sum: f64 = out.iter().sum();
    for o in out.iter_mut() {
        *o /= sum;
    }
}

/// Fills `w` with the Student-t kernel `1 / (1 + d²)` (zero on the diagonal) and
/// returns its total over all `i != j`.
#[allow(clippy::many_single_char_names)]
fn student_weights(y: &[f64], n: usize, w: &mut [f64]) -> f64 {
    let mut z = 0.0;
    for i in 0..n {
        w[i * n + i] = 0.0;
        for j in (i + 1)..n {
            let dx = y[2 * i] - y[2 * j];
            let dy = y[2 * i + 1] - y[2 * j + 1];
            let v = 1.0 / (1.0 + dx * dx + dy * dy);
            w[i * n + j] = v;
            w[j * n + i] = v;
            z += 2.0 * v;
        }
    }
    z
}

/// Subtracts the mean so the layout stays centred at the origin.
fn center(y: &mut [f64]) {
    let n = (y.len() / 2) as f64;
    let (sx, sy) = y.chunks_exact(2).fold((0.0, 0.0), |(x, v), c| (x + c[0], v + c[1]));
    let (mx, my) = (sx / n, sy / n);
    for c in y.chunks_exact_mut(2) {
        c[0] -= mx;
        c[1] -= my;
    }
}

fn sq_dist(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b).map(|(x, y)| (x - y) * (x - y)).sum()
}

#[cfg(test)]
mod tests {
    use super::super::rng::test_support::{blob_separation, two_blobs};
    use super::*;

    #[test]
    fn separates_two_gaussian_blobs_in_10d() {
        let n_per = 30;
        let data = two_blobs(n_per, 10, 10.0, 7);
        let params =
            TsneParams { perplexity: 10.0, iterations: 500, learning_rate: 100.0, seed: 3 };
        let result = tsne(&data, &params).unwrap();
        assert_eq!((result.embedding.rows(), result.embedding.cols()), (2 * n_per, 2));
        let (between, within) = blob_separation(&result.embedding, n_per);
        assert!(between > within, "between={between} within={within}");
        assert!(result.kl_divergence.is_finite() && result.kl_divergence >= 0.0);
    }

    #[test]
    fn same_seed_gives_same_layout() {
        let data = two_blobs(10, 4, 6.0, 1);
        let params = TsneParams { perplexity: 5.0, iterations: 50, ..TsneParams::default() };
        assert_eq!(tsne(&data, &params).unwrap(), tsne(&data, &params).unwrap());
    }

    #[test]
    fn rejects_bad_perplexity() {
        let data = two_blobs(5, 3, 4.0, 2);
        let too_big = TsneParams { perplexity: 9.0, ..TsneParams::default() };
        assert!(matches!(tsne(&data, &too_big), Err(StatsError::InvalidArgument(_))));
        let too_small = TsneParams { perplexity: 0.5, ..TsneParams::default() };
        assert!(matches!(tsne(&data, &too_small), Err(StatsError::InvalidArgument(_))));
    }

    #[test]
    fn calibrated_rows_sum_to_one_and_favour_near_points() {
        let d = [1.0, 4.0, 9.0, 16.0, 25.0];
        let mut out = [0.0; 5];
        conditional_row(&d, 2.0, &mut out);
        assert!((out.iter().sum::<f64>() - 1.0).abs() < 1e-12);
        assert!(out[0] > out[4]);
    }
}
