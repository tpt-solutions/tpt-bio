//! UMAP: a fuzzy k-nearest-neighbour graph laid out by stochastic gradient descent.
//!
//! This is a simplified implementation, not a reference one. The neighbour search is
//! brute force, so it costs O(n²) distance evaluations. Each point's fuzzy membership
//! scale is found by bisection so that its memberships sum to `log2(k)`. Memberships
//! are symmetrised with the fuzzy union. The layout is optimised by stochastic
//! gradient descent with negative sampling, driven by a seeded generator.
//!
//! The curve parameters `a` and `b` in `1 / (1 + a x^(2b))` are fitted by a coarse grid
//! search to the target curve for the chosen `min_dist`, with spread fixed at 1. This
//! replaces the reference curve fit. Output has not been compared with umap-learn or
//! any other reference implementation.

use alloc::vec::Vec;

use super::rng::SplitMix64;
use super::{Matrix, StatsError};

/// Bisection steps for the per-point membership scale.
const SIGMA_STEPS: usize = 64;
/// Gradient components are clipped to this magnitude, as in the reference method.
const GRAD_CLIP: f64 = 4.0;
/// Learning rate at the start of the epoch schedule, decayed linearly to zero.
const INITIAL_ALPHA: f64 = 1.0;

/// Settings for [`umap`].
#[derive(Debug, Clone, PartialEq)]
pub struct UmapParams {
    /// Number of nearest neighbours. Must be at least 2 and below the row count.
    pub n_neighbors: usize,
    /// Minimum distance between embedded points, in `(0, 1]`.
    pub min_dist: f64,
    /// Number of passes over the edge list. Must be at least 1.
    pub n_epochs: usize,
    /// Repulsive samples drawn for each attractive update. Zero disables repulsion.
    pub negative_sample_rate: usize,
    /// Seed for the initial layout and for negative sampling.
    pub seed: u64,
}

impl Default for UmapParams {
    fn default() -> Self {
        Self { n_neighbors: 15, min_dist: 0.1, n_epochs: 200, negative_sample_rate: 5, seed: 0 }
    }
}

/// Embeds the rows of `data` in two dimensions with UMAP.
///
/// # Errors
/// Returns [`StatsError::InvalidArgument`] when there are fewer than three rows, when
/// `n_neighbors` is outside `[2, n)`, or when `min_dist` or `n_epochs` is not valid.
pub fn umap(data: &Matrix, params: &UmapParams) -> Result<Matrix, StatsError> {
    let n = data.rows();
    if n < 3 {
        return Err(StatsError::InvalidArgument("UMAP needs at least 3 rows"));
    }
    if !(params.n_neighbors >= 2 && params.n_neighbors < n) {
        return Err(StatsError::InvalidArgument("n_neighbors must be in [2, n)"));
    }
    if !(params.min_dist > 0.0 && params.min_dist <= 1.0) {
        return Err(StatsError::InvalidArgument("min_dist must be in (0, 1]"));
    }
    if params.n_epochs == 0 {
        return Err(StatsError::InvalidArgument("n_epochs must be at least 1"));
    }

    let knn = nearest_neighbours(data, params.n_neighbors);
    let edges = fuzzy_edges(&knn, params.n_neighbors);
    let (a, b) = fit_curve(params.min_dist);

    let mut rng = SplitMix64::new(params.seed);
    let mut y: Vec<f64> = (0..2 * n).map(|_| 2.0 * rng.next_unit() - 1.0).collect();

    // Each edge is sampled once every 1/weight epochs, so strong edges move more often.
    let period: Vec<f64> = edges.iter().map(|e| 1.0 / e.2).collect();
    let mut next_sample = period.clone();
    let epochs = params.n_epochs as f64;
    for epoch in 0..params.n_epochs {
        let alpha = INITIAL_ALPHA * (1.0 - epoch as f64 / epochs);
        let now = (epoch + 1) as f64;
        for (e, &(src, dst, _)) in edges.iter().enumerate() {
            if now < next_sample[e] {
                continue;
            }
            next_sample[e] += period[e];
            attract(&mut y, src, dst, a, b, alpha);
            for _ in 0..params.negative_sample_rate {
                let other = rng.below(n);
                if other != src {
                    repel(&mut y, src, other, a, b, alpha);
                }
            }
        }
    }

    Matrix::new(n, 2, y).map_err(|_| StatsError::Shape("UMAP layout"))
}

/// For each row, the indices and distances of its `k` nearest other rows, nearest first.
fn nearest_neighbours(data: &Matrix, k: usize) -> Vec<Vec<(usize, f64)>> {
    let n = data.rows();
    let mut out = Vec::with_capacity(n);
    for i in 0..n {
        let mut cand: Vec<(usize, f64)> =
            (0..n).filter(|&j| j != i).map(|j| (j, sq_dist(data.row(i), data.row(j)))).collect();
        cand.sort_by(|x, y| x.1.total_cmp(&y.1).then(x.0.cmp(&y.0)));
        cand.truncate(k);
        for entry in &mut cand {
            entry.1 = libm::sqrt(entry.1);
        }
        out.push(cand);
    }
    out
}

/// Symmetrised fuzzy edges `(low, high, weight)` with one entry per unordered pair.
fn fuzzy_edges(knn: &[Vec<(usize, f64)>], k: usize) -> Vec<(usize, usize, f64)> {
    let target = libm::log2(k as f64);
    // Directed memberships mu_ij for each kNN entry.
    let directed: Vec<Vec<(usize, f64)>> = knn
        .iter()
        .map(|row| {
            let rho = row.first().map_or(0.0, |&(_, d)| d);
            let sigma = membership_scale(row, rho, target);
            row.iter().map(|&(j, d)| (j, libm::exp(-(d - rho).max(0.0) / sigma))).collect()
        })
        .collect();

    let lookup = |from: usize, to: usize| -> f64 {
        directed[from].iter().find(|&&(j, _)| j == to).map_or(0.0, |&(_, w)| w)
    };
    let mut edges = Vec::new();
    for (i, row) in directed.iter().enumerate() {
        for &(j, wij) in row {
            let wji = lookup(j, i);
            // Emit each unordered pair once: from the lower index, or from the higher
            // index only when the lower index does not list it as a neighbour.
            if i < j || wji == 0.0 {
                let w = wij + wji - wij * wji;
                if w > 0.0 {
                    let (lo, hi) = if i < j { (i, j) } else { (j, i) };
                    edges.push((lo, hi, w));
                }
            }
        }
    }
    edges
}

/// Finds `sigma` so that `sum_j exp(-max(0, d_j - rho) / sigma)` equals `target`.
fn membership_scale(row: &[(usize, f64)], rho: f64, target: f64) -> f64 {
    let (mut lo, mut hi, mut sigma) = (0.0_f64, f64::INFINITY, 1.0_f64);
    for _ in 0..SIGMA_STEPS {
        let total: f64 = row.iter().map(|&(_, d)| libm::exp(-(d - rho).max(0.0) / sigma)).sum();
        if libm::fabs(total - target) < 1e-5 {
            break;
        }
        if total > target {
            hi = sigma;
            sigma = f64::midpoint(lo, hi);
        } else {
            lo = sigma;
            sigma = if hi.is_infinite() { sigma * 2.0 } else { f64::midpoint(lo, hi) };
        }
    }
    sigma.max(1e-12)
}

/// Fits `a` and `b` in `1 / (1 + a x^(2b))` to the target curve for `min_dist`.
fn fit_curve(min_dist: f64) -> (f64, f64) {
    let samples: Vec<(f64, f64)> = (1..=300)
        .map(|s| {
            let x = f64::from(s) * 0.01;
            let target = if x <= min_dist { 1.0 } else { libm::exp(-(x - min_dist)) };
            (x, target)
        })
        .collect();
    let loss = |a: f64, b: f64| -> f64 {
        samples
            .iter()
            .map(|&(x, t)| {
                let fit = 1.0 / (1.0 + a * libm::pow(x, 2.0 * b));
                (fit - t) * (fit - t)
            })
            .sum()
    };
    let mut best = (1.0, 1.0, f64::INFINITY);
    let (mut a_lo, mut a_hi, mut b_lo, mut b_hi) = (0.1, 4.0, 0.3, 2.0);
    // Two passes: a coarse grid over the full ranges, then a finer grid around the best.
    for _ in 0..2 {
        let steps = 80;
        let (a_step, b_step) = ((a_hi - a_lo) / steps as f64, (b_hi - b_lo) / steps as f64);
        for i in 0..=steps {
            for j in 0..=steps {
                let a = a_lo + a_step * i as f64;
                let b = b_lo + b_step * j as f64;
                let l = loss(a, b);
                if l < best.2 {
                    best = (a, b, l);
                }
            }
        }
        a_lo = (best.0 - a_step).max(1e-3);
        a_hi = best.0 + a_step;
        b_lo = (best.1 - b_step).max(1e-3);
        b_hi = best.1 + b_step;
    }
    (best.0, best.1)
}

/// Attractive update along an edge.
fn attract(y: &mut [f64], src: usize, dst: usize, a: f64, b: f64, alpha: f64) {
    let (dx, dy) = (y[2 * src] - y[2 * dst], y[2 * src + 1] - y[2 * dst + 1]);
    let d2 = dx * dx + dy * dy;
    if d2 <= 0.0 {
        return;
    }
    let coeff = -2.0 * a * b * libm::pow(d2, b - 1.0) / (1.0 + a * libm::pow(d2, b));
    let gx = clip(coeff * dx) * alpha;
    let gy = clip(coeff * dy) * alpha;
    y[2 * src] += gx;
    y[2 * src + 1] += gy;
    y[2 * dst] -= gx;
    y[2 * dst + 1] -= gy;
}

/// Repulsive update of `src` away from a sampled non-neighbour `other`.
fn repel(y: &mut [f64], src: usize, other: usize, a: f64, b: f64, alpha: f64) {
    let (dx, dy) = (y[2 * src] - y[2 * other], y[2 * src + 1] - y[2 * other + 1]);
    let d2 = dx * dx + dy * dy;
    let coeff = 2.0 * b / ((0.001 + d2) * (1.0 + a * libm::pow(d2, b)));
    y[2 * src] += clip(coeff * dx) * alpha;
    y[2 * src + 1] += clip(coeff * dy) * alpha;
}

fn clip(v: f64) -> f64 {
    v.clamp(-GRAD_CLIP, GRAD_CLIP)
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
        let data = two_blobs(n_per, 10, 10.0, 11);
        let params =
            UmapParams { n_neighbors: 10, n_epochs: 300, seed: 5, ..UmapParams::default() };
        let emb = umap(&data, &params).unwrap();
        assert_eq!((emb.rows(), emb.cols()), (2 * n_per, 2));
        let (between, within) = blob_separation(&emb, n_per);
        assert!(between > within, "between={between} within={within}");
    }

    #[test]
    fn same_seed_gives_same_layout() {
        let data = two_blobs(12, 5, 8.0, 4);
        let params = UmapParams { n_neighbors: 5, n_epochs: 50, ..UmapParams::default() };
        assert_eq!(umap(&data, &params).unwrap(), umap(&data, &params).unwrap());
    }

    #[test]
    fn curve_fit_is_close_to_reference_values_for_default_min_dist() {
        // The reference curve for min_dist = 0.1 is about a = 1.58, b = 0.90.
        let (a, b) = fit_curve(0.1);
        assert!((a - 1.58).abs() < 0.2, "a={a}");
        assert!((b - 0.90).abs() < 0.1, "b={b}");
    }

    #[test]
    fn rejects_bad_parameters() {
        let data = two_blobs(5, 3, 4.0, 2);
        let bad_k = UmapParams { n_neighbors: 10, ..UmapParams::default() };
        assert!(matches!(umap(&data, &bad_k), Err(StatsError::InvalidArgument(_))));
        let bad_dist = UmapParams { n_neighbors: 3, min_dist: 0.0, ..UmapParams::default() };
        assert!(matches!(umap(&data, &bad_dist), Err(StatsError::InvalidArgument(_))));
    }
}
