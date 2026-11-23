//! Clustering: k-means with k-means++ seeding, and agglomerative hierarchical clustering.
//!
//! Both are deterministic. k-means uses a `SplitMix64` generator seeded by the caller,
//! so a given seed always gives the same labels. Hierarchical clustering has no
//! randomness.

use alloc::vec;
use alloc::vec::Vec;

use super::rng::SplitMix64;
use super::{Matrix, StatsError};

/// Result of [`kmeans`].
#[derive(Debug, Clone, PartialEq)]
pub struct KMeans {
    /// Cluster index (`0..k`) for each row.
    pub labels: Vec<usize>,
    /// `k x features` centroid matrix.
    pub centroids: Matrix,
    /// Sum of squared distances from each row to its centroid.
    pub inertia: f64,
    /// Number of assignment passes performed.
    pub iterations: usize,
}

/// Linkage rule for [`agglomerative`]. Distances between clusters are Euclidean.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Linkage {
    /// Distance between the closest pair of members.
    Single,
    /// Distance between the farthest pair of members.
    Complete,
    /// Mean distance over all member pairs (UPGMA).
    Average,
}

fn sq_dist(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b).map(|(x, y)| (x - y) * (x - y)).sum()
}

/// Index of the nearest centroid. Ties go to the lower index.
fn nearest(row: &[f64], centroids: &Matrix) -> (usize, f64) {
    let mut best = (0, f64::INFINITY);
    for c in 0..centroids.rows() {
        let d = sq_dist(row, centroids.row(c));
        if d < best.1 {
            best = (c, d);
        }
    }
    best
}

/// Lloyd's algorithm with k-means++ seeding.
///
/// Stops when no label changes or after `max_iter` assignment passes. A cluster that
/// loses all its members keeps its previous centroid.
///
/// # Errors
/// Returns [`StatsError::InvalidArgument`] when `k` is zero, larger than the row
/// count, or `max_iter` is zero.
pub fn kmeans(data: &Matrix, k: usize, seed: u64, max_iter: usize) -> Result<KMeans, StatsError> {
    let (n, d) = (data.rows(), data.cols());
    if k == 0 || k > n {
        return Err(StatsError::InvalidArgument("k must be between 1 and the row count"));
    }
    if max_iter == 0 {
        return Err(StatsError::InvalidArgument("max_iter must be at least 1"));
    }
    let mut rng = SplitMix64::new(seed);

    // k-means++: first centre uniform, then each next centre sampled with
    // probability proportional to its squared distance from the nearest chosen centre.
    let mut centroids = Matrix::zeros(k, d);
    let first = rng.below(n);
    copy_row(&mut centroids, 0, data.row(first));
    let mut min_d2: Vec<f64> = (0..n).map(|i| sq_dist(data.row(i), data.row(first))).collect();
    for c in 1..k {
        let total: f64 = min_d2.iter().sum();
        let pick = if total > 0.0 {
            let target = rng.next_unit() * total;
            let mut acc = 0.0;
            let mut chosen = n - 1;
            for (i, &w) in min_d2.iter().enumerate() {
                acc += w;
                if acc > target {
                    chosen = i;
                    break;
                }
            }
            chosen
        } else {
            // Every row coincides with a centre already; any row will do.
            rng.below(n)
        };
        copy_row(&mut centroids, c, data.row(pick));
        for (i, w) in min_d2.iter_mut().enumerate() {
            let d2 = sq_dist(data.row(i), data.row(pick));
            if d2 < *w {
                *w = d2;
            }
        }
    }

    let mut labels = vec![usize::MAX; n];
    let mut iterations = 0;
    for _ in 0..max_iter {
        iterations += 1;
        let mut changed = false;
        for (i, label) in labels.iter_mut().enumerate() {
            let (c, _) = nearest(data.row(i), &centroids);
            if *label != c {
                *label = c;
                changed = true;
            }
        }
        if !changed {
            break;
        }
        let mut sums = Matrix::zeros(k, d);
        let mut counts = vec![0usize; k];
        for (i, &c) in labels.iter().enumerate() {
            counts[c] += 1;
            for (j, &x) in data.row(i).iter().enumerate() {
                let v = sums.get(c, j) + x;
                sums.set(c, j, v);
            }
        }
        for (c, &count) in counts.iter().enumerate() {
            if count > 0 {
                for j in 0..d {
                    centroids.set(c, j, sums.get(c, j) / count as f64);
                }
            }
        }
    }

    let inertia =
        labels.iter().enumerate().map(|(i, &c)| sq_dist(data.row(i), centroids.row(c))).sum();
    Ok(KMeans { labels, centroids, inertia, iterations })
}

fn copy_row(m: &mut Matrix, r: usize, src: &[f64]) {
    for (j, &x) in src.iter().enumerate() {
        m.set(r, j, x);
    }
}

/// Agglomerative hierarchical clustering, merging until `n_clusters` remain.
///
/// Uses Lance-Williams updates on the Euclidean distance matrix, so each merge is
/// exact for the chosen linkage. Runs in O(n^3) time, which suits the modest sample
/// counts of plate-level analysis. Ties are broken by the lowest index pair.
///
/// Returned labels are numbered `0..n_clusters` in order of each cluster's first row.
///
/// # Errors
/// Returns [`StatsError::InvalidArgument`] when `n_clusters` is zero or exceeds the row count.
pub fn agglomerative(
    data: &Matrix,
    linkage: Linkage,
    n_clusters: usize,
) -> Result<Vec<usize>, StatsError> {
    let n = data.rows();
    if n_clusters == 0 || n_clusters > n {
        return Err(StatsError::InvalidArgument("n_clusters must be between 1 and the row count"));
    }
    let mut dist = Matrix::zeros(n, n);
    for i in 0..n {
        for j in 0..n {
            let d2 = sq_dist(data.row(i), data.row(j));
            dist.set(i, j, libm::sqrt(d2));
        }
    }
    let mut size = vec![1usize; n];
    let mut active = vec![true; n];
    let mut owner: Vec<usize> = (0..n).collect();
    let mut remaining = n;

    while remaining > n_clusters {
        let mut best: Option<(usize, usize, f64)> = None;
        for (i, &alive_i) in active.iter().enumerate() {
            if !alive_i {
                continue;
            }
            for (j, &alive_j) in active.iter().enumerate().skip(i + 1) {
                if alive_j && best.is_none_or(|(_, _, d)| dist.get(i, j) < d) {
                    best = Some((i, j, dist.get(i, j)));
                }
            }
        }
        let Some((a, b, _)) = best else { break };
        for (k, &alive_k) in active.iter().enumerate() {
            if !alive_k || k == a || k == b {
                continue;
            }
            let (da, db) = (dist.get(a, k), dist.get(b, k));
            let merged = match linkage {
                Linkage::Single => da.min(db),
                Linkage::Complete => da.max(db),
                Linkage::Average => {
                    let (na, nb) = (size[a] as f64, size[b] as f64);
                    (na * da + nb * db) / (na + nb)
                }
            };
            dist.set(a, k, merged);
            dist.set(k, a, merged);
        }
        size[a] += size[b];
        active[b] = false;
        for o in &mut owner {
            if *o == b {
                *o = a;
            }
        }
        remaining -= 1;
    }

    let mut ids = vec![usize::MAX; n];
    let mut next = 0;
    let mut labels = vec![0usize; n];
    for i in 0..n {
        let root = owner[i];
        if ids[root] == usize::MAX {
            ids[root] = next;
            next += 1;
        }
        labels[i] = ids[root];
    }
    Ok(labels)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Three tight blobs of four points each, far apart. Rows 0-3, 4-7 and 8-11 belong together.
    fn blobs() -> Matrix {
        let centres = [(0.0, 0.0), (10.0, 0.0), (0.0, 10.0)];
        let offsets = [(0.1, 0.0), (-0.1, 0.0), (0.0, 0.1), (0.0, -0.1)];
        let mut rows = Vec::new();
        for &(cx, cy) in &centres {
            for &(ox, oy) in &offsets {
                rows.push(vec![cx + ox, cy + oy]);
            }
        }
        Matrix::from_rows(&rows).unwrap()
    }

    /// Checks that labels form exactly the blob partition, regardless of label numbering.
    fn is_blob_partition(labels: &[usize]) -> bool {
        let group = |g: usize| &labels[g * 4..g * 4 + 4];
        let each_uniform = (0..3).all(|g| group(g).iter().all(|&l| l == group(g)[0]));
        let distinct =
            group(0)[0] != group(1)[0] && group(0)[0] != group(2)[0] && group(1)[0] != group(2)[0];
        each_uniform && distinct
    }

    #[test]
    fn kmeans_separates_blobs_for_several_seeds() {
        let data = blobs();
        for seed in 0..20 {
            let res = kmeans(&data, 3, seed, 100).unwrap();
            assert!(is_blob_partition(&res.labels), "seed {seed}: {:?}", res.labels);
            // Each centroid sits at its blob centre to within the offset.
            assert!(res.inertia < 0.2, "inertia {}", res.inertia);
        }
    }

    #[test]
    fn kmeans_is_deterministic_for_a_seed() {
        let data = blobs();
        let a = kmeans(&data, 3, 42, 100).unwrap();
        let b = kmeans(&data, 3, 42, 100).unwrap();
        assert_eq!(a, b);
    }

    #[test]
    fn kmeans_rejects_bad_k() {
        let data = blobs();
        assert!(matches!(kmeans(&data, 0, 1, 10), Err(StatsError::InvalidArgument(_))));
        assert!(matches!(kmeans(&data, 13, 1, 10), Err(StatsError::InvalidArgument(_))));
        assert!(matches!(kmeans(&data, 2, 1, 0), Err(StatsError::InvalidArgument(_))));
    }

    #[test]
    fn agglomerative_recovers_blobs_for_every_linkage() {
        let data = blobs();
        for linkage in [Linkage::Single, Linkage::Complete, Linkage::Average] {
            let labels = agglomerative(&data, linkage, 3).unwrap();
            assert!(is_blob_partition(&labels), "{linkage:?}: {labels:?}");
            // Numbering is by first appearance.
            assert_eq!(labels[0], 0);
        }
    }

    #[test]
    fn agglomerative_extremes() {
        let data = blobs();
        assert_eq!(
            agglomerative(&data, Linkage::Average, 12).unwrap(),
            (0..12).collect::<Vec<_>>()
        );
        assert!(agglomerative(&data, Linkage::Single, 1).unwrap().iter().all(|&l| l == 0));
    }

    #[test]
    fn average_linkage_matches_hand_computed_merge() {
        // Points 0, 1 at x=0 and x=1; point 2 at x=5. Average-linkage with k=2
        // merges {0,1} first, then keeps 2 alone.
        let data = Matrix::from_rows(&[vec![0.0], vec![1.0], vec![5.0]]).unwrap();
        assert_eq!(agglomerative(&data, Linkage::Average, 2).unwrap(), vec![0, 0, 1]);
        // Single linkage on the chain 0 - 1 - 2.5 with spacings 1 and 1.5 merges 0 and 1 first.
        let chain = Matrix::from_rows(&[vec![0.0], vec![1.0], vec![2.5]]).unwrap();
        assert_eq!(agglomerative(&chain, Linkage::Single, 2).unwrap(), vec![0, 0, 1]);
        // Complete linkage also merges the closest pair (distance 1) first.
        assert_eq!(agglomerative(&chain, Linkage::Complete, 2).unwrap(), vec![0, 0, 1]);
    }
}
