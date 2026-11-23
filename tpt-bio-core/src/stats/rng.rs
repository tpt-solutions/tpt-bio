//! `SplitMix64`, the seeded generator shared by the clustering and embedding code.
//!
//! It is used only for seeding and for sampling in the stochastic routines. It is not
//! a cryptographic generator.

/// `SplitMix64`: a small, well-known generator with a 64-bit state.
#[derive(Debug, Clone)]
pub(crate) struct SplitMix64(u64);

impl SplitMix64 {
    pub(crate) const fn new(seed: u64) -> Self {
        Self(seed)
    }

    pub(crate) fn next_u64(&mut self) -> u64 {
        self.0 = self.0.wrapping_add(0x9E37_79B9_7F4A_7C15);
        let mut z = self.0;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
        z ^ (z >> 31)
    }

    /// Uniform in `[0, 1)` with 53 bits of precision.
    pub(crate) fn next_unit(&mut self) -> f64 {
        (self.next_u64() >> 11) as f64 * (1.0 / (1u64 << 53) as f64)
    }

    /// Uniform index in `0..n`. The modulo bias is negligible for the sizes used here.
    pub(crate) fn below(&mut self, n: usize) -> usize {
        (self.next_u64() % n as u64) as usize
    }

    /// Standard normal deviate by the Box-Muller transform.
    #[cfg(test)]
    pub(crate) fn next_gaussian(&mut self) -> f64 {
        // Map to (0, 1] so the logarithm is finite.
        let u1 = 1.0 - self.next_unit();
        let u2 = self.next_unit();
        libm::sqrt(-2.0 * libm::log(u1)) * libm::cos(core::f64::consts::TAU * u2)
    }
}

/// Test fixtures shared by the embedding tests.
#[cfg(test)]
pub(crate) mod test_support {
    use super::SplitMix64;
    use crate::stats::Matrix;
    use alloc::vec::Vec;

    /// Two standard-normal blobs in `dim` dimensions with `n_per` points each.
    /// The second blob is shifted by `shift` along the first axis.
    pub(crate) fn two_blobs(n_per: usize, dim: usize, shift: f64, seed: u64) -> Matrix {
        let mut rng = SplitMix64::new(seed);
        let mut rows: Vec<Vec<f64>> = Vec::new();
        for blob in 0..2 {
            for _ in 0..n_per {
                let mut row: Vec<f64> = (0..dim).map(|_| rng.next_gaussian()).collect();
                if blob == 1 {
                    row[0] += shift;
                }
                rows.push(row);
            }
        }
        Matrix::from_rows(&rows).unwrap()
    }

    /// Returns `(distance between blob centroids, mean distance to own centroid)` for a
    /// 2D embedding whose first `n_per` rows are one blob and the rest the other.
    pub(crate) fn blob_separation(emb: &Matrix, n_per: usize) -> (f64, f64) {
        let centroid = |range: core::ops::Range<usize>| -> [f64; 2] {
            let m = range.len() as f64;
            let (mut x, mut y) = (0.0, 0.0);
            for r in range {
                x += emb.get(r, 0);
                y += emb.get(r, 1);
            }
            [x / m, y / m]
        };
        let dist = |a: [f64; 2], b: [f64; 2]| {
            let (dx, dy) = (a[0] - b[0], a[1] - b[1]);
            libm::sqrt(dx * dx + dy * dy)
        };
        let c0 = centroid(0..n_per);
        let c1 = centroid(n_per..2 * n_per);
        let mut within = 0.0;
        for r in 0..2 * n_per {
            let c = if r < n_per { c0 } else { c1 };
            within += dist([emb.get(r, 0), emb.get(r, 1)], c);
        }
        (dist(c0, c1), within / (2 * n_per) as f64)
    }
}
