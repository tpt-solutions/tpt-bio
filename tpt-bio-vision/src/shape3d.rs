//! Shape features of a 3D binary object: volume, surface area, bounding box and
//! principal axis lengths.
//!
//! Voxel `(x, y, z)` occupies the physical point `(x * sx, y * sy, z * sz)` for voxel
//! spacing `(sx, sy, sz)`.
//!
//! * `volume` is the foreground voxel count times the voxel volume.
//! * `surface_area` is a face-count estimate, not a marching-cubes surface: it sums the
//!   area of every voxel face that borders background or lies outside the image. A
//!   staircase of voxels has this area, which overestimates the area of smooth curved
//!   surfaces (a digital sphere comes out larger than `4 pi r^2`).
//! * `principal_axis_lengths` are `4 * sqrt(lambda)` for the eigenvalues `lambda` of the
//!   population covariance matrix of the voxel centres, sorted in descending order. This
//!   is the scikit-image convention, which gives the full axis length of a uniform
//!   ellipsoid in 2D. The 3D value for a uniform solid differs from the semi-axis by a
//!   constant factor.

use crate::image::{FACE_OFFSETS, Image};

/// Errors from shape feature extraction.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum ShapeError {
    /// The mask has no foreground voxels.
    EmptyObject,
    /// A voxel spacing is not finite and positive.
    InvalidSpacing,
}

impl std::fmt::Display for ShapeError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::EmptyObject => f.write_str("mask has no foreground voxels"),
            Self::InvalidSpacing => f.write_str("voxel spacing must be finite and positive"),
        }
    }
}

impl std::error::Error for ShapeError {}

/// Shape features of one 3D object.
#[derive(Debug, Clone, PartialEq)]
pub struct ShapeFeatures3d {
    /// Volume in physical units.
    pub volume: f64,
    /// Surface area estimate in physical units (exposed voxel faces).
    pub surface_area: f64,
    /// Inclusive bounding box `(min, max)` in voxel coordinates `(x, y, z)`.
    pub bounding_box: ([usize; 3], [usize; 3]),
    /// Principal axis lengths in physical units, in descending order.
    pub principal_axis_lengths: [f64; 3],
}

/// Shape features of the foreground of `mask` with voxel spacing `spacing = (sx, sy, sz)`.
///
/// # Errors
/// Returns [`ShapeError::EmptyObject`] if there is no foreground, and
/// [`ShapeError::InvalidSpacing`] if a spacing is not finite and positive.
pub fn shape_features_3d(
    mask: &Image<bool>,
    spacing: (f64, f64, f64),
) -> Result<ShapeFeatures3d, ShapeError> {
    let (sx, sy, sz) = spacing;
    if ![sx, sy, sz].iter().all(|s| s.is_finite() && *s > 0.0) {
        return Err(ShapeError::InvalidSpacing);
    }
    // Area of a face whose normal is along x, y or z, in the order of FACE_OFFSETS.
    let face_area = [sy * sz, sy * sz, sx * sz, sx * sz, sx * sy, sx * sy];
    let voxel_volume = sx * sy * sz;

    let mut count = 0usize;
    let mut surface = 0.0f64;
    let mut lo = [usize::MAX; 3];
    let mut hi = [0usize; 3];
    let mut sum = [0.0f64; 3];
    let mut sum_sq = [[0.0f64; 3]; 3];

    for idx in 0..mask.len() {
        if !mask.data()[idx] {
            continue;
        }
        count += 1;
        let (x, y, z) = mask.coords(idx);
        let c = [x, y, z];
        for axis in 0..3 {
            lo[axis] = lo[axis].min(c[axis]);
            hi[axis] = hi[axis].max(c[axis]);
        }
        for (k, &d) in FACE_OFFSETS.iter().enumerate() {
            let exposed = mask.offset(idx, d).is_none_or(|n| !mask.data()[n]);
            if exposed {
                surface += face_area[k];
            }
        }
        let p = [x as f64 * sx, y as f64 * sy, z as f64 * sz];
        for a in 0..3 {
            sum[a] += p[a];
            for b in 0..3 {
                sum_sq[a][b] += p[a] * p[b];
            }
        }
    }
    if count == 0 {
        return Err(ShapeError::EmptyObject);
    }

    let n = count as f64;
    let mean = [sum[0] / n, sum[1] / n, sum[2] / n];
    let mut cov = [[0.0f64; 3]; 3];
    for a in 0..3 {
        for b in 0..3 {
            cov[a][b] = sum_sq[a][b] / n - mean[a] * mean[b];
        }
    }
    let eig = symmetric_eigenvalues_3x3(cov);
    let axes = eig.map(|l| 4.0 * l.max(0.0).sqrt());

    Ok(ShapeFeatures3d {
        volume: n * voxel_volume,
        surface_area: surface,
        bounding_box: (lo, hi),
        principal_axis_lengths: axes,
    })
}

/// Eigenvalues of a symmetric 3x3 matrix, by cyclic Jacobi rotations, in descending order.
///
/// Only the upper triangle is read; the matrix is treated as symmetric.
#[must_use]
#[allow(clippy::needless_range_loop)]
pub fn symmetric_eigenvalues_3x3(mut a: [[f64; 3]; 3]) -> [f64; 3] {
    // Symmetrise from the upper triangle so rounding in the input cannot matter.
    for p in 0..3 {
        for q in 0..p {
            a[p][q] = a[q][p];
        }
    }
    for _sweep in 0..100 {
        let off: f64 = (0..3)
            .flat_map(|p| (p + 1..3).map(move |q| (p, q)))
            .map(|(p, q)| a[p][q] * a[p][q])
            .sum();
        let scale: f64 = a.iter().flatten().map(|v| v * v).sum::<f64>().max(f64::MIN_POSITIVE);
        if off <= 1e-30 * scale {
            break;
        }
        for p in 0..3 {
            for q in p + 1..3 {
                if a[p][q] == 0.0 {
                    continue;
                }
                let theta = (a[q][q] - a[p][p]) / (2.0 * a[p][q]);
                let t = if theta >= 0.0 {
                    1.0 / (theta + (theta * theta + 1.0).sqrt())
                } else {
                    -1.0 / (-theta + (theta * theta + 1.0).sqrt())
                };
                let c = 1.0 / (t * t + 1.0).sqrt();
                let s = t * c;
                // A <- J^T A J with the rotation in the (p, q) plane.
                for k in 0..3 {
                    let (akp, akq) = (a[k][p], a[k][q]);
                    a[k][p] = c * akp - s * akq;
                    a[k][q] = s * akp + c * akq;
                }
                for k in 0..3 {
                    let (apk, aqk) = (a[p][k], a[q][k]);
                    a[p][k] = c * apk - s * aqk;
                    a[q][k] = s * apk + c * aqk;
                }
            }
        }
    }
    let mut eig = [a[0][0], a[1][1], a[2][2]];
    eig.sort_by(|x, y| y.total_cmp(x));
    eig
}

#[allow(clippy::float_cmp)]
#[cfg(test)]
mod tests {
    use super::*;

    /// Solid box with `lens = (nx, ny, nz)` voxels, placed at `origin`, in a volume of `dims`.
    fn boxed(
        dims: (usize, usize, usize),
        origin: (usize, usize, usize),
        lens: (usize, usize, usize),
    ) -> Image<bool> {
        let (w, h, d) = dims;
        let mut data = vec![false; w * h * d];
        for z in origin.2..origin.2 + lens.2 {
            for y in origin.1..origin.1 + lens.1 {
                for x in origin.0..origin.0 + lens.0 {
                    data[(z * h + y) * w + x] = true;
                }
            }
        }
        Image::new_3d(w, h, d, data).unwrap()
    }

    #[test]
    fn jacobi_matches_known_eigenvalues() {
        // [[2,1,0],[1,2,0],[0,0,3]] has eigenvalues 3, 3 and 1.
        let eig = symmetric_eigenvalues_3x3([[2.0, 1.0, 0.0], [1.0, 2.0, 0.0], [0.0, 0.0, 3.0]]);
        for (got, want) in eig.iter().zip([3.0, 3.0, 1.0]) {
            assert!((got - want).abs() < 1e-10, "{eig:?}");
        }
        let diag = symmetric_eigenvalues_3x3([[5.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 2.0]]);
        assert_eq!(diag, [5.0, 2.0, -1.0]);
    }

    #[test]
    fn cube_features_match_hand_values() {
        // 3x3x3 cube, unit spacing. Volume 27, exposed faces 6 * 9 = 54.
        // Per-axis population variance of 3 equally spaced points is (9 - 1) / 12 = 2/3.
        let mask = boxed((5, 5, 5), (1, 1, 1), (3, 3, 3));
        let f = shape_features_3d(&mask, (1.0, 1.0, 1.0)).unwrap();
        assert_eq!(f.volume, 27.0);
        assert_eq!(f.surface_area, 54.0);
        assert_eq!(f.bounding_box, ([1, 1, 1], [3, 3, 3]));
        let want = 4.0 * (2.0f64 / 3.0).sqrt();
        for axis in f.principal_axis_lengths {
            assert!((axis - want).abs() < 1e-9, "{axis} vs {want}");
        }
    }

    #[test]
    fn box_features_match_hand_values_with_spacing() {
        // 3x5x7 box. Volume 105. Surface 2 * (15 + 21 + 35) = 142.
        // Variances: (9-1)/12 = 2/3, (25-1)/12 = 2, (49-1)/12 = 4 along x, y, z.
        let mask = boxed((9, 9, 9), (0, 0, 0), (3, 5, 7));
        let f = shape_features_3d(&mask, (1.0, 1.0, 1.0)).unwrap();
        assert_eq!(f.volume, 105.0);
        assert_eq!(f.surface_area, 142.0);
        assert_eq!(f.bounding_box, ([0, 0, 0], [2, 4, 6]));
        let want = [4.0 * 2.0, 4.0 * 2.0f64.sqrt(), 4.0 * (2.0f64 / 3.0).sqrt()];
        for (got, w) in f.principal_axis_lengths.iter().zip(want) {
            assert!((got - w).abs() < 1e-9, "{got} vs {w}");
        }

        // Spacing (0.5, 2, 1): volume scales by 1, surface uses the face areas.
        let g = shape_features_3d(&mask, (0.5, 2.0, 1.0)).unwrap();
        assert!((g.volume - 105.0).abs() < 1e-9);
        // Faces per normal direction, each counted on both sides:
        // x-normal faces: 2 * (sy*sz) * (ny*nz) = 2 * 2 * 35 = 140
        // y-normal faces: 2 * (sx*sz) * (nx*nz) = 2 * 0.5 * 21 = 21
        // z-normal faces: 2 * (sx*sy) * (nx*ny) = 2 * 1 * 15 = 30
        assert!((g.surface_area - (140.0 + 21.0 + 30.0)).abs() < 1e-9, "{}", g.surface_area);
    }

    #[test]
    fn empty_and_bad_spacing_are_errors() {
        let empty = Image::new_3d(2, 2, 2, vec![false; 8]).unwrap();
        assert_eq!(shape_features_3d(&empty, (1.0, 1.0, 1.0)), Err(ShapeError::EmptyObject));
        let one = boxed((2, 2, 2), (0, 0, 0), (1, 1, 1));
        assert_eq!(shape_features_3d(&one, (0.0, 1.0, 1.0)), Err(ShapeError::InvalidSpacing));
    }
}
