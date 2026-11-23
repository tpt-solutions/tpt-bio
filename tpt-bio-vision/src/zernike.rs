//! Zernike moments of a 2D binary object.
//!
//! Object pixels are mapped onto the unit disk using the object's centroid `(cx, cy)`
//! and its equivalent radius `R = sqrt(area / pi)`, where `area` is the pixel count:
//! `u = (x - cx) / R`, `v = (y - cy) / R`, `rho = sqrt(u^2 + v^2)`, `theta = atan2(v, u)`.
//! Pixels with `rho > 1` are excluded. Each pixel has area `1 / R^2` on the disk.
//!
//! The Zernike moment of order `n` and repetition `m` (`m >= 0`, `n - m` even) is
//! `Z_nm = (n + 1) / pi * sum_p R_nm(rho_p) * exp(-i m theta_p) * (1 / R^2)`,
//! where `R_nm` is the Zernike radial polynomial.

use std::f64::consts::PI;

use crate::image::{Image, ImageError};

/// Errors from Zernike moment computation.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum ZernikeError {
    /// The mask is a volume; Zernike moments are defined for 2D objects.
    NotTwoDimensional {
        /// Depth of the offending mask.
        depth: usize,
    },
    /// The mask has no foreground pixels, so there is no centroid.
    EmptyObject,
    /// An underlying image error.
    Image(ImageError),
}

impl std::fmt::Display for ZernikeError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::NotTwoDimensional { depth } => {
                write!(f, "Zernike moments need a 2D mask, got depth {depth}")
            }
            Self::EmptyObject => f.write_str("mask has no foreground pixels"),
            Self::Image(e) => write!(f, "{e}"),
        }
    }
}

impl std::error::Error for ZernikeError {}

/// One complex Zernike moment `Z_nm`.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ZernikeMoment {
    /// Order `n`.
    pub n: u32,
    /// Repetition `m`, with `0 <= m <= n` and `n - m` even.
    pub m: u32,
    /// Real part.
    pub re: f64,
    /// Imaginary part.
    pub im: f64,
}

impl ZernikeMoment {
    /// `|Z_nm|`, which is invariant under rotation of the object.
    #[must_use]
    pub fn magnitude(&self) -> f64 {
        self.re.hypot(self.im)
    }
}

/// The Zernike radial polynomial `R_nm(rho)`, for `0 <= m <= n` and `n - m` even.
#[must_use]
pub fn radial_polynomial(n: u32, m: u32, rho: f64) -> f64 {
    if m > n || !(n - m).is_multiple_of(2) {
        return 0.0;
    }
    let half_sum = u32::midpoint(n, m);
    let half_diff = (n - m) / 2;
    let mut total = 0.0;
    for s in 0..=half_diff {
        let num = factorial(n - s);
        let den = factorial(s) * factorial(half_sum - s) * factorial(half_diff - s);
        let sign = if s % 2 == 0 { 1.0 } else { -1.0 };
        total += sign * num / den * rho.powi((n - 2 * s) as i32);
    }
    total
}

fn factorial(k: u32) -> f64 {
    (1..=k).fold(1.0, |acc, i| acc * f64::from(i))
}

/// Zernike moments `Z_nm` for every `0 <= m <= n <= order` with `n - m` even.
///
/// The moments are returned in order of `n`, then `m`.
///
/// # Errors
/// Returns [`ZernikeError::NotTwoDimensional`] for a volume and
/// [`ZernikeError::EmptyObject`] for an empty mask.
pub fn zernike_moments(mask: &Image<bool>, order: u32) -> Result<Vec<ZernikeMoment>, ZernikeError> {
    if mask.depth() != 1 {
        return Err(ZernikeError::NotTwoDimensional { depth: mask.depth() });
    }
    let w = mask.width();
    let pixels: Vec<(f64, f64)> = mask
        .data()
        .iter()
        .enumerate()
        .filter(|&(_, &on)| on)
        .map(|(i, _)| ((i % w) as f64, (i / w) as f64))
        .collect();
    let area = pixels.len();
    if area == 0 {
        return Err(ZernikeError::EmptyObject);
    }
    let count = area as f64;
    let cx = pixels.iter().map(|p| p.0).sum::<f64>() / count;
    let cy = pixels.iter().map(|p| p.1).sum::<f64>() / count;
    let radius = (count / PI).sqrt();
    let pixel_area = 1.0 / (radius * radius);

    // Polar coordinates of the object pixels inside the unit disk.
    let polar: Vec<(f64, f64)> = pixels
        .iter()
        .filter_map(|&(x, y)| {
            let u = (x - cx) / radius;
            let v = (y - cy) / radius;
            let rho = u.hypot(v);
            (rho <= 1.0).then(|| (rho, v.atan2(u)))
        })
        .collect();

    let mut out = Vec::new();
    for n in 0..=order {
        let norm = f64::from(n + 1) / PI * pixel_area;
        for m in (n % 2..=n).step_by(2) {
            let (mut re, mut im) = (0.0, 0.0);
            for &(rho, theta) in &polar {
                let r = radial_polynomial(n, m, rho);
                let phase = f64::from(m) * theta;
                re += r * phase.cos();
                im -= r * phase.sin();
            }
            out.push(ZernikeMoment { n, m, re: norm * re, im: norm * im });
        }
    }
    Ok(out)
}

#[allow(clippy::float_cmp)]
#[cfg(test)]
mod tests {
    use super::*;

    /// A binary disk of the given radius centred on the pixel `(c, c)`.
    fn disk(size: usize, c: f64, r: f64) -> Image<bool> {
        let data: Vec<bool> = (0..size * size)
            .map(|i| {
                let (x, y) = ((i % size) as f64, (i / size) as f64);
                (x - c).powi(2) + (y - c).powi(2) <= r * r
            })
            .collect();
        Image::new_2d(size, size, data).unwrap()
    }

    /// Asymmetric L-shaped object, to test rotation.
    fn l_shape(size: usize) -> Image<bool> {
        let data: Vec<bool> = (0..size * size)
            .map(|i| {
                let (x, y) = (i % size, i / size);
                (x < 4 && y < 12) || (y < 4 && x < 12)
            })
            .collect();
        Image::new_2d(size, size, data).unwrap()
    }

    /// Rotates a square mask by 90 degrees clockwise.
    fn rotate_90(img: &Image<bool>) -> Image<bool> {
        let s = img.width();
        let mut out = vec![false; s * s];
        for y in 0..s {
            for x in 0..s {
                out[x * s + (s - 1 - y)] = img.data()[y * s + x];
            }
        }
        Image::new_2d(s, s, out).unwrap()
    }

    #[test]
    fn radial_polynomial_known_values() {
        // R_20(rho) = 2 rho^2 - 1, R_40(rho) = 6 rho^4 - 6 rho^2 + 1, R_11(rho) = rho.
        assert!((radial_polynomial(2, 0, 0.5) - (2.0 * 0.25 - 1.0)).abs() < 1e-12);
        let r = 0.7f64;
        assert!((radial_polynomial(4, 0, r) - (6.0 * r.powi(4) - 6.0 * r * r + 1.0)).abs() < 1e-12);
        assert!((radial_polynomial(1, 1, r) - r).abs() < 1e-12);
        assert_eq!(radial_polynomial(2, 1, 0.5), 0.0);
    }

    #[test]
    fn disk_has_only_zeroth_moment_up_to_discretisation() {
        // Centre on a pixel centre so the discrete disk is symmetric about the centroid.
        let mask = disk(181, 90.0, 60.0);
        let moments = zernike_moments(&mask, 6).unwrap();
        let z00 = moments.iter().find(|z| z.n == 0 && z.m == 0).unwrap();
        assert!((z00.re - 1.0).abs() < 0.05, "Z00 = {}", z00.re);
        assert!(z00.im.abs() < 1e-9);
        for z in moments.iter().filter(|z| z.n > 0) {
            assert!(z.magnitude() < 0.05, "Z{}{} = {:?}", z.n, z.m, z);
        }
    }

    #[test]
    fn magnitudes_are_invariant_under_rotation() {
        let mask = l_shape(16);
        let rotated = rotate_90(&mask);
        let a = zernike_moments(&mask, 8).unwrap();
        let b = zernike_moments(&rotated, 8).unwrap();
        assert_eq!(a.len(), b.len());
        for (za, zb) in a.iter().zip(&b) {
            assert_eq!((za.n, za.m), (zb.n, zb.m));
            assert!(
                (za.magnitude() - zb.magnitude()).abs() < 1e-9,
                "Z{}{}: {} vs {}",
                za.n,
                za.m,
                za.magnitude(),
                zb.magnitude()
            );
        }
        // The shape is not symmetric, so some higher-order magnitude is clearly nonzero.
        assert!(a.iter().any(|z| z.n >= 2 && z.magnitude() > 0.01));
    }

    #[test]
    fn rejects_volumes_and_empty_masks() {
        let vol = Image::new_3d(2, 2, 2, vec![true; 8]).unwrap();
        assert_eq!(zernike_moments(&vol, 2), Err(ZernikeError::NotTwoDimensional { depth: 2 }));
        let empty = Image::new_2d(3, 3, vec![false; 9]).unwrap();
        assert_eq!(zernike_moments(&empty, 2), Err(ZernikeError::EmptyObject));
    }
}
