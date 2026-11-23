//! 2D rank, mean and smoothing filters, plus morphological background subtraction.
//!
//! Windows are clipped at the image border: a pixel's median or mean uses only the
//! in-image samples of its window. Gaussian smoothing replicates the edge sample.
//! Inputs must be 2D (`depth == 1`).

use crate::image::{Image, ImageError};
use crate::morph::{StructuringElement, white_tophat};

fn require_2d(image: &Image<f64>) -> Result<(), ImageError> {
    if image.depth() == 1 {
        Ok(())
    } else {
        Err(ImageError::NotTwoDimensional { depth: image.depth() })
    }
}

/// Median over a `(2r+1)²` window, clipped at the image border.
///
/// # Errors
/// Returns [`ImageError::NotTwoDimensional`] for volumes, or [`ImageError::NonFinite`]
/// if the image contains NaN.
pub fn median_filter(image: &Image<f64>, radius: usize) -> Result<Image<f64>, ImageError> {
    require_2d(image)?;
    if image.data().iter().any(|v| v.is_nan()) {
        return Err(ImageError::NonFinite);
    }
    let (w, h) = (image.width(), image.height());
    let r = radius;
    let mut out = Vec::with_capacity(w * h);
    let mut window = Vec::with_capacity((2 * r + 1) * (2 * r + 1));
    for y in 0..h {
        for x in 0..w {
            window.clear();
            for yy in y.saturating_sub(r)..(y + r + 1).min(h) {
                for xx in x.saturating_sub(r)..(x + r + 1).min(w) {
                    window.push(image.data()[yy * w + xx]);
                }
            }
            window.sort_by(f64::total_cmp);
            let n = window.len();
            let m = if n % 2 == 1 {
                window[n / 2]
            } else {
                f64::midpoint(window[n / 2 - 1], window[n / 2])
            };
            out.push(m);
        }
    }
    Image::new_2d(w, h, out)
}

/// Mean over a `(2r+1)²` window, clipped at the image border.
///
/// # Errors
/// Returns [`ImageError::NotTwoDimensional`] for volumes.
pub fn mean_filter(image: &Image<f64>, radius: usize) -> Result<Image<f64>, ImageError> {
    require_2d(image)?;
    let (w, h) = (image.width(), image.height());
    let r = radius;
    let mut out = Vec::with_capacity(w * h);
    for y in 0..h {
        for x in 0..w {
            let (mut sum, mut n) = (0.0, 0usize);
            for yy in y.saturating_sub(r)..(y + r + 1).min(h) {
                for xx in x.saturating_sub(r)..(x + r + 1).min(w) {
                    sum += image.data()[yy * w + xx];
                    n += 1;
                }
            }
            out.push(sum / n as f64);
        }
    }
    Image::new_2d(w, h, out)
}

/// Gaussian smoothing with edge replication. The kernel radius is `ceil(3 * sigma)`.
///
/// Returns the input unchanged when `sigma` is not positive.
///
/// # Errors
/// Returns [`ImageError::NotTwoDimensional`] for volumes.
pub fn gaussian_smooth(image: &Image<f64>, sigma: f64) -> Result<Image<f64>, ImageError> {
    require_2d(image)?;
    if sigma.is_nan() || sigma <= 0.0 {
        return Ok(image.clone());
    }
    let radius = (3.0 * sigma).ceil() as isize;
    let mut kernel: Vec<f64> =
        (-radius..=radius).map(|i| (-(i * i) as f64 / (2.0 * sigma * sigma)).exp()).collect();
    let total: f64 = kernel.iter().sum();
    for k in &mut kernel {
        *k /= total;
    }
    let (w, h) = (image.width(), image.height());
    let clamp = |v: isize, n: usize| v.clamp(0, n as isize - 1) as usize;
    // Horizontal pass.
    let mut tmp = vec![0.0; w * h];
    for y in 0..h {
        for x in 0..w {
            tmp[y * w + x] = kernel
                .iter()
                .enumerate()
                .map(|(k, &c)| {
                    let sx = clamp(x as isize + k as isize - radius, w);
                    c * image.data()[y * w + sx]
                })
                .sum();
        }
    }
    // Vertical pass.
    let mut out = vec![0.0; w * h];
    for y in 0..h {
        for x in 0..w {
            out[y * w + x] = kernel
                .iter()
                .enumerate()
                .map(|(k, &c)| {
                    let sy = clamp(y as isize + k as isize - radius, h);
                    c * tmp[sy * w + x]
                })
                .sum();
        }
    }
    Image::new_2d(w, h, out)
}

/// Morphological background subtraction with a disk of the given radius.
///
/// The background is the grayscale opening of the image with a disk, so the result
/// is `image - opening(image)`: non-negative, and zero wherever the image is locally
/// smooth on the scale of the disk. This is the morphological equivalent of the
/// rolling-ball method (Sternberg, 1983); it is not bit-identical to `ImageJ`'s
/// rolling-ball implementation, which uses a sphere-based background estimate.
///
/// # Errors
/// Returns [`ImageError::NotTwoDimensional`] for volumes.
pub fn subtract_background(image: &Image<f64>, radius: usize) -> Result<Image<f64>, ImageError> {
    require_2d(image)?;
    Ok(white_tophat(image, &StructuringElement::disk(radius)))
}

#[allow(clippy::float_cmp)]
#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn median_removes_salt_noise() {
        let mut data = vec![1.0; 25];
        data[12] = 100.0;
        let img = Image::new_2d(5, 5, data).unwrap();
        let out = median_filter(&img, 1).unwrap();
        assert_eq!(*out.get(2, 2, 0).unwrap(), 1.0);
    }

    #[test]
    fn median_on_known_window() {
        // 3x1 row [1, 5, 2] with radius 1: pixel 1 sees {1,5,2} -> median 2.
        // Pixel 0 sees {1,5} (clipped) -> median 3. Pixel 2 sees {5,2} -> median 3.5.
        let img = Image::new_2d(3, 1, vec![1.0, 5.0, 2.0]).unwrap();
        let out = median_filter(&img, 1).unwrap();
        assert_eq!(out.data(), &[3.0, 2.0, 3.5]);
    }

    #[test]
    fn mean_on_clipped_window() {
        let img = Image::new_2d(3, 3, (1..=9).map(f64::from).collect()).unwrap();
        let out = mean_filter(&img, 1).unwrap();
        // Centre sees all nine values: mean 5. Corner (0,0) sees {1,2,4,5}: mean 3.
        assert_eq!(*out.get(1, 1, 0).unwrap(), 5.0);
        assert_eq!(*out.get(0, 0, 0).unwrap(), 3.0);
    }

    #[test]
    fn gaussian_preserves_constant_and_spreads_impulse() {
        let flat = Image::new_2d(7, 7, vec![4.0; 49]).unwrap();
        let out = gaussian_smooth(&flat, 1.0).unwrap();
        assert!(out.data().iter().all(|v| (v - 4.0).abs() < 1e-12));

        let mut data = vec![0.0; 49];
        data[24] = 1.0;
        let imp = Image::new_2d(7, 7, data).unwrap();
        let out = gaussian_smooth(&imp, 1.0).unwrap();
        let total: f64 = out.data().iter().sum();
        assert!((total - 1.0).abs() < 1e-9, "mass is preserved, got {total}");
        assert!(*out.get(3, 3, 0).unwrap() > *out.get(2, 3, 0).unwrap());
    }

    #[test]
    fn background_subtraction_keeps_small_bright_spot() {
        // Smooth ramp background plus a 1-pixel spot of height 10 on a 9x9 image.
        let mut data: Vec<f64> = (0..81).map(|i| (i % 9) as f64).collect();
        data[40] += 10.0;
        let img = Image::new_2d(9, 9, data).unwrap();
        let out = subtract_background(&img, 2).unwrap();
        assert!(*out.get(4, 4, 0).unwrap() >= 9.0);
        assert!(out.data().iter().all(|&v| v >= 0.0));
    }

    #[test]
    fn rejects_volumes() {
        let vol = Image::new_3d(2, 2, 2, vec![0.0; 8]).unwrap();
        assert!(median_filter(&vol, 1).is_err());
        assert!(mean_filter(&vol, 1).is_err());
        assert!(gaussian_smooth(&vol, 1.0).is_err());
        assert!(subtract_background(&vol, 1).is_err());
    }
}
