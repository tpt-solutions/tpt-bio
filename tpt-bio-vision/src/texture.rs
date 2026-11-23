//! Texture and scale features: grey-level co-occurrence (Haralick), granularity
//! spectrum, and Gabor filter responses.
//!
//! Everything here is 2D. Volumes are rejected with [`ImageError::NotTwoDimensional`]
//! where a 2D operation is needed.
//!
//! - [`quantize`] maps intensities to `levels` grey levels over the image's min-max range.
//! - [`haralick_features`] builds one symmetric grey-level co-occurrence matrix (GLCM)
//!   per object and direction, at distance 1, and averages the Haralick features over
//!   the four directions `0`, `45`, `90` and `135` degrees. Pairs are counted only when
//!   both pixels belong to the same object. A direction with no pairs in an object
//!   carries no information and is left out of the average; an object with no pairs at
//!   all gets zeros.
//!   Correlation is 0 when a marginal has zero variance (for example a constant patch),
//!   where the formula is undefined.
//! - [`granularity_spectrum`] opens the image with disks of radius `1..=max_radius` and
//!   reports the fraction of total intensity removed at each radius.
//! - [`gabor_kernel`] and [`gabor_bank`] build real-valued Gabor kernels with the DC
//!   component removed, [`gabor_response`] convolves with them (zero padding), and
//!   [`mean_abs_response`] averages the absolute response per object.

use std::f64::consts::PI;

use crate::features::measure;
use crate::image::{Image, ImageError};
use crate::morph::{StructuringElement, open_gray};

/// Haralick texture features of one co-occurrence matrix, or their average over
/// directions.
#[derive(Debug, Clone, Copy, Default, PartialEq)]
pub struct Haralick {
    /// Angular second moment (energy): `sum P(i,j)^2`.
    pub angular_second_moment: f64,
    /// Contrast: `sum (i-j)^2 P(i,j)`.
    pub contrast: f64,
    /// Correlation: `sum (i-mu)(j-mu) P(i,j) / var`, or 0 when `var` is zero.
    pub correlation: f64,
    /// Homogeneity (inverse difference moment): `sum P(i,j) / (1 + (i-j)^2)`.
    pub homogeneity: f64,
    /// Entropy: `-sum P(i,j) ln P(i,j)`, in nats.
    pub entropy: f64,
}

/// Haralick features of one 2D object, averaged over the four directions.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ObjectHaralick {
    /// Consecutive object id, `1..=count`.
    pub label: u32,
    /// The averaged features.
    pub features: Haralick,
}

/// The four directions at distance 1, as `(dx, dy)` with `y` pointing down:
/// 0 degrees, 45, 90 and 135.
pub const DIRECTIONS: [(isize, isize); 4] = [(1, 0), (1, -1), (0, -1), (-1, -1)];

/// Quantises intensities to `levels` grey levels `0..levels` over the image's min-max
/// range. Each value maps to `floor((v - min) / (max - min) * levels)`, clamped to
/// `levels - 1`. A constant image maps to level 0 everywhere.
///
/// # Errors
/// Returns [`ImageError::NonFinite`] if a sample is NaN.
pub fn quantize(image: &Image<f64>, levels: usize) -> Result<Image<u32>, ImageError> {
    let levels = levels.max(1);
    let values = image.data();
    if values.iter().any(|v| v.is_nan()) {
        return Err(ImageError::NonFinite);
    }
    let lo = values.iter().copied().fold(f64::INFINITY, f64::min);
    let hi = values.iter().copied().fold(f64::NEG_INFINITY, f64::max);
    let top = (levels - 1) as u32;
    let data = if hi > lo {
        let scale = levels as f64 / (hi - lo);
        values.iter().map(|&v| (((v - lo) * scale) as u32).min(top)).collect()
    } else {
        vec![0; values.len()]
    };
    Image::new_3d(image.width(), image.height(), image.depth(), data)
}

/// Co-occurrence counts for the pixels `px` of one object, at offset `d`.
///
/// The matrix is symmetric: each pair `(p, n)` adds one to `(i, j)` and one to `(j, i)`,
/// where `i` and `j` are the grey levels of `p` and `n`. Levels at or above `levels` are
/// clamped to `levels - 1`. `raw` holds the original object labels, and a pair counts
/// only when both pixels carry the same label.
fn pair_counts(
    grey: &Image<u32>,
    raw: &[u32],
    px: &[usize],
    levels: usize,
    d: (isize, isize),
) -> Vec<u64> {
    let mut counts = vec![0u64; levels * levels];
    let top = levels - 1;
    let values = grey.data();
    for &p in px {
        let Some(n) = grey.offset(p, (d.0, d.1, 0)) else {
            continue;
        };
        if raw[n] != raw[p] {
            continue;
        }
        let i = (values[p] as usize).min(top);
        let j = (values[n] as usize).min(top);
        counts[i * levels + j] += 1;
        counts[j * levels + i] += 1;
    }
    counts
}

/// Haralick features of one co-occurrence matrix given as `levels * levels` counts,
/// row-major, with `counts[i * levels + j]` the count of grey-level pair `(i, j)`.
///
/// All-zero counts give all-zero features. The matrix is normalised to sum 1 first.
#[must_use]
pub fn haralick_from_glcm(counts: &[u64], levels: usize) -> Haralick {
    let total: u64 = counts.iter().sum();
    if total == 0 {
        return Haralick::default();
    }
    let n = total as f64;
    let mut out = Haralick::default();
    let mut marginal = vec![0.0f64; levels];
    for i in 0..levels {
        for j in 0..levels {
            let c = counts[i * levels + j];
            if c == 0 {
                continue;
            }
            let p = c as f64 / n;
            let d = i as f64 - j as f64;
            out.angular_second_moment += p * p;
            out.contrast += d * d * p;
            out.homogeneity += p / (1.0 + d * d);
            out.entropy -= p * p.ln();
            marginal[i] += p;
        }
    }
    // The GLCM is symmetric, so the column marginal equals the row marginal.
    let mu: f64 = marginal.iter().enumerate().map(|(i, &p)| i as f64 * p).sum();
    let var: f64 = marginal.iter().enumerate().map(|(i, &p)| (i as f64 - mu).powi(2) * p).sum();
    if var > 1e-12 {
        let mut cov = 0.0;
        for i in 0..levels {
            for j in 0..levels {
                let c = counts[i * levels + j];
                if c > 0 {
                    let p = c as f64 / n;
                    cov += (i as f64 - mu) * (j as f64 - mu) * p;
                }
            }
        }
        out.correlation = cov / var;
    }
    out
}

/// The symmetric co-occurrence counts of one object (the pixels with label `object`) at
/// offset `d`, as `levels * levels` row-major counts.
///
/// # Errors
/// Returns a shape error if `grey` and `labels` differ in shape, and
/// [`ImageError::NotTwoDimensional`] for volumes.
pub fn glcm_counts(
    grey: &Image<u32>,
    labels: &Image<u32>,
    object: u32,
    levels: usize,
    d: (isize, isize),
) -> Result<Vec<u64>, ImageError> {
    labels.check_same_shape(grey)?;
    require_2d(labels)?;
    let px: Vec<usize> =
        labels.data().iter().enumerate().filter_map(|(p, &l)| (l == object).then_some(p)).collect();
    Ok(pair_counts(grey, labels.data(), &px, levels.max(1), d))
}

/// Haralick features of every object, from grey levels in `grey` (values below `levels`,
/// as from [`quantize`]).
///
/// # Errors
/// Returns [`ImageError::ShapeMismatch`] if `grey` and `labels` differ in shape, and
/// [`ImageError::NotTwoDimensional`] for volumes.
pub fn haralick_features(
    grey: &Image<u32>,
    labels: &Image<u32>,
    levels: usize,
) -> Result<Vec<ObjectHaralick>, ImageError> {
    labels.check_same_shape(grey)?;
    require_2d(labels)?;
    let levels = levels.max(1);
    let raw = labels.data();
    Ok(measure(labels, |label, px| {
        let mut sum = Haralick::default();
        let mut used = 0usize;
        for &d in &DIRECTIONS {
            let counts = pair_counts(grey, raw, px, levels, d);
            if counts.iter().all(|&c| c == 0) {
                continue;
            }
            let h = haralick_from_glcm(&counts, levels);
            sum.angular_second_moment += h.angular_second_moment;
            sum.contrast += h.contrast;
            sum.correlation += h.correlation;
            sum.homogeneity += h.homogeneity;
            sum.entropy += h.entropy;
            used += 1;
        }
        if used > 0 {
            let k = used as f64;
            sum.angular_second_moment /= k;
            sum.contrast /= k;
            sum.correlation /= k;
            sum.homogeneity /= k;
            sum.entropy /= k;
        }
        ObjectHaralick { label, features: sum }
    }))
}

/// Granularity spectrum of a 2D intensity image.
///
/// Element `i` of the result is the fraction of the image's total intensity removed by
/// opening with a disk of radius `i + 1` and not by the disk of radius `i`:
/// `(S(i) - S(i + 1)) / total`, where `S(r)` is the sum of the image opened with
/// `disk(r)`, and `S(0)` is the total. Values are clamped at 0. An image with zero total
/// intensity gives zeros.
///
/// Cost: one grayscale opening per radius, each `O(N * disk area)`.
#[must_use]
pub fn granularity_spectrum(image: &Image<f64>, max_radius: usize) -> Vec<f64> {
    let total: f64 = image.data().iter().sum();
    let mut spectrum = Vec::with_capacity(max_radius);
    let mut previous = total;
    for r in 1..=max_radius {
        let opened = open_gray(image, &StructuringElement::disk(r));
        let remaining: f64 = opened.data().iter().sum();
        let removed = if total > 0.0 { ((previous - remaining) / total).max(0.0) } else { 0.0 };
        spectrum.push(removed);
        previous = remaining;
    }
    spectrum
}

/// A real-valued 2D Gabor kernel with its DC component removed.
///
/// The kernel is square with side `2 * half_width() + 1`, stored row-major.
#[derive(Debug, Clone, PartialEq)]
pub struct GaborKernel {
    half: usize,
    weights: Vec<f64>,
}

impl GaborKernel {
    /// Half the side length: the kernel covers `-half..=half` in each axis.
    #[must_use]
    pub fn half_width(&self) -> usize {
        self.half
    }

    /// Side length `2 * half_width() + 1`.
    #[must_use]
    pub fn size(&self) -> usize {
        2 * self.half + 1
    }

    /// Kernel weights, row-major.
    #[must_use]
    pub fn weights(&self) -> &[f64] {
        &self.weights
    }
}

/// Builds a Gabor kernel: a Gaussian envelope times a cosine carrier.
///
/// With `x' = x cos(theta) + y sin(theta)` and `y' = -x sin(theta) + y cos(theta)`, the
/// weight is `exp(-(x'^2 + aspect^2 y'^2) / (2 sigma^2)) * cos(2 pi x' / wavelength)`.
/// The carrier runs along `x'`, so a stripe pattern matches when `theta` points across
/// its stripes. The kernel extends to `ceil(3 sigma)` pixels from its centre, and its
/// mean is subtracted so that flat regions give no response.
///
/// `wavelength`, `sigma` and `aspect` must be positive; `sigma` is in pixels.
#[must_use]
pub fn gabor_kernel(wavelength: f64, theta: f64, sigma: f64, aspect: f64) -> GaborKernel {
    let half = (3.0 * sigma).ceil().max(1.0) as usize;
    let (sin_t, cos_t) = theta.sin_cos();
    let h = half as isize;
    let mut weights = Vec::with_capacity((2 * half + 1) * (2 * half + 1));
    for dy in -h..=h {
        for dx in -h..=h {
            let (x, y) = (dx as f64, dy as f64);
            let xr = x * cos_t + y * sin_t;
            let yr = -x * sin_t + y * cos_t;
            let envelope = (-(xr * xr + aspect * aspect * yr * yr) / (2.0 * sigma * sigma)).exp();
            weights.push(envelope * (2.0 * PI * xr / wavelength).cos());
        }
    }
    let mean = weights.iter().sum::<f64>() / weights.len() as f64;
    for w in &mut weights {
        *w -= mean;
    }
    GaborKernel { half, weights }
}

/// One Gabor kernel per orientation (radians), all with the same wavelength, sigma and
/// aspect ratio.
#[must_use]
pub fn gabor_bank(wavelength: f64, sigma: f64, aspect: f64, thetas: &[f64]) -> Vec<GaborKernel> {
    thetas.iter().map(|&t| gabor_kernel(wavelength, t, sigma, aspect)).collect()
}

/// Convolves a 2D image with a Gabor kernel, treating pixels outside the image as 0.
///
/// Cost: `O(N * kernel size)`. The kernel is even, so correlation and convolution agree.
///
/// # Errors
/// Returns [`ImageError::NotTwoDimensional`] for volumes.
pub fn gabor_response(image: &Image<f64>, kernel: &GaborKernel) -> Result<Image<f64>, ImageError> {
    require_2d(image)?;
    let (w, h) = (image.width(), image.height());
    let size = kernel.size();
    let half = kernel.half_width() as isize;
    let data = image.data();
    let mut out = vec![0.0f64; w * h];
    for y in 0..h {
        for x in 0..w {
            let mut acc = 0.0;
            for (ky, row) in kernel.weights().chunks(size).enumerate() {
                let sy = y as isize + ky as isize - half;
                if sy < 0 || sy >= h as isize {
                    continue;
                }
                for (kx, &weight) in row.iter().enumerate() {
                    let sx = x as isize + kx as isize - half;
                    if sx < 0 || sx >= w as isize {
                        continue;
                    }
                    acc += weight * data[sy as usize * w + sx as usize];
                }
            }
            out[y * w + x] = acc;
        }
    }
    Image::new_2d(w, h, out)
}

/// Mean of the absolute response over each object's pixels, one value per object in
/// consecutive id order (see [`crate::features::measure`]).
///
/// # Errors
/// Returns [`ImageError::ShapeMismatch`] if `response` and `labels` differ in shape.
pub fn mean_abs_response(
    response: &Image<f64>,
    labels: &Image<u32>,
) -> Result<Vec<f64>, ImageError> {
    labels.check_same_shape(response)?;
    let values = response.data();
    Ok(measure(labels, |_, px| {
        if px.is_empty() {
            0.0
        } else {
            px.iter().map(|&p| values[p].abs()).sum::<f64>() / px.len() as f64
        }
    }))
}

fn require_2d<T>(image: &Image<T>) -> Result<(), ImageError> {
    if image.depth() == 1 {
        Ok(())
    } else {
        Err(ImageError::NotTwoDimensional { depth: image.depth() })
    }
}

#[allow(clippy::float_cmp)]
#[cfg(test)]
mod tests {
    use super::*;

    fn approx(a: f64, b: f64) {
        assert!((a - b).abs() < 1e-9, "{a} != {b}");
    }

    #[test]
    fn checkerboard_glcm_counts_by_hand() {
        // Grey levels: [[0, 1], [1, 0]] with one object covering all four pixels.
        //   0 degrees:   pairs (0,1) in row 0 and (1,0) in row 1 -> counts (0,1)=2, (1,0)=2.
        //   45 degrees:  only (0,1) -> (1,0) is in range, both level 1 -> counts (1,1)=2.
        //   90 degrees:  pairs (1,0) at x=0 and (0,1) at x=1 -> counts (0,1)=2, (1,0)=2.
        //   135 degrees: only (1,1) -> (0,0) is in range, both level 0 -> counts (0,0)=2.
        let grey = Image::new_2d(2, 2, vec![0u32, 1, 1, 0]).unwrap();
        let labels = Image::new_2d(2, 2, vec![1u32; 4]).unwrap();
        let per_direction = [[0u64, 2, 2, 0], [0, 0, 0, 2], [0, 2, 2, 0], [2, 0, 0, 0]];
        for (d, expected) in DIRECTIONS.iter().zip(per_direction) {
            let counts = glcm_counts(&grey, &labels, 1, 2, *d).unwrap();
            assert_eq!(counts, expected, "direction {d:?}");
        }
    }

    #[test]
    fn checkerboard_haralick_hand_values() {
        // 0 and 90 degrees: P = [[0, .5], [.5, 0]]. ASM = .5, contrast = 1,
        //   homogeneity = .5 / 2 * 2 = .5, entropy = ln 2, correlation = -1
        //   (mu = .5, var = .25, cov = -.25).
        // 45 and 135 degrees: P = [[0, 0], [0, 1]] or [[1, 0], [0, 0]]. ASM = 1,
        //   contrast = 0, homogeneity = 1, entropy = 0, correlation = 0 (zero variance).
        // Averages over four directions: ASM .75, contrast .5, homogeneity .75,
        // entropy ln2 / 2, correlation -.5.
        let grey = Image::new_2d(2, 2, vec![0u32, 1, 1, 0]).unwrap();
        let labels = Image::new_2d(2, 2, vec![1u32; 4]).unwrap();
        let f = haralick_features(&grey, &labels, 2).unwrap();
        assert_eq!(f.len(), 1);
        let h = f[0].features;
        approx(h.angular_second_moment, 0.75);
        approx(h.contrast, 0.5);
        approx(h.homogeneity, 0.75);
        approx(h.entropy, 2f64.ln() / 2.0);
        approx(h.correlation, -0.5);
    }

    #[test]
    fn constant_patch_has_zero_contrast_and_unit_energy() {
        let grey = Image::new_2d(3, 3, vec![2u32; 9]).unwrap();
        let labels = Image::new_2d(3, 3, vec![1u32; 9]).unwrap();
        let h = haralick_features(&grey, &labels, 4).unwrap()[0].features;
        approx(h.contrast, 0.0);
        approx(h.angular_second_moment, 1.0);
        approx(h.homogeneity, 1.0);
        approx(h.entropy, 0.0);
        approx(h.correlation, 0.0);
    }

    #[test]
    fn pairs_do_not_cross_objects() {
        // Two objects side by side: the level-0 pixel of object 1 must not pair with
        // the level-1 pixel of object 2 across the boundary.
        let grey = Image::new_2d(2, 1, vec![0u32, 1]).unwrap();
        let labels = Image::new_2d(2, 1, vec![1u32, 2]).unwrap();
        let f = haralick_features(&grey, &labels, 2).unwrap();
        assert_eq!(f.len(), 2);
        approx(f[0].features.contrast, 0.0);
        approx(f[1].features.contrast, 0.0);
    }

    #[test]
    fn quantize_maps_range_to_levels() {
        let img = Image::new_2d(4, 1, vec![0.0, 1.0, 2.0, 4.0]).unwrap();
        let q = quantize(&img, 4).unwrap();
        // scale = 4 / 4 = 1: 0 -> 0, 1 -> 1, 2 -> 2, 4 -> 4 clamped to 3.
        assert_eq!(q.data(), &[0, 1, 2, 3]);
        let flat = Image::new_2d(2, 1, vec![5.0, 5.0]).unwrap();
        assert_eq!(quantize(&flat, 4).unwrap().data(), &[0, 0]);
    }

    #[test]
    fn granularity_of_a_single_spot() {
        // A disk of radius 3 with value 1 on a zero background. Opening with disk(r)
        // keeps a disk of radius r exactly, so removal is 0 for r <= 3, and disk(4) is
        // larger than the spot, so it removes everything.
        let r = 3usize;
        let size = 11usize;
        let c = 5isize;
        let mut data = vec![0.0; size * size];
        for y in 0..size {
            for x in 0..size {
                let (dx, dy) = (x as isize - c, y as isize - c);
                if dx * dx + dy * dy <= (r * r) as isize {
                    data[y * size + x] = 1.0;
                }
            }
        }
        let img = Image::new_2d(size, size, data).unwrap();
        let spec = granularity_spectrum(&img, 5);
        // spec[i] is the removal at radius i + 1. Radius r - 1 and radius r keep the spot.
        // Radius 1 is not checked: disk(1) is a plus shape, which cannot cover the
        // diagonal corners of a discrete disk, so it removes about 14% at this size.
        assert!(spec[r - 2].abs() < 1e-12, "radius {} removed {}", r - 1, spec[r - 2]);
        assert!(spec[r - 1].abs() < 1e-12, "radius {r} should keep the spot");
        assert!(spec[r] > 0.99, "radius {} removed only {}", r + 1, spec[r]);
    }

    #[test]
    fn gabor_kernel_is_zero_mean_and_matching_orientation_responds_more() {
        let k = gabor_kernel(4.0, 0.0, 2.0, 1.0);
        let sum: f64 = k.weights().iter().sum();
        assert!(sum.abs() < 1e-9);
        assert_eq!(k.size(), 2 * k.half_width() + 1);

        // Vertical stripes: the intensity varies along x with period 4.
        let (w, h) = (16usize, 16usize);
        let data: Vec<f64> =
            (0..w * h).map(|i| 1.0 + (2.0 * PI * (i % w) as f64 / 4.0).cos()).collect();
        let img = Image::new_2d(w, h, data).unwrap();
        let labels = Image::new_2d(w, h, vec![1u32; w * h]).unwrap();

        let matched = gabor_response(&img, &gabor_kernel(4.0, 0.0, 2.0, 1.0)).unwrap();
        let orthogonal = gabor_response(&img, &gabor_kernel(4.0, PI / 2.0, 2.0, 1.0)).unwrap();
        let m = mean_abs_response(&matched, &labels).unwrap()[0];
        let o = mean_abs_response(&orthogonal, &labels).unwrap()[0];
        assert!(m > 5.0 * o, "matched {m} should dominate orthogonal {o}");
    }

    #[test]
    fn gabor_bank_has_one_kernel_per_orientation() {
        let bank = gabor_bank(4.0, 2.0, 0.5, &[0.0, PI / 4.0, PI / 2.0, 3.0 * PI / 4.0]);
        assert_eq!(bank.len(), 4);
        assert_ne!(bank[0], bank[1]);
    }

    #[test]
    fn volumes_are_rejected() {
        let vol = Image::new_3d(2, 2, 2, vec![0u32; 8]).unwrap();
        let labels = Image::new_3d(2, 2, 2, vec![1u32; 8]).unwrap();
        assert!(haralick_features(&vol, &labels, 2).is_err());
        let fimg = Image::new_3d(2, 2, 2, vec![0.0f64; 8]).unwrap();
        assert!(gabor_response(&fimg, &gabor_kernel(4.0, 0.0, 1.0, 1.0)).is_err());
    }
}
