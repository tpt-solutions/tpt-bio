//! Image filters built on `tpt-simd` (ADR-002).
//!
//! Images are row-major `f32` planes. Convolution keeps the input size and uses zero
//! padding, with each kernel centred at `len / 2`.

use std::fmt;

use crate::image::{Image, ImageError};

/// Errors from filter calls with inconsistent inputs.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum FilterError {
    /// `data.len()` does not equal `width * height`.
    SizeMismatch {
        /// Expected number of samples.
        expected: usize,
        /// Number of samples supplied.
        actual: usize,
    },
    /// A kernel has no taps.
    EmptyKernel,
    /// A square kernel must have an odd, positive side length.
    InvalidKernelSize {
        /// The side length supplied.
        size: usize,
    },
    /// A square kernel's length is not `size * size`.
    KernelLength {
        /// Expected number of taps, `size * size`.
        expected: usize,
        /// Number of taps supplied.
        actual: usize,
    },
    /// An operation that is defined only for 2D images got a volume.
    NotTwoDimensional {
        /// Depth of the offending image.
        depth: usize,
    },
    /// Building the output image failed.
    Image(ImageError),
}

impl fmt::Display for FilterError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::SizeMismatch { expected, actual } => {
                write!(f, "image has {actual} samples, expected {expected}")
            }
            Self::EmptyKernel => f.write_str("convolution kernel is empty"),
            Self::InvalidKernelSize { size } => {
                write!(f, "square kernel side must be odd and positive, got {size}")
            }
            Self::KernelLength { expected, actual } => {
                write!(f, "square kernel has {actual} taps, expected {expected}")
            }
            Self::NotTwoDimensional { depth } => {
                write!(f, "operation needs a 2D image, got depth {depth}")
            }
            Self::Image(e) => write!(f, "{e}"),
        }
    }
}

impl std::error::Error for FilterError {}

impl From<ImageError> for FilterError {
    fn from(e: ImageError) -> Self {
        Self::Image(e)
    }
}

/// How samples outside the image are defined along each axis.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Boundary {
    /// Outside samples are zero.
    Zero,
    /// Outside samples repeat the nearest edge sample.
    Clamp,
    /// Outside samples reflect about the edge sample, which is not repeated:
    /// `-1` reads sample `1`, and `n` reads sample `n - 2`.
    Mirror,
}

impl Boundary {
    /// Source index for sample `i` along an axis of length `n`, or `None` for a zero sample.
    /// Returns `None` for every `i` when `n == 0`.
    pub(crate) fn index(self, i: isize, n: usize) -> Option<usize> {
        if n == 0 {
            return None;
        }
        let len = n as isize;
        if (0..len).contains(&i) {
            return Some(i as usize);
        }
        match self {
            Self::Zero => None,
            Self::Clamp => Some(i.clamp(0, len - 1) as usize),
            Self::Mirror => {
                if n == 1 {
                    return Some(0);
                }
                let period = 2 * (len - 1);
                let m = i.rem_euclid(period);
                Some(if m < len { m } else { period - m } as usize)
            }
        }
    }
}

/// Reusable buffers for [`convolve_line`].
#[derive(Debug, Default)]
pub(crate) struct LineScratch {
    padded: Vec<f32>,
    full: Vec<f32>,
}

/// Convolves one line of samples with a 1D kernel, with boundary handling.
///
/// Writes (or, when `accumulate`, adds) `out[i] = sum_j kernel[j] * line[i + c - j]`
/// with `c = kernel.len() / 2`, where out-of-range samples come from `boundary`.
/// The padded line is convolved with `tpt-simd`'s 1D kernel.
///
/// `out` must have the same length as `line`.
pub(crate) fn convolve_line(
    line: &[f32],
    kernel: &[f32],
    boundary: Boundary,
    scratch: &mut LineScratch,
    out: &mut [f32],
    accumulate: bool,
) {
    let n = line.len();
    let k = kernel.len();
    debug_assert!(k > 0 && out.len() == n, "convolve_line: bad lengths");
    if n == 0 || k == 0 {
        return;
    }
    let pad = k - 1;
    let c = k / 2;
    scratch.padded.clear();
    scratch.padded.extend(
        (0..n + 2 * pad)
            .map(|m| boundary.index(m as isize - pad as isize, n).map_or(0.0, |s| line[s])),
    );
    scratch.full.resize(scratch.padded.len() + pad, 0.0);
    tpt_simd::convolve::convolve_1d_f32(&scratch.padded, kernel, &mut scratch.full);
    for (i, o) in out.iter_mut().enumerate() {
        let v = scratch.full[i + c + pad];
        if accumulate {
            *o += v;
        } else {
            *o = v;
        }
    }
}

/// General (non-separable) 2D convolution with a square `size x size` kernel.
///
/// `kernel` is row-major with `kernel[j * size + i]` the tap at offset `(i, j)` from the
/// centre `(size / 2, size / 2)`, so
/// `out[y][x] = sum_j sum_i kernel[j * size + i] * in[y + c - j][x + c - i]`.
/// Each kernel row is applied with a SIMD 1D convolution (`tpt-simd`), so the cost is
/// `size` 1D passes per row rather than `size * size` scalar multiplies.
///
/// # Errors
/// Returns [`FilterError::NotTwoDimensional`] for volumes,
/// [`FilterError::InvalidKernelSize`] when `size` is zero or even, and
/// [`FilterError::KernelLength`] when `kernel.len() != size * size`.
pub fn convolve_2d(
    input: &Image<f32>,
    kernel: &[f32],
    size: usize,
    boundary: Boundary,
) -> Result<Image<f32>, FilterError> {
    if input.depth() != 1 {
        return Err(FilterError::NotTwoDimensional { depth: input.depth() });
    }
    if size == 0 || size.is_multiple_of(2) {
        return Err(FilterError::InvalidKernelSize { size });
    }
    if kernel.len() != size * size {
        return Err(FilterError::KernelLength { expected: size * size, actual: kernel.len() });
    }
    let (w, h) = (input.width(), input.height());
    let c = (size / 2) as isize;
    let mut out = vec![0.0f32; w * h];
    let mut scratch = LineScratch::default();
    for y in 0..h {
        let out_row = &mut out[y * w..(y + 1) * w];
        for j in 0..size {
            let sy = y as isize + c - j as isize;
            // A zero boundary row contributes nothing; clamp and mirror map to a real row.
            let Some(sy) = boundary.index(sy, h) else {
                continue;
            };
            let row = &input.data()[sy * w..(sy + 1) * w];
            convolve_line(
                row,
                &kernel[j * size..(j + 1) * size],
                boundary,
                &mut scratch,
                out_row,
                true,
            );
        }
    }
    Ok(Image::new_2d(w, h, out)?)
}

/// Separable 2D convolution with zero padding. `kernel_h` filters each row and
/// `kernel_v` filters each column.
///
/// # Errors
/// Returns [`FilterError`] if the image size is wrong or a kernel is empty.
pub fn convolve_separable(
    input: &[f32],
    width: usize,
    height: usize,
    kernel_h: &[f32],
    kernel_v: &[f32],
) -> Result<Vec<f32>, FilterError> {
    let expected = width * height;
    if input.len() != expected {
        return Err(FilterError::SizeMismatch { expected, actual: input.len() });
    }
    if kernel_h.is_empty() || kernel_v.is_empty() {
        return Err(FilterError::EmptyKernel);
    }
    let mut output = vec![0.0f32; expected];
    let mut scratch = vec![0.0f32; expected];
    tpt_simd::convolve::convolve_2d_separable_f32(
        input,
        width,
        height,
        kernel_h,
        kernel_v,
        &mut output,
        &mut scratch,
    );
    Ok(output)
}

/// A normalized Gaussian kernel with `radius * 2 + 1` taps.
///
/// Returns an empty vector when `sigma` is not positive.
#[must_use]
pub fn gaussian_kernel(sigma: f32, radius: usize) -> Vec<f32> {
    if sigma.is_nan() || sigma <= 0.0 {
        return Vec::new();
    }
    let s2 = 2.0 * sigma * sigma;
    let mut k: Vec<f32> = (0..=2 * radius)
        .map(|i| {
            let x = i as f32 - radius as f32;
            (-x * x / s2).exp()
        })
        .collect();
    let sum: f32 = k.iter().sum();
    for v in &mut k {
        *v /= sum;
    }
    k
}

/// Gaussian blur, applied separably.
///
/// # Errors
/// Returns [`FilterError`] for a size mismatch, or an empty kernel when `sigma` is not positive.
pub fn gaussian_blur(
    input: &[f32],
    width: usize,
    height: usize,
    sigma: f32,
) -> Result<Vec<f32>, FilterError> {
    let radius = (3.0 * sigma).ceil().max(1.0) as usize;
    let k = gaussian_kernel(sigma, radius);
    convolve_separable(input, width, height, &k, &k)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Naive reference with the same convention as `tpt-simd`: zero padding and
    /// true convolution, centred at `len / 2`.
    fn reference(input: &[f32], w: usize, h: usize, kh: &[f32], kv: &[f32]) -> Vec<f32> {
        let (ch, cv) = (kh.len() as isize / 2, kv.len() as isize / 2);
        let mut tmp = vec![0.0f32; w * h];
        for y in 0..h {
            for x in 0..w {
                let mut acc = 0.0f32;
                for (i, &kv_i) in kh.iter().enumerate() {
                    let sx = x as isize + ch - i as isize;
                    if (0..w as isize).contains(&sx) {
                        acc += input[y * w + sx as usize] * kv_i;
                    }
                }
                tmp[y * w + x] = acc;
            }
        }
        let mut out = vec![0.0f32; w * h];
        for y in 0..h {
            for x in 0..w {
                let mut acc = 0.0f32;
                for (i, &kv_i) in kv.iter().enumerate() {
                    let sy = y as isize + cv - i as isize;
                    if (0..h as isize).contains(&sy) {
                        acc += tmp[sy as usize * w + x] * kv_i;
                    }
                }
                out[y * w + x] = acc;
            }
        }
        out
    }

    #[test]
    fn matches_reference_on_irregular_sizes() {
        // Sizes that do not divide the SIMD width, to exercise the edge paths.
        for &(w, h) in &[(1usize, 1usize), (3, 3), (13, 7), (37, 19)] {
            let input: Vec<f32> = (0..w * h).map(|i| ((i * 7919) % 97) as f32 / 13.0).collect();
            let kh = [0.25f32, 0.5, 0.25];
            let kv = [1.0f32, -2.0, 1.0, 0.5, 0.0];
            let got = convolve_separable(&input, w, h, &kh, &kv).unwrap();
            let want = reference(&input, w, h, &kh, &kv);
            for (g, r) in got.iter().zip(&want) {
                assert!((g - r).abs() < 1e-4, "{w}x{h}: {g} vs {r}");
            }
        }
    }

    #[test]
    fn identity_kernel_leaves_image_unchanged() {
        let input = [1.0f32, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0];
        assert_eq!(convolve_separable(&input, 3, 3, &[1.0], &[1.0]).unwrap(), input);
    }

    #[test]
    fn gaussian_blur_preserves_constant_interior() {
        let input = vec![5.0f32; 20 * 20];
        let out = gaussian_blur(&input, 20, 20, 1.5).unwrap();
        // Interior pixels are far enough from the zero-padded border to stay at 5.
        assert!((out[10 * 20 + 10] - 5.0).abs() < 1e-4);
    }

    #[test]
    fn boundary_index_maps_outside_samples() {
        assert_eq!(Boundary::Zero.index(-1, 4), None);
        assert_eq!(Boundary::Zero.index(4, 4), None);
        assert_eq!(Boundary::Clamp.index(-3, 4), Some(0));
        assert_eq!(Boundary::Clamp.index(9, 4), Some(3));
        // Mirror: -1 -> 1, 4 -> 2 for n = 4 (period 6).
        assert_eq!(Boundary::Mirror.index(-1, 4), Some(1));
        assert_eq!(Boundary::Mirror.index(4, 4), Some(2));
        assert_eq!(Boundary::Mirror.index(-6, 4), Some(0));
        assert_eq!(Boundary::Mirror.index(-5, 4), Some(1));
        assert_eq!(Boundary::Mirror.index(7, 1), Some(0));
        assert_eq!(Boundary::Clamp.index(0, 0), None);
    }

    /// Naive 2D reference for `convolve_2d`, using the same boundary mapping.
    fn reference_2d(
        input: &[f32],
        w: usize,
        h: usize,
        k: &[f32],
        s: usize,
        b: Boundary,
    ) -> Vec<f32> {
        let c = (s / 2) as isize;
        let mut out = vec![0.0f32; w * h];
        for y in 0..h {
            for x in 0..w {
                let mut acc = 0.0f32;
                for j in 0..s {
                    for i in 0..s {
                        let sy = b.index(y as isize + c - j as isize, h);
                        let sx = b.index(x as isize + c - i as isize, w);
                        if let (Some(sy), Some(sx)) = (sy, sx) {
                            acc += input[sy * w + sx] * k[j * s + i];
                        }
                    }
                }
                out[y * w + x] = acc;
            }
        }
        out
    }

    #[test]
    fn square_kernel_matches_brute_force_on_irregular_sizes() {
        for b in [Boundary::Zero, Boundary::Clamp, Boundary::Mirror] {
            for &(w, h) in &[(1usize, 1usize), (2, 5), (13, 7), (37, 19)] {
                for &s in &[1usize, 3, 5] {
                    let input: Vec<f32> =
                        (0..w * h).map(|i| ((i * 7919) % 97) as f32 / 13.0 - 3.0).collect();
                    let k: Vec<f32> = (0..s * s).map(|i| ((i * 31) % 11) as f32 - 5.0).collect();
                    let img = Image::new_2d(w, h, input.clone()).unwrap();
                    let got = convolve_2d(&img, &k, s, b).unwrap();
                    let want = reference_2d(&input, w, h, &k, s, b);
                    for (g, r) in got.data().iter().zip(&want) {
                        assert!((g - r).abs() < 1e-3, "{b:?} {w}x{h} s={s}: {g} vs {r}");
                    }
                }
            }
        }
    }

    #[test]
    fn outer_product_square_kernel_matches_separable_path() {
        let (w, h) = (13, 9);
        let input: Vec<f32> = (0..w * h).map(|i| ((i * 7919) % 97) as f32 / 13.0).collect();
        let kh = [0.25f32, 0.5, 0.25];
        let kv = [1.0f32, -2.0, 1.0];
        let outer: Vec<f32> = kv.iter().flat_map(|&v| kh.iter().map(move |&u| u * v)).collect();
        let img = Image::new_2d(w, h, input.clone()).unwrap();
        let general = convolve_2d(&img, &outer, 3, Boundary::Zero).unwrap();
        let separable = convolve_separable(&input, w, h, &kh, &kv).unwrap();
        for (g, s) in general.data().iter().zip(&separable) {
            assert!((g - s).abs() < 1e-4, "{g} vs {s}");
        }
    }

    #[test]
    fn square_kernel_rejects_bad_shapes() {
        let img = Image::new_2d(3, 3, vec![0.0f32; 9]).unwrap();
        assert_eq!(
            convolve_2d(&img, &[1.0; 4], 2, Boundary::Zero),
            Err(FilterError::InvalidKernelSize { size: 2 })
        );
        assert_eq!(
            convolve_2d(&img, &[1.0; 8], 3, Boundary::Zero),
            Err(FilterError::KernelLength { expected: 9, actual: 8 })
        );
        let vol = Image::new_3d(3, 3, 2, vec![0.0f32; 18]).unwrap();
        assert_eq!(
            convolve_2d(&vol, &[1.0], 1, Boundary::Zero),
            Err(FilterError::NotTwoDimensional { depth: 2 })
        );
    }

    #[test]
    fn rejects_bad_inputs() {
        assert_eq!(
            convolve_separable(&[1.0, 2.0], 3, 1, &[1.0], &[1.0]),
            Err(FilterError::SizeMismatch { expected: 3, actual: 2 })
        );
        assert_eq!(convolve_separable(&[1.0], 1, 1, &[], &[1.0]), Err(FilterError::EmptyKernel));
        assert!(gaussian_kernel(0.0, 2).is_empty());
    }
}
