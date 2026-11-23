//! Image filters built on `tpt-simd` (ADR-002).
//!
//! Images are row-major `f32` planes. Convolution keeps the input size and uses zero
//! padding, with each kernel centred at `len / 2`.

use std::fmt;

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
}

impl fmt::Display for FilterError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::SizeMismatch { expected, actual } => {
                write!(f, "image has {actual} samples, expected {expected}")
            }
            Self::EmptyKernel => f.write_str("convolution kernel is empty"),
        }
    }
}

impl std::error::Error for FilterError {}

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
    fn rejects_bad_inputs() {
        assert_eq!(
            convolve_separable(&[1.0, 2.0], 3, 1, &[1.0], &[1.0]),
            Err(FilterError::SizeMismatch { expected: 3, actual: 2 })
        );
        assert_eq!(convolve_separable(&[1.0], 1, 1, &[], &[1.0]), Err(FilterError::EmptyKernel));
        assert!(gaussian_kernel(0.0, 2).is_empty());
    }
}
