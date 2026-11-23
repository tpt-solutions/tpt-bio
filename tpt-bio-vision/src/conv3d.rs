//! Separable 3D convolution and Gaussian blur on [`Image<f32>`] (ADR-002).
//!
//! The volume is filtered by three 1D passes, along `x`, then `y`, then `z`. Each pass
//! runs every line through [`crate::filter`]'s padded line convolution, which uses the
//! `tpt-simd` 1D kernel. Out-of-range samples follow a [`Boundary`] mode.

use crate::filter::{Boundary, FilterError, LineScratch, convolve_line, gaussian_kernel};
use crate::image::Image;

/// Separable 3D convolution: `kernel_x` along `x`, `kernel_y` along `y`, `kernel_z` along `z`.
///
/// Each kernel is centred at `len / 2`, so the output at `(x, y, z)` is
/// `sum over (i, j, l) of kx[i] * ky[j] * kz[l] * in[x + cx - i][y + cy - j][z + cz - l]`,
/// with out-of-range samples taken from `boundary`. The output has the input's shape.
///
/// # Errors
/// Returns [`FilterError::EmptyKernel`] if any kernel is empty.
pub fn convolve_separable_3d(
    input: &Image<f32>,
    kernel_x: &[f32],
    kernel_y: &[f32],
    kernel_z: &[f32],
    boundary: Boundary,
) -> Result<Image<f32>, FilterError> {
    if kernel_x.is_empty() || kernel_y.is_empty() || kernel_z.is_empty() {
        return Err(FilterError::EmptyKernel);
    }
    let shape = input.shape();
    let mut data = input.data().to_vec();
    for (axis, kernel) in [(0, kernel_x), (1, kernel_y), (2, kernel_z)] {
        data = convolve_axis(&data, shape, axis, kernel, boundary);
    }
    Ok(Image::new_3d(shape.0, shape.1, shape.2, data)?)
}

/// Gaussian blur of a volume, applied through [`convolve_separable_3d`].
///
/// The kernel is [`gaussian_kernel`] with radius `ceil(truncate * sigma)` (at least 1), so
/// a `truncate` of 4.0 keeps the taps out to four standard deviations. The kernel is
/// normalised to sum to one.
///
/// # Errors
/// Returns [`FilterError::EmptyKernel`] when `sigma` is not finite and positive.
pub fn gaussian_blur_3d(
    input: &Image<f32>,
    sigma: f32,
    truncate: f32,
    boundary: Boundary,
) -> Result<Image<f32>, FilterError> {
    if !sigma.is_finite() || sigma <= 0.0 {
        return Err(FilterError::EmptyKernel);
    }
    let radius = (truncate * sigma).ceil().max(1.0) as usize;
    let k = gaussian_kernel(sigma, radius);
    convolve_separable_3d(input, &k, &k, &k, boundary)
}

/// Convolves every line of `src` along `axis` (0 = x, 1 = y, 2 = z).
fn convolve_axis(
    src: &[f32],
    shape: (usize, usize, usize),
    axis: usize,
    kernel: &[f32],
    boundary: Boundary,
) -> Vec<f32> {
    let (w, h, d) = shape;
    let (stride, len) = match axis {
        0 => (1, w),
        1 => (w, h),
        _ => (w * h, d),
    };
    let total = src.len();
    let mut out = vec![0.0f32; total];
    let mut line = vec![0.0f32; len];
    let mut res = vec![0.0f32; len];
    let mut scratch = LineScratch::default();
    for base in 0..total {
        // Only the line starts (coordinate 0 along `axis`) are visited.
        if (base / stride) % len != 0 {
            continue;
        }
        for (i, v) in line.iter_mut().enumerate() {
            *v = src[base + i * stride];
        }
        convolve_line(&line, kernel, boundary, &mut scratch, &mut res, false);
        for (i, &v) in res.iter().enumerate() {
            out[base + i * stride] = v;
        }
    }
    out
}

#[allow(clippy::float_cmp)]
#[cfg(test)]
mod tests {
    use super::*;

    /// Naive triple-loop reference: direct 3D convolution with the product kernel.
    fn brute_3d(input: &Image<f32>, kx: &[f32], ky: &[f32], kz: &[f32], b: Boundary) -> Vec<f32> {
        let (w, h, d) = input.shape();
        let (cx, cy, cz) =
            ((kx.len() / 2) as isize, (ky.len() / 2) as isize, (kz.len() / 2) as isize);
        let mut out = vec![0.0f32; w * h * d];
        for z in 0..d {
            for y in 0..h {
                for x in 0..w {
                    let mut acc = 0.0f32;
                    for (l, &vz) in kz.iter().enumerate() {
                        for (j, &vy) in ky.iter().enumerate() {
                            for (i, &vx) in kx.iter().enumerate() {
                                let sx = b.index(x as isize + cx - i as isize, w);
                                let sy = b.index(y as isize + cy - j as isize, h);
                                let sz = b.index(z as isize + cz - l as isize, d);
                                if let (Some(sx), Some(sy), Some(sz)) = (sx, sy, sz) {
                                    acc += input.data()[input.index(sx, sy, sz)] * vx * vy * vz;
                                }
                            }
                        }
                    }
                    out[input.index(x, y, z)] = acc;
                }
            }
        }
        out
    }

    #[test]
    fn separable_3d_matches_brute_force_on_irregular_sizes() {
        let kx = [0.25f32, 0.5, 0.25];
        let ky = [1.0f32, -2.0, 1.0, 0.5];
        let kz = [0.1f32, 0.2, 0.4, 0.2, 0.1];
        for b in [Boundary::Zero, Boundary::Clamp, Boundary::Mirror] {
            for &(w, h, d) in &[(1usize, 1usize, 1usize), (3, 2, 1), (7, 5, 3), (13, 9, 4)] {
                let data: Vec<f32> =
                    (0..w * h * d).map(|i| ((i * 7919) % 97) as f32 / 13.0 - 3.0).collect();
                let img = Image::new_3d(w, h, d, data).unwrap();
                let got = convolve_separable_3d(&img, &kx, &ky, &kz, b).unwrap();
                let want = brute_3d(&img, &kx, &ky, &kz, b);
                for (g, r) in got.data().iter().zip(&want) {
                    assert!((g - r).abs() < 1e-3, "{b:?} {w}x{h}x{d}: {g} vs {r}");
                }
            }
        }
    }

    #[test]
    fn identity_kernels_leave_volume_unchanged() {
        let data: Vec<f32> = (0..24).map(|i| i as f32).collect();
        let img = Image::new_3d(4, 3, 2, data.clone()).unwrap();
        let out = convolve_separable_3d(&img, &[1.0], &[1.0], &[1.0], Boundary::Zero).unwrap();
        assert_eq!(out.data(), &data[..]);
    }

    #[test]
    fn rejects_empty_kernel() {
        let img = Image::new_3d(2, 2, 2, vec![0.0f32; 8]).unwrap();
        assert_eq!(
            convolve_separable_3d(&img, &[1.0], &[], &[1.0], Boundary::Zero),
            Err(FilterError::EmptyKernel)
        );
        assert_eq!(gaussian_blur_3d(&img, 0.0, 4.0, Boundary::Zero), Err(FilterError::EmptyKernel));
    }

    #[test]
    fn constant_volume_stays_constant_under_clamp_and_mirror() {
        let img = Image::new_3d(7, 5, 6, vec![3.5f32; 7 * 5 * 6]).unwrap();
        for b in [Boundary::Clamp, Boundary::Mirror] {
            let out = gaussian_blur_3d(&img, 1.5, 4.0, b).unwrap();
            for &v in out.data() {
                assert!((v - 3.5).abs() < 1e-4, "{b:?}: {v}");
            }
        }
    }

    #[test]
    fn constant_volume_interior_stays_constant_under_zero_boundary() {
        let img = Image::new_3d(21, 21, 21, vec![3.5f32; 21 * 21 * 21]).unwrap();
        let out = gaussian_blur_3d(&img, 1.0, 4.0, Boundary::Zero).unwrap();
        assert!((out.get(10, 10, 10).unwrap() - 3.5).abs() < 1e-4);
    }

    #[test]
    fn unit_impulse_spreads_to_normalised_gaussian_weights() {
        let n = 21;
        let mut data = vec![0.0f32; n * n * n];
        let c = 10;
        let img_tmp = Image::new_3d(n, n, n, vec![0.0f32; n * n * n]).unwrap();
        data[img_tmp.index(c, c, c)] = 1.0;
        let img = Image::new_3d(n, n, n, data).unwrap();
        let sigma = 1.0f32;
        let out = gaussian_blur_3d(&img, sigma, 4.0, Boundary::Zero).unwrap();
        let k = gaussian_kernel(sigma, 4);
        assert_eq!(k.len(), 9);
        // Impulse at the centre: each voxel gets the product of 1D weights.
        for dz in -4isize..=4 {
            for dy in -4isize..=4 {
                for dx in -4isize..=4 {
                    let (x, y, z) = (
                        (c as isize + dx) as usize,
                        (c as isize + dy) as usize,
                        (c as isize + dz) as usize,
                    );
                    let want = k[(dx + 4) as usize] * k[(dy + 4) as usize] * k[(dz + 4) as usize];
                    let got = *out.get(x, y, z).unwrap();
                    assert!((got - want).abs() < 1e-6, "({dx},{dy},{dz}): {got} vs {want}");
                }
            }
        }
        // Total mass is one: the truncated kernel is normalised and the impulse is interior.
        let total: f32 = out.data().iter().sum();
        assert!((total - 1.0).abs() < 1e-5, "total {total}");
        // Centre weight is k[r]^3 and is symmetric in every axis.
        assert!((out.get(c, c, c).unwrap() - k[4] * k[4] * k[4]).abs() < 1e-6);
        assert_eq!(out.get(c - 1, c, c), out.get(c + 1, c, c));
    }
}
