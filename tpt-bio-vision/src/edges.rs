//! Edge detection in 2D: Sobel and Scharr gradients, and Canny.
//!
//! Gradients use 3x3 kernels with edge replication at the border. Canny here does
//! not smooth the input; call [`crate::spatial::gaussian_smooth`] first if needed.

use std::collections::VecDeque;

use crate::image::{Image, ImageError};

/// A 3x3 gradient kernel pair.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum GradientKernel {
    /// Sobel: smoothing weights `1, 2, 1`.
    Sobel,
    /// Scharr: smoothing weights `3, 10, 3`, a more rotation-invariant gradient.
    Scharr,
}

impl GradientKernel {
    /// Smoothing weights of the 3-tap perpendicular direction, `[w_-1, w_0, w_1]`.
    fn smooth(self) -> [f64; 3] {
        match self {
            Self::Sobel => [1.0, 2.0, 1.0],
            Self::Scharr => [3.0, 10.0, 3.0],
        }
    }
}

/// Horizontal and vertical gradients `(gx, gy)` with edge replication.
///
/// `gx` is the difference across `x` (`right - left`), so a ramp `p = x` gives
/// `gx = 8` for Sobel and `32` for Scharr away from the border.
///
/// # Errors
/// Returns [`ImageError::NotTwoDimensional`] for volumes.
pub fn gradients(
    image: &Image<f64>,
    kernel: GradientKernel,
) -> Result<(Image<f64>, Image<f64>), ImageError> {
    if image.depth() != 1 {
        return Err(ImageError::NotTwoDimensional { depth: image.depth() });
    }
    let (w, h) = (image.width(), image.height());
    let s = kernel.smooth();
    let at = |x: isize, y: isize| {
        let cx = x.clamp(0, w as isize - 1) as usize;
        let cy = y.clamp(0, h as isize - 1) as usize;
        image.data()[cy * w + cx]
    };
    let mut gx = Vec::with_capacity(w * h);
    let mut gy = Vec::with_capacity(w * h);
    for y in 0..h as isize {
        for x in 0..w as isize {
            // Horizontal derivative: smooth vertically, difference horizontally.
            let dx: f64 = (-1..=1)
                .map(|j| {
                    let wgt = s[(j + 1) as usize];
                    wgt * (at(x + 1, y + j) - at(x - 1, y + j))
                })
                .sum();
            // Vertical derivative: smooth horizontally, difference vertically.
            let dy: f64 = (-1..=1)
                .map(|i| {
                    let wgt = s[(i + 1) as usize];
                    wgt * (at(x + i, y + 1) - at(x + i, y - 1))
                })
                .sum();
            gx.push(dx);
            gy.push(dy);
        }
    }
    Ok((Image::new_2d(w, h, gx)?, Image::new_2d(w, h, gy)?))
}

/// Gradient magnitude `sqrt(gx² + gy²)` with the given kernel.
///
/// # Errors
/// Returns [`ImageError::NotTwoDimensional`] for volumes.
pub fn gradient_magnitude(
    image: &Image<f64>,
    kernel: GradientKernel,
) -> Result<Image<f64>, ImageError> {
    let (gx, gy) = gradients(image, kernel)?;
    let data = gx.data().iter().zip(gy.data()).map(|(a, b)| a.hypot(*b)).collect();
    Image::new_2d(image.width(), image.height(), data)
}

/// Sobel gradient magnitude.
///
/// # Errors
/// Returns [`ImageError::NotTwoDimensional`] for volumes.
pub fn sobel(image: &Image<f64>) -> Result<Image<f64>, ImageError> {
    gradient_magnitude(image, GradientKernel::Sobel)
}

/// Scharr gradient magnitude.
///
/// # Errors
/// Returns [`ImageError::NotTwoDimensional`] for volumes.
pub fn scharr(image: &Image<f64>) -> Result<Image<f64>, ImageError> {
    gradient_magnitude(image, GradientKernel::Scharr)
}

/// Canny edges: Sobel gradients, non-maximum suppression along the gradient direction,
/// then hysteresis. A pixel is an edge if it is a local maximum with magnitude at
/// least `low` and is connected (8-neighbourhood, through such pixels) to a pixel of
/// magnitude at least `high`.
///
/// Suppression keeps a pixel when its magnitude is at least that of the neighbour on
/// one side and strictly greater than the neighbour on the other, so a two-pixel-wide
/// step produces a one-pixel edge.
///
/// # Errors
/// Returns [`ImageError::NotTwoDimensional`] for volumes, or [`ImageError::NonFinite`]
/// for NaN input.
pub fn canny(image: &Image<f64>, low: f64, high: f64) -> Result<Image<bool>, ImageError> {
    if image.data().iter().any(|v| v.is_nan()) {
        return Err(ImageError::NonFinite);
    }
    let (gx, gy) = gradients(image, GradientKernel::Sobel)?;
    let (w, h) = (image.width(), image.height());
    let mag: Vec<f64> = gx.data().iter().zip(gy.data()).map(|(a, b)| a.hypot(*b)).collect();
    let magnitude_at = |x: isize, y: isize| {
        if x < 0 || y < 0 || x >= w as isize || y >= h as isize {
            0.0
        } else {
            mag[y as usize * w + x as usize]
        }
    };

    let mut nms = vec![false; w * h];
    for y in 0..h {
        for x in 0..w {
            let idx = y * w + x;
            let (a, b) = (gx.data()[idx], gy.data()[idx]);
            let m = mag[idx];
            if m <= 0.0 {
                continue;
            }
            // Quantise the gradient direction to one of four neighbour pairs.
            let (ax, ay) = (a.abs(), b.abs());
            let (dx, dy): (isize, isize) = if ay <= 0.414_213_56 * ax {
                (1, 0)
            } else if ax <= 0.414_213_56 * ay {
                (0, 1)
            } else if (a > 0.0) == (b > 0.0) {
                (1, 1)
            } else {
                (1, -1)
            };
            let (xi, yi) = (x as isize, y as isize);
            let before = magnitude_at(xi - dx, yi - dy);
            let after = magnitude_at(xi + dx, yi + dy);
            nms[idx] = m >= before && m > after;
        }
    }

    let strong: Vec<bool> = (0..w * h).map(|i| nms[i] && mag[i] >= high).collect();
    let weak: Vec<bool> = (0..w * h).map(|i| nms[i] && mag[i] >= low).collect();
    hysteresis(&Image::new_2d(w, h, strong)?, &Image::new_2d(w, h, weak)?)
}

/// Pixels of `weak` that are 8-connected, through `weak`, to a `strong` pixel.
/// `strong` is expected to be a subset of `weak`; extra strong pixels are kept as is.
///
/// # Errors
/// Returns [`ImageError::ShapeMismatch`] if the two masks differ in shape.
pub fn hysteresis(strong: &Image<bool>, weak: &Image<bool>) -> Result<Image<bool>, ImageError> {
    strong.check_same_shape(weak)?;
    let mut keep = vec![false; strong.len()];
    let mut queue = VecDeque::new();
    for (i, &s) in strong.data().iter().enumerate() {
        if s {
            keep[i] = true;
            queue.push_back(i);
        }
    }
    while let Some(p) = queue.pop_front() {
        for &d in &crate::image::OFFSETS_8 {
            if let Some(n) = weak.offset(p, d)
                && weak.data()[n]
                && !keep[n]
            {
                keep[n] = true;
                queue.push_back(n);
            }
        }
    }
    Image::new_3d(strong.width(), strong.height(), strong.depth(), keep)
}

#[allow(clippy::float_cmp)]
#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn sobel_and_scharr_on_a_ramp() {
        let img = Image::new_2d(5, 3, (0..15).map(|i| f64::from(i % 5)).collect()).unwrap();
        let (gx, gy) = gradients(&img, GradientKernel::Sobel).unwrap();
        // Interior: smoothed right column 4*(x+1) minus smoothed left column 4*(x-1) = 8.
        assert_eq!(*gx.get(2, 1, 0).unwrap(), 8.0);
        assert_eq!(*gy.get(2, 1, 0).unwrap(), 0.0);
        // Left border replicates x=0: right column 4*1 minus left column 4*0 = 4.
        assert_eq!(*gx.get(0, 1, 0).unwrap(), 4.0);
        let (sx, _) = gradients(&img, GradientKernel::Scharr).unwrap();
        assert_eq!(*sx.get(2, 1, 0).unwrap(), 32.0);
    }

    #[test]
    fn vertical_ramp_gives_vertical_gradient() {
        let img = Image::new_2d(3, 5, (0..15).map(|i| f64::from(i / 3)).collect()).unwrap();
        let (gx, gy) = gradients(&img, GradientKernel::Sobel).unwrap();
        assert_eq!(*gy.get(1, 2, 0).unwrap(), 8.0);
        assert_eq!(*gx.get(1, 2, 0).unwrap(), 0.0);
        assert_eq!(*sobel(&img).unwrap().get(1, 2, 0).unwrap(), 8.0);
    }

    #[test]
    fn canny_finds_a_vertical_step_as_one_column() {
        // Step from 0 to 10 between x=2 and x=3. Sobel magnitude is 40 at x=2 and x=3.
        let data: Vec<f64> = (0..6 * 5).map(|i| if i % 6 < 3 { 0.0 } else { 10.0 }).collect();
        let img = Image::new_2d(6, 5, data).unwrap();
        let edges = canny(&img, 1.0, 1.0).unwrap();
        for y in 0..5 {
            for x in 0..6 {
                let expect = x == 3;
                assert_eq!(*edges.get(x, y, 0).unwrap(), expect, "pixel ({x},{y})");
            }
        }
    }

    #[test]
    fn canny_flat_image_has_no_edges() {
        let img = Image::new_2d(4, 4, vec![3.0; 16]).unwrap();
        let edges = canny(&img, 0.0, 0.0).unwrap();
        assert!(edges.data().iter().all(|&b| !b));
    }

    #[test]
    fn hysteresis_keeps_weak_pixels_only_when_connected_to_strong() {
        // Row of 5 weak pixels; the strong pixel is at x=0. Weak run 0..=2 is kept,
        // and the isolated weak pixel at x=4 is dropped.
        let weak = Image::new_2d(5, 1, vec![true, true, true, false, true]).unwrap();
        let strong = Image::new_2d(5, 1, vec![true, false, false, false, false]).unwrap();
        let kept = hysteresis(&strong, &weak).unwrap();
        assert_eq!(kept.data(), &[true, true, true, false, false]);
    }

    #[test]
    fn rejects_volumes() {
        let vol = Image::new_3d(2, 2, 2, vec![0.0; 8]).unwrap();
        assert!(sobel(&vol).is_err());
        assert!(canny(&vol, 0.0, 1.0).is_err());
    }
}
