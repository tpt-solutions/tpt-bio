//! Seed detection: local maxima with a minimum separation, and Laplacian-of-Gaussian
//! blob seeds. Seeds are returned as linear pixel indices, strongest first.
//! [`seeds_to_markers`] turns them into a marker image for [`crate::watershed`].

use crate::image::{Image, ImageError};
use crate::spatial::gaussian_smooth;

/// A seed: linear pixel index and the response value at that pixel.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Seed {
    /// Linear index into the image (see [`Image::index`]).
    pub index: usize,
    /// Response value used for ranking.
    pub value: f64,
}

/// Local maxima of `image` with a suppression radius of `min_distance` (Chebyshev).
///
/// A pixel is a seed if its value is at least `min_value` and no pixel within
/// `min_distance` (in all three axes; the z range is 0 for 2D images) is larger.
/// Plateaus yield one seed, the first in raster order. Results are sorted by value,
/// strongest first, with ties in raster order.
///
/// # Errors
/// Returns [`ImageError::NonFinite`] if the image contains NaN.
pub fn local_maxima(
    image: &Image<f64>,
    min_distance: usize,
    min_value: f64,
) -> Result<Vec<Seed>, ImageError> {
    if image.data().iter().any(|v| v.is_nan()) {
        return Err(ImageError::NonFinite);
    }
    let (w, h, d) = image.shape();
    let r = min_distance as isize;
    let rz = if d > 1 { r } else { 0 };
    let mut seeds = Vec::new();
    for idx in 0..image.len() {
        let v = image.data()[idx];
        if v < min_value {
            continue;
        }
        let (x, y, z) = image.coords(idx);
        let mut is_max = true;
        'scan: for dz in -rz..=rz {
            for dy in -r..=r {
                for dx in -r..=r {
                    if dx == 0 && dy == 0 && dz == 0 {
                        continue;
                    }
                    let (nx, ny, nz) = (x as isize + dx, y as isize + dy, z as isize + dz);
                    if nx < 0 || ny < 0 || nz < 0 {
                        continue;
                    }
                    let (nx, ny, nz) = (nx as usize, ny as usize, nz as usize);
                    if nx >= w || ny >= h || nz >= d {
                        continue;
                    }
                    let n = image.data()[image.index(nx, ny, nz)];
                    // Strictly larger neighbours disqualify; equal neighbours before this
                    // pixel in raster order win the plateau.
                    if n > v || (n.total_cmp(&v).is_eq() && image.index(nx, ny, nz) < idx) {
                        is_max = false;
                        break 'scan;
                    }
                }
            }
        }
        if is_max {
            seeds.push(Seed { index: idx, value: v });
        }
    }
    seeds.sort_by(|a, b| b.value.total_cmp(&a.value).then(a.index.cmp(&b.index)));
    Ok(seeds)
}

/// Blob response `-sigma² * Laplacian(G_sigma * image)` for 2D images.
///
/// Bright blobs of radius about `sigma * sqrt(2)` give a positive peak at their centre.
/// The Laplacian is the 4-neighbour stencil with edge replication.
///
/// # Errors
/// Returns [`ImageError::NotTwoDimensional`] for volumes.
pub fn log_response(image: &Image<f64>, sigma: f64) -> Result<Image<f64>, ImageError> {
    let smooth = gaussian_smooth(image, sigma)?;
    let (w, h) = (smooth.width(), smooth.height());
    let s = smooth.data();
    let at = |x: isize, y: isize| {
        let cx = x.clamp(0, w as isize - 1) as usize;
        let cy = y.clamp(0, h as isize - 1) as usize;
        s[cy * w + cx]
    };
    let mut out = Vec::with_capacity(w * h);
    for y in 0..h as isize {
        for x in 0..w as isize {
            let lap = at(x + 1, y) + at(x - 1, y) + at(x, y + 1) + at(x, y - 1) - 4.0 * at(x, y);
            out.push(-sigma * sigma * lap);
        }
    }
    Image::new_2d(w, h, out)
}

/// Seeds from local maxima of the `LoG` blob response at scale `sigma`.
///
/// # Errors
/// Returns [`ImageError::NotTwoDimensional`] for volumes, or [`ImageError::NonFinite`]
/// for NaN input.
pub fn log_seeds(
    image: &Image<f64>,
    sigma: f64,
    min_distance: usize,
    min_response: f64,
) -> Result<Vec<Seed>, ImageError> {
    let response = log_response(image, sigma)?;
    local_maxima(&response, min_distance, min_response)
}

/// Marker image for watershed: seed `i` (in the given order) gets label `i + 1`.
/// Labels are `u32`, so more than `u32::MAX` seeds are not representable.
///
/// # Errors
/// Returns [`ImageError::SizeMismatch`] if a seed index lies outside the shape, which
/// is checked by the index bound against `width * height * depth`.
pub fn seeds_to_markers(
    width: usize,
    height: usize,
    depth: usize,
    seeds: &[Seed],
) -> Result<Image<u32>, ImageError> {
    let n = width * height * depth;
    let mut data = vec![0u32; n];
    for (i, s) in seeds.iter().enumerate() {
        if s.index >= n {
            return Err(ImageError::SizeMismatch { expected: n, actual: s.index + 1 });
        }
        data[s.index] = i as u32 + 1;
    }
    Image::new_3d(width, height, depth, data)
}

#[allow(clippy::float_cmp)]
#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn two_peaks_with_separation() {
        // Peaks at x=1 (value 5) and x=5 (value 3) in a 1x7 row.
        let img = Image::new_2d(7, 1, vec![0.0, 5.0, 1.0, 0.0, 2.0, 3.0, 0.0]).unwrap();
        let seeds = local_maxima(&img, 1, 0.0).unwrap();
        let idx: Vec<usize> = seeds.iter().map(|s| s.index).collect();
        assert_eq!(idx, vec![1, 5]);
        // Requiring separation 4 suppresses the weaker peak.
        let far = local_maxima(&img, 4, 0.0).unwrap();
        assert_eq!(far.len(), 1);
        assert_eq!(far[0].index, 1);
    }

    #[test]
    fn min_value_filters_weak_maxima() {
        let img = Image::new_2d(5, 1, vec![0.0, 1.0, 0.0, 4.0, 0.0]).unwrap();
        let seeds = local_maxima(&img, 1, 2.0).unwrap();
        assert_eq!(seeds.len(), 1);
        assert_eq!(seeds[0].index, 3);
    }

    #[test]
    fn plateau_yields_one_seed() {
        let img = Image::new_2d(4, 1, vec![0.0, 3.0, 3.0, 0.0]).unwrap();
        let seeds = local_maxima(&img, 1, 0.0).unwrap();
        assert_eq!(seeds.len(), 1);
        assert_eq!(seeds[0].index, 1);
    }

    #[test]
    fn log_finds_gaussian_blob_centres() {
        let (w, h) = (21usize, 21usize);
        let mut data = vec![0.0; w * h];
        for &(cx, cy) in &[(5.0f64, 5.0f64), (15.0, 14.0)] {
            for y in 0..h {
                for x in 0..w {
                    let d2 = (x as f64 - cx).powi(2) + (y as f64 - cy).powi(2);
                    data[y * w + x] += (-d2 / (2.0 * 2.0 * 2.0)).exp();
                }
            }
        }
        let img = Image::new_2d(w, h, data).unwrap();
        let seeds = log_seeds(&img, 2.0, 5, 0.0).unwrap();
        let coords: Vec<(usize, usize)> =
            seeds.iter().map(|s| (s.index % w, s.index / w)).collect();
        assert_eq!(coords.len(), 2, "found {coords:?}");
        assert!(coords.contains(&(5, 5)) && coords.contains(&(15, 14)), "{coords:?}");
    }

    #[test]
    fn markers_label_seeds_in_order() {
        let seeds = [Seed { index: 4, value: 2.0 }, Seed { index: 0, value: 1.0 }];
        let m = seeds_to_markers(3, 2, 1, &seeds).unwrap();
        assert_eq!(m.data(), &[2, 0, 0, 0, 1, 0]);
    }

    #[test]
    fn flat_image_has_one_seed() {
        // A flat image has no strict maximum; the first pixel in raster order wins.
        let img = Image::new_2d(3, 3, vec![1.0; 9]).unwrap();
        assert_eq!(local_maxima(&img, 1, 0.0).unwrap().len(), 1);
    }
}
