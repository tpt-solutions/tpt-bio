//! Otsu thresholding, global and tile-based adaptive.
//!
//! A threshold `t` splits pixels into foreground (`value >= t`) and background
//! (`value < t`). Otsu's method picks the `t` that maximises the between-class
//! variance of the intensity histogram.

use crate::image::{Image, ImageError};

/// Number of histogram bins used when none is specified.
pub const DEFAULT_BINS: usize = 256;

/// Otsu's split of a histogram: the index of the last background bin.
///
/// Returns `None` when the histogram is empty or has only one non-empty bin
/// (there is nothing to separate).
#[must_use]
pub fn otsu_from_histogram(histogram: &[u64]) -> Option<usize> {
    let total: u64 = histogram.iter().sum();
    let non_empty = histogram.iter().filter(|&&c| c > 0).count();
    if total == 0 || non_empty < 2 {
        return None;
    }
    let sum_all: f64 = histogram.iter().enumerate().map(|(i, &c)| i as f64 * c as f64).sum();
    let total_f = total as f64;

    let mut w_b = 0.0f64;
    let mut sum_b = 0.0f64;
    let mut best = f64::NEG_INFINITY;
    let mut best_bin = 0usize;
    // Only split positions k < last non-empty bin are candidates.
    for (k, &count) in histogram.iter().enumerate().take(histogram.len() - 1) {
        w_b += count as f64;
        if w_b == 0.0 {
            continue;
        }
        let w_f = total_f - w_b;
        if w_f == 0.0 {
            break;
        }
        sum_b += k as f64 * count as f64;
        let mean_b = sum_b / w_b;
        let mean_f = (sum_all - sum_b) / w_f;
        let between = w_b * w_f * (mean_b - mean_f) * (mean_b - mean_f);
        if between > best {
            best = between;
            best_bin = k;
        }
    }
    Some(best_bin)
}

/// Global Otsu threshold of the samples, using `bins` histogram bins over `[min, max]`.
///
/// The returned value `t` satisfies: samples below `t` are background, samples at or
/// above `t` are foreground. Returns `None` when the image is empty or constant.
///
/// # Errors
/// Returns [`ImageError::NonFinite`] if a sample is NaN.
pub fn otsu_threshold<T: Copy + Into<f64>>(
    image: &Image<T>,
    bins: usize,
) -> Result<Option<f64>, ImageError> {
    let values: Vec<f64> = image.data().iter().map(|&v| v.into()).collect();
    otsu_of_values(&values, bins)
}

fn otsu_of_values(values: &[f64], bins: usize) -> Result<Option<f64>, ImageError> {
    let bins = bins.max(2);
    let mut lo = f64::INFINITY;
    let mut hi = f64::NEG_INFINITY;
    for &v in values {
        if v.is_nan() {
            return Err(ImageError::NonFinite);
        }
        lo = lo.min(v);
        hi = hi.max(v);
    }
    if values.is_empty() || hi <= lo {
        return Ok(None);
    }
    let scale = bins as f64 / (hi - lo);
    let mut hist = vec![0u64; bins];
    for &v in values {
        let b = (((v - lo) * scale) as usize).min(bins - 1);
        hist[b] += 1;
    }
    Ok(otsu_from_histogram(&hist).map(|k| lo + (k + 1) as f64 / scale))
}

/// Binary mask of the samples at or above `threshold`.
#[must_use]
pub fn apply_threshold<T: Copy + Into<f64>>(image: &Image<T>, threshold: f64) -> Image<bool> {
    image.map(|&v| v.into() >= threshold)
}

/// Per-pixel Otsu threshold map from non-overlapping cubic-or-square tiles.
///
/// Each tile of `tile` by `tile` pixels (clipped at the image edge) gets its own Otsu
/// threshold, and every pixel receives the threshold of its tile. Constant tiles
/// get `f64::INFINITY`, so they produce no foreground.
///
/// # Errors
/// Returns [`ImageError::NonFinite`] if a sample is NaN, or
/// [`ImageError::NotTwoDimensional`] if the image has depth above 1.
pub fn local_otsu_tiles<T: Copy + Into<f64>>(
    image: &Image<T>,
    tile: usize,
    bins: usize,
) -> Result<Image<f64>, ImageError> {
    if image.depth() != 1 {
        return Err(ImageError::NotTwoDimensional { depth: image.depth() });
    }
    let (w, h) = (image.width(), image.height());
    let tile = tile.max(1);
    let mut out = vec![f64::INFINITY; w * h];
    for ty in (0..h).step_by(tile) {
        for tx in (0..w).step_by(tile) {
            let mut values = Vec::with_capacity(tile * tile);
            for y in ty..(ty + tile).min(h) {
                for x in tx..(tx + tile).min(w) {
                    values.push(image.data()[image.index(x, y, 0)].into());
                }
            }
            let t = otsu_of_values(&values, bins)?.unwrap_or(f64::INFINITY);
            for y in ty..(ty + tile).min(h) {
                for x in tx..(tx + tile).min(w) {
                    out[y * w + x] = t;
                }
            }
        }
    }
    Image::new_2d(w, h, out)
}

/// Per-pixel Otsu threshold map from a sliding square window.
///
/// Each pixel gets the Otsu threshold of the `(2 * radius + 1)` square window centred on
/// it, clipped at the image edge, computed over that window's histogram. A window that
/// holds a single grey level has no Otsu split, so it takes the global Otsu threshold of
/// the whole image instead; if the whole image is constant, that is `f64::INFINITY`.
///
/// Cost: each pixel builds its own histogram over its window, so the work is
/// `O(N * (2r+1)^2 + N * bins)` for `N` pixels. The window is not updated
/// incrementally, which keeps the code simple at the price of a factor of `(2r+1)^2`.
///
/// # Errors
/// Returns [`ImageError::NonFinite`] if a sample is NaN, or
/// [`ImageError::NotTwoDimensional`] if the image has depth above 1.
pub fn local_otsu_window<T: Copy + Into<f64>>(
    image: &Image<T>,
    radius: usize,
    bins: usize,
) -> Result<Image<f64>, ImageError> {
    if image.depth() != 1 {
        return Err(ImageError::NotTwoDimensional { depth: image.depth() });
    }
    let (w, h) = (image.width(), image.height());
    let values: Vec<f64> = image.data().iter().map(|&v| v.into()).collect();
    let global = otsu_of_values(&values, bins)?.unwrap_or(f64::INFINITY);
    let mut out = vec![global; w * h];
    let mut window = Vec::with_capacity((2 * radius + 1) * (2 * radius + 1));
    for y in 0..h {
        for x in 0..w {
            window.clear();
            for wy in y.saturating_sub(radius)..(y + radius + 1).min(h) {
                for wx in x.saturating_sub(radius)..(x + radius + 1).min(w) {
                    window.push(values[wy * w + wx]);
                }
            }
            if let Some(t) = otsu_of_values(&window, bins)? {
                out[y * w + x] = t;
            }
        }
    }
    Image::new_2d(w, h, out)
}

/// Foreground mask from a per-pixel threshold map (`value >= threshold[pixel]`).
///
/// # Errors
/// Returns [`ImageError::ShapeMismatch`] if the image and threshold map differ in shape.
pub fn apply_threshold_map<T: Copy + Into<f64>>(
    image: &Image<T>,
    thresholds: &Image<f64>,
) -> Result<Image<bool>, ImageError> {
    image.check_same_shape(thresholds)?;
    let data = image.data().iter().zip(thresholds.data()).map(|(&v, &t)| v.into() >= t).collect();
    Image::new_3d(image.width(), image.height(), image.depth(), data)
}

#[allow(clippy::float_cmp)]
#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn histogram_split_between_two_peaks() {
        // Background at bin 1 (count 10), foreground at bin 6 (count 10).
        let mut h = vec![0u64; 8];
        h[1] = 10;
        h[6] = 10;
        // Any split in 1..6 is optimal; the first one is returned.
        assert_eq!(otsu_from_histogram(&h), Some(1));
    }

    #[test]
    fn histogram_split_is_exact_for_hand_case() {
        // Values 0 (x3), 1 (x1), 2 (x1), 9 (x3). Background {0,1,2} vs {9}.
        // Split after bin 2 gives w_b=5, mean_b=0.6, w_f=3, mean_f=9.
        // Between-class variance = 5*3*(8.4)^2 = 1058.4, the maximum.
        let mut h = vec![0u64; 10];
        h[0] = 3;
        h[1] = 1;
        h[2] = 1;
        h[9] = 3;
        assert_eq!(otsu_from_histogram(&h), Some(2));
    }

    #[test]
    fn single_populated_bin_has_no_split() {
        assert_eq!(otsu_from_histogram(&[0, 5, 0]), None);
        assert_eq!(otsu_from_histogram(&[]), None);
    }

    #[test]
    fn global_threshold_separates_two_levels() {
        let data: Vec<u8> = (0..64).map(|i| if i < 32 { 10 } else { 200 }).collect();
        let img = Image::new_2d(8, 8, data).unwrap();
        let t = otsu_threshold(&img, DEFAULT_BINS).unwrap().unwrap();
        assert!(t > 10.0 && t <= 200.0, "threshold {t}");
        let mask = apply_threshold(&img, t);
        assert_eq!(mask.data().iter().filter(|&&b| b).count(), 32);
        assert!(!mask.data()[0] && mask.data()[63]);
    }

    #[test]
    fn constant_image_has_no_threshold() {
        let img = Image::new_2d(2, 2, vec![7u16; 4]).unwrap();
        assert_eq!(otsu_threshold(&img, 256).unwrap(), None);
    }

    #[test]
    fn nan_is_rejected() {
        let img = Image::new_2d(2, 1, vec![0.0f32, f32::NAN]).unwrap();
        assert_eq!(otsu_threshold(&img, 16), Err(ImageError::NonFinite));
    }

    #[test]
    fn tiles_use_their_own_threshold() {
        // Left tile is dark (values 0 and 10), right tile is bright (100 and 200).
        // Each tile gets its own threshold, so the bright tile's background (100)
        // is separated from its foreground (200) independently of the dark tile.
        let row = [0u8, 10, 100, 200];
        let data: Vec<u8> = (0..4).flat_map(|_| row).collect();
        let img = Image::new_2d(4, 4, data).unwrap();
        let map = local_otsu_tiles(&img, 2, 256).unwrap();
        // Tile at x 0..2 threshold lies in (0, 10]; tile at x 2..4 lies in (100, 200].
        assert!(map.get(0, 0, 0).unwrap() > &0.0 && map.get(0, 0, 0).unwrap() <= &10.0);
        assert!(map.get(2, 0, 0).unwrap() > &100.0 && map.get(2, 0, 0).unwrap() <= &200.0);
        let mask = apply_threshold_map(&img, &map).unwrap();
        assert_eq!(mask.data().iter().filter(|&&b| b).count(), 8);
        assert!(!mask.get(0, 0, 0).unwrap() && *mask.get(3, 0, 0).unwrap());
    }

    #[test]
    fn sliding_window_separates_a_step_everywhere() {
        // Noise-free step: columns 0..4 are 10, columns 4..8 are 200. Windows that
        // straddle the step get their own Otsu split; windows wholly inside one level
        // fall back to the global threshold, which lies in (10, 200]. Every pixel must
        // be classed by its level.
        let data: Vec<u8> = (0..64).map(|i| if i % 8 < 4 { 10 } else { 200 }).collect();
        let img = Image::new_2d(8, 8, data.clone()).unwrap();
        for radius in [1, 2, 3, 10] {
            let map = local_otsu_window(&img, radius, 256).unwrap();
            let mask = apply_threshold_map(&img, &map).unwrap();
            let expected: Vec<bool> = data.iter().map(|&v| v == 200).collect();
            assert_eq!(mask.data(), &expected[..], "radius {radius}");
        }
    }

    #[test]
    fn sliding_window_threshold_is_exact_for_straddling_window() {
        // A 3x1 image [10, 10, 200] with radius 1: the window at x=1 holds {10, 10, 200}.
        // Histogram over [10, 200] with 2 bins: scale = 2/190; bin(10) = 0, bin(200) = 1.
        // The split after bin 0 is the only candidate, so the threshold is
        // lo + (0 + 1) / scale = 10 + 95 = 105.
        let img = Image::new_2d(3, 1, vec![10.0f64, 10.0, 200.0]).unwrap();
        let map = local_otsu_window(&img, 1, 2).unwrap();
        assert!((map.data()[1] - 105.0).abs() < 1e-9, "{}", map.data()[1]);
    }

    #[test]
    fn constant_tile_is_background() {
        let img = Image::new_2d(2, 2, vec![5u8; 4]).unwrap();
        let map = local_otsu_tiles(&img, 2, 256).unwrap();
        assert!(map.data().iter().all(|t| t.is_infinite()));
    }
}
