//! Neighbourhood and colocalisation features, per object.
//!
//! - [`neighbour_features`]: for each object, the number of distinct other objects that
//!   touch it through a 4-neighbour, and the distance from its centroid to the nearest
//!   other object's centroid (`None` when it is the only object).
//! - [`colocalization`]: Pearson correlation and Manders' coefficients between two
//!   channels, computed over each object's pixels. Manders' `M1` is the fraction of
//!   channel A's intensity lying where channel B reaches its threshold, and `M2` is the
//!   converse. Thresholds are inclusive (`value >= threshold`), as in
//!   [`crate::threshold::apply_threshold`]. A coefficient with a zero denominator is 0.

use std::collections::HashSet;

use crate::features::measure;
use crate::image::{Image, ImageError, OFFSETS_4};

/// Neighbourhood features of one 2D object.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct NeighbourFeatures {
    /// Consecutive object id, `1..=count`.
    pub label: u32,
    /// Number of distinct other objects with a pixel 4-adjacent to this object.
    pub touching: usize,
    /// Distance from this object's centroid to the nearest other object's centroid,
    /// in pixels. `None` if there is no other object.
    pub nearest_centroid_distance: Option<f64>,
}

/// Colocalisation features of one object.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ColocalizationFeatures {
    /// Consecutive object id, `1..=count`.
    pub label: u32,
    /// Pearson correlation of the two channels over the object's pixels, in `[-1, 1]`.
    /// 0 when either channel is constant within the object.
    pub pearson: f64,
    /// Manders' `M1`: the fraction of channel A's intensity at pixels where channel B is
    /// at or above its threshold.
    pub manders_m1: f64,
    /// Manders' `M2`: the fraction of channel B's intensity at pixels where channel A is
    /// at or above its threshold.
    pub manders_m2: f64,
}

/// Neighbour features of every object in a 2D label image.
///
/// Label `0` is background. Objects are identified by their label value, and the
/// results are in consecutive id order (see [`crate::features::measure`]).
///
/// # Errors
/// Returns [`ImageError::NotTwoDimensional`] for volumes.
pub fn neighbour_features(labels: &Image<u32>) -> Result<Vec<NeighbourFeatures>, ImageError> {
    require_2d(labels)?;
    let w = labels.width();
    let raw = labels.data();
    let centroids: Vec<(f64, f64)> = measure(labels, |_, px| {
        let n = px.len() as f64;
        let (sx, sy) =
            px.iter().fold((0.0, 0.0), |(sx, sy), &p| (sx + (p % w) as f64, sy + (p / w) as f64));
        (sx / n, sy / n)
    });
    let touching: Vec<usize> = measure(labels, |_, px| {
        let mut others = HashSet::new();
        for &p in px {
            let own = raw[p];
            for &d in &OFFSETS_4 {
                if let Some(n) = labels.offset(p, d) {
                    let other = raw[n];
                    if other != 0 && other != own {
                        others.insert(other);
                    }
                }
            }
        }
        others.len()
    });
    Ok(centroids
        .iter()
        .enumerate()
        .map(|(i, &(cx, cy))| {
            let nearest = centroids
                .iter()
                .enumerate()
                .filter(|&(j, _)| j != i)
                .map(|(_, &(x, y))| ((x - cx).powi(2) + (y - cy).powi(2)).sqrt())
                .reduce(f64::min);
            NeighbourFeatures {
                label: i as u32 + 1,
                touching: touching[i],
                nearest_centroid_distance: nearest,
            }
        })
        .collect())
}

/// Colocalisation features of every object, between channels `channel_a` and `channel_b`.
///
/// `threshold_a` and `threshold_b` are the Manders thresholds for each channel. Pass 0
/// for both to count every pixel of the other channel.
///
/// # Errors
/// Returns [`ImageError::ShapeMismatch`] if a channel and `labels` differ in shape.
pub fn colocalization(
    labels: &Image<u32>,
    channel_a: &Image<f64>,
    channel_b: &Image<f64>,
    threshold_a: f64,
    threshold_b: f64,
) -> Result<Vec<ColocalizationFeatures>, ImageError> {
    labels.check_same_shape(channel_a)?;
    labels.check_same_shape(channel_b)?;
    let (va, vb) = (channel_a.data(), channel_b.data());
    Ok(measure(labels, |label, px| {
        let n = px.len() as f64;
        let (mut sum_a, mut sum_b) = (0.0f64, 0.0f64);
        for &p in px {
            sum_a += va[p];
            sum_b += vb[p];
        }
        let (mean_a, mean_b) = (sum_a / n, sum_b / n);
        let (mut cov, mut var_a, mut var_b) = (0.0f64, 0.0f64, 0.0f64);
        let (mut m1_num, mut m2_num) = (0.0f64, 0.0f64);
        for &p in px {
            let da = va[p] - mean_a;
            let db = vb[p] - mean_b;
            cov += da * db;
            var_a += da * da;
            var_b += db * db;
            if vb[p] >= threshold_b {
                m1_num += va[p];
            }
            if va[p] >= threshold_a {
                m2_num += vb[p];
            }
        }
        let pearson = if var_a > 0.0 && var_b > 0.0 { cov / (var_a * var_b).sqrt() } else { 0.0 };
        ColocalizationFeatures {
            label,
            pearson,
            manders_m1: if sum_a > 0.0 { m1_num / sum_a } else { 0.0 },
            manders_m2: if sum_b > 0.0 { m2_num / sum_b } else { 0.0 },
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
    fn three_objects_in_a_row_with_known_spacing() {
        // Single-pixel objects at x = 0, 10 and 25. Nearest-centroid distances are
        // 10 (0 to 10), 10 (10 to 0) and 15 (25 to 10). None touch.
        let mut data = vec![0u32; 26];
        data[0] = 1;
        data[10] = 2;
        data[25] = 3;
        let labels = Image::new_2d(26, 1, data).unwrap();
        let f = neighbour_features(&labels).unwrap();
        assert_eq!(f.len(), 3);
        approx(f[0].nearest_centroid_distance.unwrap(), 10.0);
        approx(f[1].nearest_centroid_distance.unwrap(), 10.0);
        approx(f[2].nearest_centroid_distance.unwrap(), 15.0);
        assert!(f.iter().all(|o| o.touching == 0));
    }

    #[test]
    fn touching_counts_distinct_neighbours() {
        // Row [1, 1, 2, 2, 0, 3]. Object 1 (centroid x = 0.5) touches 2 only; object 2
        // (centroid x = 2.5) touches 1 only; object 3 (x = 5) touches nothing.
        // Nearest centroids: 1 -> 2 at 2.0; 2 -> 1 at 2.0; 3 -> 2 at 2.5.
        let labels = Image::new_2d(6, 1, vec![1u32, 1, 2, 2, 0, 3]).unwrap();
        let f = neighbour_features(&labels).unwrap();
        assert_eq!(f.iter().map(|o| o.touching).collect::<Vec<_>>(), vec![1, 1, 0]);
        approx(f[0].nearest_centroid_distance.unwrap(), 2.0);
        approx(f[1].nearest_centroid_distance.unwrap(), 2.0);
        approx(f[2].nearest_centroid_distance.unwrap(), 2.5);
    }

    #[test]
    fn lone_object_has_no_nearest_distance() {
        let labels = Image::new_2d(3, 1, vec![0u32, 4, 0]).unwrap();
        let f = neighbour_features(&labels).unwrap();
        assert_eq!(f.len(), 1);
        assert_eq!(f[0].nearest_centroid_distance, None);
    }

    #[test]
    fn identical_channels_have_pearson_one() {
        let labels = Image::new_2d(2, 2, vec![1u32; 4]).unwrap();
        let a = Image::new_2d(2, 2, vec![1.0, 2.0, 3.0, 4.0]).unwrap();
        let f = colocalization(&labels, &a, &a, 0.0, 0.0).unwrap();
        approx(f[0].pearson, 1.0);
        // With thresholds of 0 every pixel counts, so M1 = M2 = 1.
        approx(f[0].manders_m1, 1.0);
        approx(f[0].manders_m2, 1.0);
    }

    #[test]
    fn uncorrelated_channels_have_zero_pearson_and_known_manders() {
        // A = [1, 1, 0, 0] and B = [1, 0, 1, 0] over a 2x2 object (row-major).
        // Deviations from the means (0.5): A = [.5, .5, -.5, -.5], B = [.5, -.5, .5, -.5].
        // Products sum to 0, so Pearson = 0.
        // Thresholds 0.5 each. M1: pixels with B >= 0.5 are 0 and 2, where A = 1 and 0,
        // so M1 = 1 / sum(A) = 1 / 2 = 0.5. M2: pixels with A >= 0.5 are 0 and 1, where
        // B = 1 and 0, so M2 = 1 / sum(B) = 1 / 2 = 0.5.
        let labels = Image::new_2d(2, 2, vec![1u32; 4]).unwrap();
        let a = Image::new_2d(2, 2, vec![1.0, 1.0, 0.0, 0.0]).unwrap();
        let b = Image::new_2d(2, 2, vec![1.0, 0.0, 1.0, 0.0]).unwrap();
        let f = colocalization(&labels, &a, &b, 0.5, 0.5).unwrap();
        approx(f[0].pearson, 0.0);
        approx(f[0].manders_m1, 0.5);
        approx(f[0].manders_m2, 0.5);
    }

    #[test]
    fn colocalization_is_per_object() {
        // Object 1 is the left column, object 2 the right column. Each object has
        // identical channels, so each gets Pearson 1 independently of the other.
        let labels = Image::new_2d(2, 2, vec![1u32, 2, 1, 2]).unwrap();
        let a = Image::new_2d(2, 2, vec![1.0, 5.0, 2.0, 9.0]).unwrap();
        let f = colocalization(&labels, &a, &a, 0.0, 0.0).unwrap();
        assert_eq!(f.len(), 2);
        approx(f[0].pearson, 1.0);
        approx(f[1].pearson, 1.0);
    }
}
