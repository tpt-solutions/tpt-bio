//! Per-object feature extraction over a label image.
//!
//! [`measure`] is the framework: it groups pixels by object and calls a closure once per
//! object. [`shape_features_2d`] and [`intensity_features`] are built on it. Objects are
//! returned in order of their consecutive ids `1..=count`, as from
//! [`crate::label::relabel_sequential`]. Processing is serial; this crate has no
//! parallelism dependency.
//!
//! Definitions (2D, pixel coordinates are pixel centres):
//! - `area`: number of pixels.
//! - `perimeter`: length of the object boundary in pixel-edge units. Each side of an
//!   object pixel that faces background or the image exterior counts 1. This is not
//!   `CellProfiler`'s contour-based perimeter.
//! - `eccentricity`: `sqrt(1 - l2 / l1)` from the eigenvalues `l1 >= l2` of the
//!   second central moment matrix. 0 for a disk, 1 for a line.
//! - `solidity`: area divided by the area of the convex hull of the object's pixel
//!   corners. Always in `(0, 1]`.
//! - Intensity statistics use population standard deviation and linear-interpolation
//!   quantiles (the same convention as `NumPy`'s default).
//! - `mad`: median of `|x - median|`, with the same quantile convention.
//! - `edge_intensity`: mean intensity over boundary pixels, where a boundary pixel is an
//!   object pixel with a background or out-of-image 4-neighbour (in-plane for volumes).

use crate::image::{Image, ImageError, OFFSETS_4};
use crate::label::{object_pixels, relabel_sequential};

/// Shape features of one 2D object.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ShapeFeatures {
    /// Consecutive object id, `1..=count`.
    pub label: u32,
    /// Area in pixels.
    pub area: usize,
    /// Boundary length in pixel-edge units.
    pub perimeter: f64,
    /// Eccentricity in `[0, 1]`.
    pub eccentricity: f64,
    /// Area over convex-hull area, in `(0, 1]`.
    pub solidity: f64,
}

/// Intensity features of one object.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct IntensityFeatures {
    /// Consecutive object id, `1..=count`.
    pub label: u32,
    /// Mean intensity.
    pub mean: f64,
    /// Population standard deviation.
    pub std: f64,
    /// Minimum intensity.
    pub min: f64,
    /// Maximum intensity.
    pub max: f64,
    /// 25th percentile.
    pub q25: f64,
    /// Median.
    pub median: f64,
    /// 75th percentile.
    pub q75: f64,
    /// Sum of intensities (integrated intensity).
    pub integrated: f64,
    /// Median absolute deviation: the median of `|x - median|`.
    pub mad: f64,
    /// Mean intensity of the boundary pixels. A boundary pixel is an object pixel with
    /// a 4-neighbour that is background (label 0) or outside the image. Zero if the
    /// object has no boundary pixels.
    pub edge_intensity: f64,
}

/// Calls `f(label, pixel_indices)` for each object, in consecutive-id order.
///
/// Labels in `labels` need not be consecutive; they are normalised first.
pub fn measure<R>(labels: &Image<u32>, mut f: impl FnMut(u32, &[usize]) -> R) -> Vec<R> {
    let normalized = relabel_sequential(labels);
    let groups = object_pixels(&normalized.labels, normalized.count);
    groups.iter().enumerate().map(|(i, px)| f(i as u32 + 1, px)).collect()
}

/// Shape features of every object in a 2D label image.
///
/// # Errors
/// Returns [`ImageError::NotTwoDimensional`] for volumes.
pub fn shape_features_2d(labels: &Image<u32>) -> Result<Vec<ShapeFeatures>, ImageError> {
    if labels.depth() != 1 {
        return Err(ImageError::NotTwoDimensional { depth: labels.depth() });
    }
    let normalized = relabel_sequential(labels);
    let lab = normalized.labels.data();
    let (w, h) = (labels.width(), labels.height());
    Ok(measure(labels, |label, px| {
        let area = px.len();
        // Perimeter: count pixel sides that face a different label or the exterior.
        let mut perimeter = 0usize;
        for &p in px {
            let (x, y) = (p % w, p / w);
            let same = |nx: isize, ny: isize| {
                nx >= 0
                    && ny >= 0
                    && (nx as usize) < w
                    && (ny as usize) < h
                    && lab[ny as usize * w + nx as usize] == label
            };
            for (dx, dy) in [(1isize, 0isize), (-1, 0), (0, 1), (0, -1)] {
                if !same(x as isize + dx, y as isize + dy) {
                    perimeter += 1;
                }
            }
        }
        let eccentricity = eccentricity(px, w);
        let solidity = if area == 0 { 0.0 } else { area as f64 / hull_area(px, w) };
        ShapeFeatures { label, area, perimeter: perimeter as f64, eccentricity, solidity }
    }))
}

/// Eccentricity from the eigenvalues of the second central moment matrix.
fn eccentricity(px: &[usize], w: usize) -> f64 {
    let n = px.len() as f64;
    if px.is_empty() {
        return 0.0;
    }
    let (mut cx, mut cy) = (0.0, 0.0);
    for &p in px {
        cx += (p % w) as f64;
        cy += (p / w) as f64;
    }
    cx /= n;
    cy /= n;
    let (mut mu20, mut mu02, mut mu11) = (0.0, 0.0, 0.0);
    for &p in px {
        let dx = (p % w) as f64 - cx;
        let dy = (p / w) as f64 - cy;
        mu20 += dx * dx;
        mu02 += dy * dy;
        mu11 += dx * dy;
    }
    mu20 /= n;
    mu02 /= n;
    mu11 /= n;
    let half_trace = f64::midpoint(mu20, mu02);
    let det = mu20 * mu02 - mu11 * mu11;
    let disc = (half_trace * half_trace - det).max(0.0).sqrt();
    let l1 = half_trace + disc;
    let l2 = (half_trace - disc).max(0.0);
    if l1 > 0.0 { (1.0 - l2 / l1).max(0.0).sqrt() } else { 0.0 }
}

/// Area of the convex hull of the pixel corners of `px` (each pixel is a unit square).
fn hull_area(px: &[usize], w: usize) -> f64 {
    let mut pts: Vec<(i64, i64)> = Vec::with_capacity(px.len() * 4);
    for &p in px {
        let (x, y) = ((p % w) as i64, (p / w) as i64);
        pts.extend([(x, y), (x + 1, y), (x, y + 1), (x + 1, y + 1)]);
    }
    pts.sort_unstable();
    pts.dedup();
    if pts.len() < 3 {
        return 0.0;
    }
    let cross = |o: (i64, i64), a: (i64, i64), b: (i64, i64)| {
        (a.0 - o.0) * (b.1 - o.1) - (a.1 - o.1) * (b.0 - o.0)
    };
    let mut hull: Vec<(i64, i64)> = Vec::with_capacity(pts.len() * 2);
    for &p in &pts {
        while hull.len() >= 2 && cross(hull[hull.len() - 2], hull[hull.len() - 1], p) <= 0 {
            hull.pop();
        }
        hull.push(p);
    }
    let lower_len = hull.len() + 1;
    for &p in pts.iter().rev().skip(1) {
        while hull.len() >= lower_len && cross(hull[hull.len() - 2], hull[hull.len() - 1], p) <= 0 {
            hull.pop();
        }
        hull.push(p);
    }
    hull.pop();
    // Shoelace formula over the counter-clockwise hull.
    let mut twice = 0i64;
    for i in 0..hull.len() {
        let a = hull[i];
        let b = hull[(i + 1) % hull.len()];
        twice += a.0 * b.1 - b.0 * a.1;
    }
    (twice.abs() as f64) / 2.0
}

/// Intensity features of every object, using `intensity` sampled at the label pixels.
///
/// # Errors
/// Returns [`ImageError::ShapeMismatch`] if `labels` and `intensity` differ in shape.
pub fn intensity_features(
    labels: &Image<u32>,
    intensity: &Image<f64>,
) -> Result<Vec<IntensityFeatures>, ImageError> {
    labels.check_same_shape(intensity)?;
    let vals = intensity.data();
    let raw = labels.data();
    Ok(measure(labels, |label, px| {
        let mut v: Vec<f64> = px.iter().map(|&p| vals[p]).collect();
        v.sort_by(f64::total_cmp);
        let n = v.len();
        if n == 0 {
            return IntensityFeatures {
                label,
                mean: 0.0,
                std: 0.0,
                min: 0.0,
                max: 0.0,
                q25: 0.0,
                median: 0.0,
                q75: 0.0,
                integrated: 0.0,
                mad: 0.0,
                edge_intensity: 0.0,
            };
        }
        let sum: f64 = v.iter().sum();
        let mean = sum / n as f64;
        let var = v.iter().map(|x| (x - mean) * (x - mean)).sum::<f64>() / n as f64;
        let median = quantile(&v, 0.5);
        let mut deviations: Vec<f64> = v.iter().map(|x| (x - median).abs()).collect();
        deviations.sort_by(f64::total_cmp);
        let mad = quantile(&deviations, 0.5);

        let (mut edge_sum, mut edge_n) = (0.0f64, 0usize);
        for &p in px {
            let is_edge = OFFSETS_4.iter().any(|&d| match labels.offset(p, d) {
                None => true,
                Some(q) => raw[q] == 0,
            });
            if is_edge {
                edge_sum += vals[p];
                edge_n += 1;
            }
        }
        let edge_intensity = if edge_n == 0 { 0.0 } else { edge_sum / edge_n as f64 };

        IntensityFeatures {
            label,
            mean,
            std: var.sqrt(),
            min: v[0],
            max: v[n - 1],
            q25: quantile(&v, 0.25),
            median,
            q75: quantile(&v, 0.75),
            integrated: sum,
            mad,
            edge_intensity,
        }
    }))
}

/// Linear-interpolation quantile of sorted, non-empty `v`.
fn quantile(v: &[f64], q: f64) -> f64 {
    let pos = q * (v.len() - 1) as f64;
    let lo = pos.floor() as usize;
    let hi = pos.ceil() as usize;
    let frac = pos - lo as f64;
    v[lo] + (v[hi] - v[lo]) * frac
}

#[allow(clippy::float_cmp)]
#[cfg(test)]
mod tests {
    use super::*;

    fn approx(a: f64, b: f64) {
        assert!((a - b).abs() < 1e-9, "{a} != {b}");
    }

    /// Builds a 2D label image from rows of label ids.
    fn labels(rows: &[&[u32]]) -> Image<u32> {
        let w = rows[0].len();
        let data = rows.iter().flat_map(|r| r.iter().copied()).collect();
        Image::new_2d(w, rows.len(), data).unwrap()
    }

    #[test]
    fn square_three_by_three() {
        let img = labels(&[&[1, 1, 1], &[1, 1, 1], &[1, 1, 1]]);
        let f = shape_features_2d(&img).unwrap();
        assert_eq!(f.len(), 1);
        assert_eq!(f[0].area, 9);
        approx(f[0].perimeter, 12.0);
        approx(f[0].eccentricity, 0.0);
        approx(f[0].solidity, 1.0);
    }

    #[test]
    fn horizontal_domino_is_a_line() {
        let img = labels(&[&[1, 1]]);
        let f = shape_features_2d(&img).unwrap();
        approx(f[0].perimeter, 6.0);
        approx(f[0].eccentricity, 1.0);
        approx(f[0].solidity, 1.0);
    }

    #[test]
    fn l_tromino_hand_values() {
        // Pixels (0,0), (1,0), (0,1). Verified by hand:
        // perimeter = 3*4 - 2*2 adjacencies = 8; eigenvalues 1/3 and 1/9;
        // eccentricity = sqrt(1 - (1/9)/(1/3)) = sqrt(2/3);
        // hull of corners is (0,0),(2,0),(2,1),(1,2),(0,2) with area 3.5,
        // so solidity = 3 / 3.5 = 6/7.
        let img = labels(&[&[1, 1], &[1, 0]]);
        let f = shape_features_2d(&img).unwrap();
        assert_eq!(f[0].area, 3);
        approx(f[0].perimeter, 8.0);
        approx(f[0].eccentricity, (2.0f64 / 3.0).sqrt());
        approx(f[0].solidity, 6.0 / 7.0);
    }

    #[test]
    fn two_objects_are_ordered_and_separate() {
        let img = labels(&[&[7, 0, 3], &[7, 0, 3]]);
        let f = shape_features_2d(&img).unwrap();
        assert_eq!(f.len(), 2);
        assert_eq!((f[0].label, f[0].area), (1, 2));
        assert_eq!((f[1].label, f[1].area), (2, 2));
        approx(f[0].perimeter, 6.0);
    }

    #[test]
    fn intensity_statistics_hand_values() {
        let lab = labels(&[&[1, 1, 1, 1, 0]]);
        let int = Image::new_2d(5, 1, vec![1.0, 2.0, 3.0, 4.0, 100.0]).unwrap();
        let f = intensity_features(&lab, &int).unwrap();
        assert_eq!(f.len(), 1);
        let s = f[0];
        approx(s.mean, 2.5);
        approx(s.std, 1.25f64.sqrt());
        approx(s.min, 1.0);
        approx(s.max, 4.0);
        approx(s.median, 2.5);
        approx(s.q25, 1.75);
        approx(s.q75, 3.25);
        approx(s.integrated, 10.0);
    }

    #[test]
    fn mad_and_edge_intensity_on_three_by_three_pattern() {
        // Values 1..=9 in row-major order. Median is 5, absolute deviations are
        // 4,3,2,1,0,1,2,3,4, whose median is 2, so MAD = 2.
        // Only the centre (value 5) has four object neighbours; the eight others
        // touch the image exterior, so edge mean = (1+2+3+4+6+7+8+9)/8 = 5.
        let lab = labels(&[&[1, 1, 1], &[1, 1, 1], &[1, 1, 1]]);
        let int = Image::new_2d(3, 3, (1..=9).map(f64::from).collect()).unwrap();
        let f = intensity_features(&lab, &int).unwrap();
        approx(f[0].median, 5.0);
        approx(f[0].mad, 2.0);
        approx(f[0].edge_intensity, 5.0);
        approx(f[0].integrated, 45.0);
    }

    #[test]
    fn edge_pixels_include_background_neighbours() {
        // Object 1 is the two left pixels of a row. The pixel at x=1 touches background
        // (x=2), and the pixel at x=0 touches the exterior, so both are edge pixels:
        // edge mean = (10 + 30) / 2 = 20.
        let lab = labels(&[&[1, 1, 0]]);
        let int = Image::new_2d(3, 1, vec![10.0, 30.0, 0.0]).unwrap();
        let f = intensity_features(&lab, &int).unwrap();
        approx(f[0].edge_intensity, 20.0);
    }

    #[test]
    fn intensity_rejects_shape_mismatch() {
        let lab = labels(&[&[1, 1]]);
        let int = Image::new_2d(1, 2, vec![0.0, 0.0]).unwrap();
        assert!(matches!(intensity_features(&lab, &int), Err(ImageError::ShapeMismatch { .. })));
    }

    #[test]
    fn measure_passes_pixels_per_object() {
        let lab = labels(&[&[2, 0, 2]]);
        let counts = measure(&lab, |label, px| (label, px.len()));
        assert_eq!(counts, vec![(1, 2)]);
    }

    #[test]
    fn shape_rejects_volumes() {
        let vol = Image::new_3d(1, 1, 2, vec![1u32, 1]).unwrap();
        assert!(shape_features_2d(&vol).is_err());
    }
}
