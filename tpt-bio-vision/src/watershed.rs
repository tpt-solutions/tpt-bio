//! Marker-controlled watershed by priority-queue flooding.
//!
//! Each non-zero marker label seeds a region. Regions grow outward through
//! face-adjacent pixels (4-connectivity in 2D, 6 in 3D), always expanding the
//! unlabelled neighbour with the lowest value of the input image first. Ties are
//! broken by insertion order, so the result is deterministic.

use std::cmp::Ordering;
use std::collections::BinaryHeap;

use crate::image::{FACE_OFFSETS, Image, ImageError, OFFSETS_4};

/// Heap entry ordered so that the smallest value pops first; FIFO on ties.
#[derive(Debug, Clone, Copy)]
struct Entry {
    value: f64,
    seq: u64,
    idx: usize,
}

impl PartialEq for Entry {
    fn eq(&self, other: &Self) -> bool {
        self.cmp(other) == Ordering::Equal
    }
}

impl Eq for Entry {}

impl PartialOrd for Entry {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for Entry {
    fn cmp(&self, other: &Self) -> Ordering {
        // BinaryHeap is a max-heap: reverse the natural order so the lowest value wins.
        other.value.total_cmp(&self.value).then_with(|| other.seq.cmp(&self.seq))
    }
}

/// Floods `image` from the labelled `markers`, returning a label image.
///
/// Pixels that are not reached (outside `mask`, or in a region with no marker) stay `0`.
/// Marker pixels keep their labels. `mask`, when given, restricts flooding to `true` pixels.
///
/// # Errors
/// Returns [`ImageError::ShapeMismatch`] if `markers` or `mask` differ in shape from
/// `image`, or [`ImageError::NonFinite`] if `image` contains NaN.
pub fn watershed(
    image: &Image<f64>,
    markers: &Image<u32>,
    mask: Option<&Image<bool>>,
) -> Result<Image<u32>, ImageError> {
    image.check_same_shape(markers)?;
    if let Some(m) = mask {
        image.check_same_shape(m)?;
    }
    if image.data().iter().any(|v| v.is_nan()) {
        return Err(ImageError::NonFinite);
    }
    let allowed = |i: usize| mask.is_none_or(|m| m.data()[i]);

    let mut labels = markers.data().to_vec();
    let mut heap = BinaryHeap::new();
    let mut seq = 0u64;
    // Seed the queue with every marker pixel so each one expands into its neighbours.
    for (idx, &l) in labels.iter().enumerate() {
        if l != 0 && allowed(idx) {
            heap.push(Entry { value: image.data()[idx], seq, idx });
            seq += 1;
        }
    }

    while let Some(Entry { idx: p, .. }) = heap.pop() {
        let label = labels[p];
        for &d in neighbourhood(image) {
            if let Some(n) = image.offset(p, d)
                && labels[n] == 0
                && allowed(n)
            {
                labels[n] = label;
                heap.push(Entry { value: image.data()[n], seq, idx: n });
                seq += 1;
            }
        }
    }

    Image::new_3d(image.width(), image.height(), image.depth(), labels)
}

/// Face neighbourhood for the image's dimensionality.
fn neighbourhood(image: &Image<f64>) -> &'static [(isize, isize, isize)] {
    if image.depth() == 1 { &OFFSETS_4 } else { &FACE_OFFSETS }
}

#[allow(clippy::float_cmp)]
#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn two_basins_split_at_the_ridge() {
        // A 1x7 valley profile: the ridge at x=3 has the highest value and must divide
        // the two markers' regions.
        let values = vec![0.0, 1.0, 2.0, 9.0, 2.0, 1.0, 0.0];
        let img = Image::new_2d(7, 1, values).unwrap();
        let mut m = vec![0u32; 7];
        m[0] = 1;
        m[6] = 2;
        let markers = Image::new_2d(7, 1, m).unwrap();
        let out = watershed(&img, &markers, None).unwrap();
        // The ridge pixel (value 9) is claimed by label 1 first; both fronts reach it
        // at the same priority-queue level, and the earlier insertion wins.
        assert_eq!(out.data(), &[1, 1, 1, 1, 2, 2, 2]);
    }

    #[test]
    fn pixels_are_claimed_in_order_of_value() {
        // Label 2 starts at value 3 and must wait until label 1 has expanded through
        // the values 1 and 2; by then pixel 2 (value 5) is reached from the left.
        let img = Image::new_2d(5, 1, vec![0.0, 1.0, 5.0, 2.0, 3.0]).unwrap();
        let markers = Image::new_2d(5, 1, vec![1u32, 0, 0, 0, 2]).unwrap();
        let out = watershed(&img, &markers, None).unwrap();
        assert_eq!(out.data(), &[1, 1, 1, 2, 2]);
    }

    #[test]
    fn mask_blocks_flooding() {
        let img = Image::new_2d(3, 1, vec![0.0, 0.0, 0.0]).unwrap();
        let markers = Image::new_2d(3, 1, vec![4u32, 0, 0]).unwrap();
        let mask = Image::new_2d(3, 1, vec![true, false, true]).unwrap();
        let out = watershed(&img, &markers, Some(&mask)).unwrap();
        assert_eq!(out.data(), &[4, 0, 0]);
    }

    #[test]
    fn rejects_shape_mismatch() {
        let img = Image::new_2d(2, 2, vec![0.0; 4]).unwrap();
        let markers = Image::new_2d(1, 2, vec![0u32; 2]).unwrap();
        assert!(matches!(watershed(&img, &markers, None), Err(ImageError::ShapeMismatch { .. })));
    }
}
