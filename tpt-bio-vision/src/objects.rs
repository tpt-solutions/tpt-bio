//! Object filtering: remove labelled objects by size or by touching the image border,
//! then renumber the survivors to consecutive ids.

use crate::image::Image;
use crate::label::{Labeling, object_pixels, relabel_sequential};

/// Criteria for keeping an object. The default keeps everything.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub struct ObjectFilter {
    /// Smallest kept area, in voxels (inclusive).
    pub min_area: usize,
    /// Largest kept area, in voxels (inclusive). `None` means no upper limit.
    pub max_area: Option<usize>,
    /// Drop objects with any voxel on the image boundary.
    pub remove_border: bool,
}

/// Keeps the objects that satisfy `filter` and renumbers them `1..=kept`.
///
/// Input ids need not be consecutive: the input is first relabelled, so each distinct
/// nonzero value is one object. Output ids follow raster order of each object's first
/// pixel, as in [`relabel_sequential`].
#[must_use]
pub fn filter_objects(labels: &Image<u32>, filter: &ObjectFilter) -> Labeling {
    let (w, h, d) = labels.shape();
    let normalized = relabel_sequential(labels);
    let groups = object_pixels(&normalized.labels, normalized.count);
    let mut out = normalized.labels.data().to_vec();
    for pixels in &groups {
        let area = pixels.len();
        let too_small = area < filter.min_area;
        let too_large = filter.max_area.is_some_and(|m| area > m);
        let touches = filter.remove_border
            && pixels.iter().any(|&p| {
                let (x, y, z) = labels.coords(p);
                x == 0 || y == 0 || x + 1 == w || y + 1 == h || (d > 1 && (z == 0 || z + 1 == d))
            });
        if too_small || too_large || touches {
            for &p in pixels {
                out[p] = 0;
            }
        }
    }
    let kept = Image::new_3d(w, h, d, out).expect("shape preserved");
    relabel_sequential(&kept)
}

#[allow(clippy::float_cmp)]
#[cfg(test)]
mod tests {
    use super::*;

    fn labelled(w: usize, data: Vec<u32>) -> Image<u32> {
        Image::new_2d(w, data.len() / w, data).unwrap()
    }

    #[test]
    fn removes_small_and_large_objects_and_renumbers() {
        // Object 1: area 1. Object 2: area 3. Object 3: area 5.
        let img = labelled(
            5,
            vec![
                1, 0, 2, 2, 0, //
                0, 0, 2, 0, 0, //
                3, 3, 3, 3, 3, //
            ],
        );
        let f = ObjectFilter { min_area: 2, max_area: Some(4), remove_border: false };
        let out = filter_objects(&img, &f);
        assert_eq!(out.count, 1);
        assert_eq!(out.labels.data(), &[0, 0, 1, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0]);
    }

    #[test]
    fn border_touching_objects_are_removed() {
        // Object 1 sits on the left edge; object 2 is interior.
        let img = labelled(
            4,
            vec![
                1, 0, 0, 0, //
                1, 0, 2, 0, //
                0, 0, 2, 0, //
                0, 0, 0, 0, //
            ],
        );
        let f = ObjectFilter { remove_border: true, ..ObjectFilter::default() };
        let out = filter_objects(&img, &f);
        assert_eq!(out.count, 1);
        assert_eq!(*out.labels.get(2, 1, 0).unwrap(), 1);
        assert_eq!(*out.labels.get(0, 0, 0).unwrap(), 0);
    }

    #[test]
    fn default_filter_keeps_everything_and_compacts_ids() {
        let img = labelled(3, vec![0, 5, 0, 9, 0, 5]);
        let out = filter_objects(&img, &ObjectFilter::default());
        assert_eq!(out.count, 2);
        assert_eq!(out.labels.data(), &[0, 1, 0, 2, 0, 1]);
    }

    #[test]
    fn three_d_border_detection() {
        // A voxel at the centre of a 3x3x3 volume is interior; one at the corner is not.
        let mut data = vec![0u32; 27];
        data[13] = 1;
        data[0] = 2;
        let img = Image::new_3d(3, 3, 3, data).unwrap();
        let f = ObjectFilter { remove_border: true, ..ObjectFilter::default() };
        let out = filter_objects(&img, &f);
        assert_eq!(out.count, 1);
        assert_eq!(out.labels.data()[13], 1);
    }
}
