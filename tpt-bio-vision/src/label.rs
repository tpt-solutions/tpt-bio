//! Connected-component labelling in 2D (4/8-connectivity) and 3D (6/26-connectivity).
//!
//! Labels are `u32`, with `0` as background. Components are found by breadth-first
//! flood fill, and ids are assigned in raster order of each component's first pixel,
//! so the output labels are always consecutive `1..=count`.

use std::collections::VecDeque;

use crate::image::{FACE_OFFSETS, FULL_OFFSETS, Image, ImageError, OFFSETS_4, OFFSETS_8};

/// Neighbourhood for 2D labelling.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Connectivity2 {
    /// Edge-sharing neighbours only.
    Four,
    /// Edge- and corner-sharing neighbours.
    Eight,
}

/// Neighbourhood for 3D labelling.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Connectivity3 {
    /// Face-sharing neighbours only.
    Six,
    /// Face-, edge- and corner-sharing neighbours.
    TwentySix,
}

/// A label image and the number of components it contains.
#[derive(Debug, Clone, PartialEq)]
pub struct Labeling {
    /// Label per pixel; `0` is background, components are `1..=count`.
    pub labels: Image<u32>,
    /// Number of components.
    pub count: u32,
}

/// Labels the 2D foreground of `mask` with the given connectivity.
///
/// # Errors
/// Returns [`ImageError::NotTwoDimensional`] if `mask.depth() != 1`.
pub fn label_2d(mask: &Image<bool>, connectivity: Connectivity2) -> Result<Labeling, ImageError> {
    if mask.depth() != 1 {
        return Err(ImageError::NotTwoDimensional { depth: mask.depth() });
    }
    let offsets: &[(isize, isize, isize)] = match connectivity {
        Connectivity2::Four => &OFFSETS_4,
        Connectivity2::Eight => &OFFSETS_8,
    };
    Ok(flood_label(mask, offsets))
}

/// Labels the 3D foreground of `mask` with the given connectivity.
#[must_use]
pub fn label_3d(mask: &Image<bool>, connectivity: Connectivity3) -> Labeling {
    let offsets: &[(isize, isize, isize)] = match connectivity {
        Connectivity3::Six => &FACE_OFFSETS,
        Connectivity3::TwentySix => &FULL_OFFSETS,
    };
    flood_label(mask, offsets)
}

fn flood_label(mask: &Image<bool>, offsets: &[(isize, isize, isize)]) -> Labeling {
    let mut labels = vec![0u32; mask.len()];
    let mut count = 0u32;
    let mut queue = VecDeque::new();
    for seed in 0..mask.len() {
        if !mask.data()[seed] || labels[seed] != 0 {
            continue;
        }
        count += 1;
        labels[seed] = count;
        queue.push_back(seed);
        while let Some(p) = queue.pop_front() {
            for &d in offsets {
                if let Some(n) = mask.offset(p, d)
                    && mask.data()[n]
                    && labels[n] == 0
                {
                    labels[n] = count;
                    queue.push_back(n);
                }
            }
        }
    }
    let labels = Image::new_3d(mask.width(), mask.height(), mask.depth(), labels)
        .expect("label buffer matches mask shape");
    Labeling { labels, count }
}

/// Renumbers the nonzero labels of `labels` to `1..=count` in order of first appearance
/// in raster order. Background stays `0`.
#[must_use]
pub fn relabel_sequential(labels: &Image<u32>) -> Labeling {
    let mut map = std::collections::HashMap::new();
    let mut next = 0u32;
    let data: Vec<u32> = labels
        .data()
        .iter()
        .map(|&l| {
            if l == 0 {
                0
            } else {
                *map.entry(l).or_insert_with(|| {
                    next += 1;
                    next
                })
            }
        })
        .collect();
    let labels = Image::new_3d(labels.width(), labels.height(), labels.depth(), data)
        .expect("relabel buffer matches input shape");
    Labeling { labels, count: next }
}

/// Pixel indices of each object, indexed by label: `result[l - 1]` holds the pixels of
/// label `l`, for `l` in `1..=count`. Labels above `count` are ignored.
#[must_use]
pub fn object_pixels(labels: &Image<u32>, count: u32) -> Vec<Vec<usize>> {
    let mut out = vec![Vec::new(); count as usize];
    for (idx, &l) in labels.data().iter().enumerate() {
        if l >= 1 && l <= count {
            out[l as usize - 1].push(idx);
        }
    }
    out
}

#[allow(clippy::float_cmp)]
#[cfg(test)]
mod tests {
    use super::*;

    fn mask_2d(w: usize, rows: &[&str]) -> Image<bool> {
        let data: Vec<bool> = rows.iter().flat_map(|r| r.chars().map(|c| c == '#')).collect();
        Image::new_2d(w, rows.len(), data).unwrap()
    }

    #[test]
    fn diagonal_pixels_split_under_four_and_merge_under_eight() {
        let m = mask_2d(3, &["#..", ".#.", "..#"]);
        let four = label_2d(&m, Connectivity2::Four).unwrap();
        assert_eq!(four.count, 3);
        let eight = label_2d(&m, Connectivity2::Eight).unwrap();
        assert_eq!(eight.count, 1);
    }

    #[test]
    fn labels_are_consecutive_in_raster_order() {
        let m = mask_2d(5, &["##..#", "##..#", "....."]);
        let l = label_2d(&m, Connectivity2::Four).unwrap();
        assert_eq!(l.count, 2);
        assert_eq!(l.labels.data(), &[1, 1, 0, 0, 2, 1, 1, 0, 0, 2, 0, 0, 0, 0, 0]);
    }

    #[test]
    fn u_shape_is_one_component() {
        // A U shape connects only through the bottom row; flood fill must follow it.
        let m = mask_2d(3, &["#.#", "#.#", "###"]);
        let l = label_2d(&m, Connectivity2::Four).unwrap();
        assert_eq!(l.count, 1);
    }

    #[test]
    fn rejects_volume_for_2d() {
        let m = Image::new_3d(2, 2, 2, vec![true; 8]).unwrap();
        assert_eq!(
            label_2d(&m, Connectivity2::Four),
            Err(ImageError::NotTwoDimensional { depth: 2 })
        );
    }

    #[test]
    fn three_d_connectivity_counts() {
        // Two voxels that touch only along an edge (dx=dy=1, dz=0).
        let mut data = vec![false; 8];
        data[0] = true; // (0,0,0)
        data[3] = true; // (1,1,0)
        let m = Image::new_3d(2, 2, 2, data).unwrap();
        assert_eq!(label_3d(&m, Connectivity3::Six).count, 2);
        assert_eq!(label_3d(&m, Connectivity3::TwentySix).count, 1);

        // Two voxels touching only at a corner (dx=dy=dz=1).
        let mut data = vec![false; 8];
        data[0] = true; // (0,0,0)
        data[7] = true; // (1,1,1)
        let m = Image::new_3d(2, 2, 2, data).unwrap();
        assert_eq!(label_3d(&m, Connectivity3::TwentySix).count, 1);
        assert_eq!(label_3d(&m, Connectivity3::Six).count, 2);
    }

    #[test]
    fn three_d_face_connected_column() {
        let m = Image::new_3d(1, 1, 5, vec![true, true, false, true, true]).unwrap();
        let l = label_3d(&m, Connectivity3::Six);
        assert_eq!(l.count, 2);
        assert_eq!(l.labels.data(), &[1, 1, 0, 2, 2]);
    }

    #[test]
    fn relabel_compacts_arbitrary_ids() {
        let img = Image::new_2d(4, 1, vec![0u32, 7, 7, 3]).unwrap();
        let r = relabel_sequential(&img);
        assert_eq!(r.count, 2);
        assert_eq!(r.labels.data(), &[0, 1, 1, 2]);
    }

    #[test]
    fn object_pixels_groups_by_label() {
        let img = Image::new_2d(3, 1, vec![2u32, 0, 1]).unwrap();
        let groups = object_pixels(&img, 2);
        assert_eq!(groups, vec![vec![2], vec![0]]);
    }
}
