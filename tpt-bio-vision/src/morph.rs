//! Morphology in 2D and 3D: erosion, dilation, opening, closing and top-hat.
//!
//! All structuring elements are symmetric and contain the origin, so erosion and
//! dilation need no reflection. Pixels outside the image are background for binary
//! operations (erosion therefore also erodes against the image border) and are
//! ignored for grayscale operations.

use crate::image::Image;

/// A symmetric structuring element given as a list of voxel offsets.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct StructuringElement {
    offsets: Vec<(isize, isize, isize)>,
}

impl StructuringElement {
    /// Offsets of the element, each `(dx, dy, dz)` relative to the origin.
    #[must_use]
    pub fn offsets(&self) -> &[(isize, isize, isize)] {
        &self.offsets
    }

    /// Number of offsets in the element.
    #[must_use]
    pub fn len(&self) -> usize {
        self.offsets.len()
    }

    /// True if the element is empty (never true for the constructors below).
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.offsets.is_empty()
    }

    /// `(2r+1) x (2r+1)` square in the plane, `dz = 0`.
    #[must_use]
    pub fn square(radius: usize) -> Self {
        Self::from_predicate(radius, false, |_, _, _| true)
    }

    /// Disk of Euclidean `radius` in the plane: all `dx² + dy² <= r²`, `dz = 0`.
    #[must_use]
    pub fn disk(radius: usize) -> Self {
        Self::from_predicate(radius, false, |dx, dy, _| dx * dx + dy * dy <= sq(radius))
    }

    /// `(2r+1)³` cube.
    #[must_use]
    pub fn cube(radius: usize) -> Self {
        Self::from_predicate(radius, true, |_, _, _| true)
    }

    /// Ball of Euclidean `radius`: all `dx² + dy² + dz² <= r²`.
    #[must_use]
    pub fn ball(radius: usize) -> Self {
        Self::from_predicate(radius, true, |dx, dy, dz| dx * dx + dy * dy + dz * dz <= sq(radius))
    }

    fn from_predicate(
        radius: usize,
        three_d: bool,
        keep: impl Fn(isize, isize, isize) -> bool,
    ) -> Self {
        let r = radius as isize;
        let rz = if three_d { r } else { 0 };
        let mut offsets = Vec::new();
        for dz in -rz..=rz {
            for dy in -r..=r {
                for dx in -r..=r {
                    if keep(dx, dy, dz) {
                        offsets.push((dx, dy, dz));
                    }
                }
            }
        }
        Self { offsets }
    }
}

fn sq(r: usize) -> isize {
    let r = r as isize;
    r * r
}

/// Binary erosion: a pixel stays foreground only if every element offset lands
/// inside the image on a foreground pixel.
#[must_use]
pub fn erode(mask: &Image<bool>, se: &StructuringElement) -> Image<bool> {
    let data = (0..mask.len())
        .map(|p| se.offsets().iter().all(|&d| mask.offset(p, d).is_some_and(|n| mask.data()[n])))
        .collect();
    rebuild(mask, data)
}

/// Binary dilation: a pixel is foreground if any element offset lands on a foreground pixel.
#[must_use]
pub fn dilate(mask: &Image<bool>, se: &StructuringElement) -> Image<bool> {
    let data = (0..mask.len())
        .map(|p| {
            mask.data()[p]
                || se.offsets().iter().any(|&d| mask.offset(p, d).is_some_and(|n| mask.data()[n]))
        })
        .collect();
    rebuild(mask, data)
}

/// Binary opening: erosion followed by dilation. Removes features smaller than the element.
#[must_use]
pub fn open_binary(mask: &Image<bool>, se: &StructuringElement) -> Image<bool> {
    dilate(&erode(mask, se), se)
}

/// Binary closing: dilation followed by erosion. Fills gaps smaller than the element.
#[must_use]
pub fn close_binary(mask: &Image<bool>, se: &StructuringElement) -> Image<bool> {
    erode(&dilate(mask, se), se)
}

/// Grayscale erosion: the minimum over the element, ignoring out-of-image offsets.
#[must_use]
pub fn erode_gray(image: &Image<f64>, se: &StructuringElement) -> Image<f64> {
    extremum(image, se, f64::min)
}

/// Grayscale dilation: the maximum over the element, ignoring out-of-image offsets.
#[must_use]
pub fn dilate_gray(image: &Image<f64>, se: &StructuringElement) -> Image<f64> {
    extremum(image, se, f64::max)
}

/// Grayscale opening (erosion then dilation): removes bright features smaller than the element.
#[must_use]
pub fn open_gray(image: &Image<f64>, se: &StructuringElement) -> Image<f64> {
    dilate_gray(&erode_gray(image, se), se)
}

/// Grayscale closing (dilation then erosion): removes dark features smaller than the element.
#[must_use]
pub fn close_gray(image: &Image<f64>, se: &StructuringElement) -> Image<f64> {
    erode_gray(&dilate_gray(image, se), se)
}

/// White top-hat: `image - open_gray(image)`. Keeps bright features smaller than the element.
#[must_use]
pub fn white_tophat(image: &Image<f64>, se: &StructuringElement) -> Image<f64> {
    let opened = open_gray(image, se);
    difference(image, &opened)
}

/// Black top-hat: `close_gray(image) - image`. Keeps dark features smaller than the element.
#[must_use]
pub fn black_tophat(image: &Image<f64>, se: &StructuringElement) -> Image<f64> {
    let closed = close_gray(image, se);
    difference(&closed, image)
}

/// Binary-mask counterpart of the opening-based top-hat: `mask AND NOT open(mask)`.
#[must_use]
pub fn binary_tophat(mask: &Image<bool>, se: &StructuringElement) -> Image<bool> {
    let opened = open_binary(mask, se);
    let data = mask.data().iter().zip(opened.data()).map(|(&a, &b)| a && !b).collect();
    rebuild(mask, data)
}

fn extremum(image: &Image<f64>, se: &StructuringElement, pick: fn(f64, f64) -> f64) -> Image<f64> {
    let data = (0..image.len())
        .map(|p| {
            se.offsets()
                .iter()
                .filter_map(|&d| image.offset(p, d))
                .map(|n| image.data()[n])
                .fold(image.data()[p], pick)
        })
        .collect();
    rebuild_f64(image, data)
}

fn difference(a: &Image<f64>, b: &Image<f64>) -> Image<f64> {
    let data = a.data().iter().zip(b.data()).map(|(x, y)| x - y).collect();
    rebuild_f64(a, data)
}

fn rebuild(like: &Image<bool>, data: Vec<bool>) -> Image<bool> {
    Image::new_3d(like.width(), like.height(), like.depth(), data).expect("shape preserved")
}

fn rebuild_f64(like: &Image<f64>, data: Vec<f64>) -> Image<f64> {
    Image::new_3d(like.width(), like.height(), like.depth(), data).expect("shape preserved")
}

#[allow(clippy::float_cmp)]
#[cfg(test)]
mod tests {
    use super::*;

    fn grid(w: usize, h: usize, rows: &[&str]) -> Image<bool> {
        let data = rows.iter().flat_map(|r| r.chars().map(|c| c == '#')).collect();
        Image::new_2d(w, h, data).unwrap()
    }

    #[test]
    fn element_sizes() {
        assert_eq!(StructuringElement::square(1).len(), 9);
        assert_eq!(StructuringElement::disk(1).len(), 5); // centre plus 4-neighbours
        assert_eq!(StructuringElement::disk(2).len(), 13);
        assert_eq!(StructuringElement::cube(1).len(), 27);
        assert_eq!(StructuringElement::ball(1).len(), 7);
        assert_eq!(StructuringElement::ball(2).len(), 33);
    }

    #[test]
    fn erosion_and_dilation_of_a_3x3_block() {
        let m = grid(5, 5, &[".....", ".###.", ".###.", ".###.", "....."]);
        let se = StructuringElement::square(1);
        let e = erode(&m, &se);
        // Only the centre of the block survives a 3x3 erosion.
        assert_eq!(e.data().iter().filter(|&&b| b).count(), 1);
        assert!(*e.get(2, 2, 0).unwrap());
        // The 3x3 block grows to x,y in 0..=4, which is the whole 5x5 image.
        let d = dilate(&m, &se);
        assert_eq!(d.data().iter().filter(|&&b| b).count(), 25);
    }

    #[test]
    fn opening_removes_thin_spur_closing_fills_hole() {
        // A 3x3 block with a single-pixel spur sticking out to the right.
        let m = grid(6, 5, &["......", ".###..", ".###..", ".###.#", "......"]);
        let se = StructuringElement::square(1);
        let opened = open_binary(&m, &se);
        assert!(!*opened.get(5, 3, 0).unwrap(), "spur should be removed");
        assert!(*opened.get(2, 2, 0).unwrap(), "block centre should survive");

        let holed = grid(5, 5, &[".....", ".###.", ".#.#.", ".###.", "....."]);
        let closed = close_binary(&holed, &StructuringElement::square(1));
        assert!(*closed.get(2, 2, 0).unwrap(), "hole should be filled");
    }

    #[test]
    fn disk_erosion_of_a_square_keeps_only_centre_at_radius_one() {
        let m = grid(5, 5, &[".....", ".###.", ".###.", ".###.", "....."]);
        let e = erode(&m, &StructuringElement::disk(1));
        // Disk(1) is a plus shape: the centre has all four neighbours in the block.
        assert!(*e.get(2, 2, 0).unwrap());
        assert_eq!(e.data().iter().filter(|&&b| b).count(), 1);
    }

    #[test]
    fn grayscale_tophat_isolates_a_bright_dot() {
        // Flat background of 2 with one pixel at 7. Opening with a 3x3 square removes
        // the dot, so the white top-hat is 5 at the dot and 0 elsewhere.
        let mut data = vec![2.0; 25];
        data[12] = 7.0;
        let img = Image::new_2d(5, 5, data).unwrap();
        let th = white_tophat(&img, &StructuringElement::square(1));
        assert_eq!(*th.get(2, 2, 0).unwrap(), 5.0);
        assert_eq!(*th.get(0, 0, 0).unwrap(), 0.0);
    }

    #[test]
    fn grayscale_erode_dilate_on_a_ramp() {
        let img = Image::new_2d(5, 1, vec![0.0, 1.0, 2.0, 3.0, 4.0]).unwrap();
        let se = StructuringElement::square(1);
        assert_eq!(erode_gray(&img, &se).data(), &[0.0, 0.0, 1.0, 2.0, 3.0]);
        assert_eq!(dilate_gray(&img, &se).data(), &[1.0, 2.0, 3.0, 4.0, 4.0]);
    }

    #[test]
    fn three_d_ball_dilation_reaches_face_neighbour_only_at_radius_one() {
        let mut data = vec![false; 27];
        data[13] = true; // centre of a 3x3x3 volume
        let m = Image::new_3d(3, 3, 3, data).unwrap();
        let d = dilate(&m, &StructuringElement::ball(1));
        assert_eq!(d.data().iter().filter(|&&b| b).count(), 7);
        let c = dilate(&m, &StructuringElement::cube(1));
        assert_eq!(c.data().iter().filter(|&&b| b).count(), 27);
    }
}
