//! Conversions between [`Image`] and `ndarray` arrays.
//!
//! An [`Image`] stores `x` fastest, so its buffer is the C-order layout of an
//! array indexed `[z, y, x]`. Converting to an [`Array3`] therefore reuses the
//! buffer without reordering. Converting back takes the buffer of a standard-layout
//! array, and copies only when the input is not in standard layout.

use ::ndarray::{Array2, Array3};

use super::{FormatError, FormatResult};
use crate::image::Image;

/// Converts a volume to an array indexed `[z, y, x]`.
///
/// # Errors
/// Returns [`FormatError::Invalid`] if the buffer does not match the shape, which
/// an [`Image`] built through its constructors never does.
pub fn image_to_array3<T: Clone>(img: &Image<T>) -> FormatResult<Array3<T>> {
    let (w, h, d) = img.shape();
    Array3::from_shape_vec((d, h, w), img.data().to_vec())
        .map_err(|e| FormatError::Invalid(format!("array shape: {e}")))
}

/// Converts a `[z, y, x]` array to an image.
///
/// # Errors
/// Returns [`FormatError::Image`] if the shape and buffer disagree, which cannot
/// happen for an array built by this module.
pub fn array3_to_image<T: Clone>(arr: Array3<T>) -> FormatResult<Image<T>> {
    let (d, h, w) = arr.dim();
    let data = if arr.is_standard_layout() {
        arr.into_raw_vec_and_offset().0
    } else {
        arr.as_standard_layout().into_owned().into_raw_vec_and_offset().0
    };
    Ok(Image::new_3d(w, h, d, data)?)
}

/// Converts a single 2D image to an array indexed `[y, x]`.
///
/// # Errors
/// Returns [`FormatError::Image`] with a not-2D error if the image has depth
/// greater than 1.
pub fn image_to_array2<T: Clone>(img: &Image<T>) -> FormatResult<Array2<T>> {
    if img.depth() != 1 {
        return Err(crate::image::ImageError::NotTwoDimensional { depth: img.depth() }.into());
    }
    let (w, h, _) = img.shape();
    Array2::from_shape_vec((h, w), img.data().to_vec())
        .map_err(|e| FormatError::Invalid(format!("array shape: {e}")))
}

/// Converts an array indexed `[y, x]` to a 2D image.
///
/// # Errors
/// Same as [`array3_to_image`].
pub fn array2_to_image<T: Clone>(arr: Array2<T>) -> FormatResult<Image<T>> {
    let (h, w) = arr.dim();
    let data = if arr.is_standard_layout() {
        arr.into_raw_vec_and_offset().0
    } else {
        arr.as_standard_layout().into_owned().into_raw_vec_and_offset().0
    };
    Ok(Image::new_2d(w, h, data)?)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn volume_indexing_matches_array() {
        let img = Image::new_3d(4, 3, 2, (0..24u16).collect()).unwrap();
        let arr = image_to_array3(&img).unwrap();
        assert_eq!(arr.dim(), (2, 3, 4));
        for z in 0..2 {
            for y in 0..3 {
                for x in 0..4 {
                    assert_eq!(arr[[z, y, x]], *img.get(x, y, z).unwrap());
                }
            }
        }
        assert_eq!(array3_to_image(arr).unwrap(), img);
    }

    #[test]
    fn non_standard_layout_is_converted() {
        use ::ndarray::ShapeBuilder;
        let logical = Array3::from_shape_fn((2, 3, 4), |(z, y, x)| (z * 100 + y * 10 + x) as f32);
        // The same values stored in Fortran order: the first index varies fastest.
        let mut f_buf = Vec::with_capacity(24);
        for x in 0..4 {
            for y in 0..3 {
                for z in 0..2 {
                    f_buf.push(logical[[z, y, x]]);
                }
            }
        }
        let fortran = Array3::from_shape_vec((2, 3, 4).f(), f_buf).unwrap();
        assert!(!fortran.is_standard_layout());
        let img = array3_to_image(fortran).unwrap();
        assert_eq!(img.shape(), (4, 3, 2));
        assert_eq!(img.get(3, 2, 1), Some(&123.0));
        assert_eq!(image_to_array3(&img).unwrap(), logical);
    }

    #[test]
    fn plane_round_trip() {
        let img = Image::new_2d(3, 2, vec![1u8, 2, 3, 4, 5, 6]).unwrap();
        let arr = image_to_array2(&img).unwrap();
        assert_eq!(arr[[1, 0]], 4);
        assert_eq!(array2_to_image(arr).unwrap(), img);

        let transposed =
            Array2::from_shape_vec((2, 3), vec![1u8, 2, 3, 4, 5, 6]).unwrap().reversed_axes();
        let img = array2_to_image(transposed).unwrap();
        assert_eq!(img.shape(), (2, 3, 1));
        assert_eq!(img.data(), &[1, 4, 2, 5, 3, 6]);
    }

    #[test]
    fn rejects_volume_as_plane() {
        let img = Image::new_3d(2, 2, 2, vec![0u8; 8]).unwrap();
        assert!(image_to_array2(&img).is_err());
    }
}
