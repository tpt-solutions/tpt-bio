//! A simple dense image container for 2D and 3D analysis.
//!
//! Samples are stored row-major with `x` fastest, then `y`, then `z`:
//! the sample at `(x, y, z)` lives at `(z * height + y) * width + x`. A 2D image
//! is an image with `depth == 1`.

use std::fmt;

/// Errors from image construction and from operations on mismatched images.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum ImageError {
    /// `data.len()` does not equal `width * height * depth`.
    SizeMismatch {
        /// Expected number of samples.
        expected: usize,
        /// Number of samples supplied.
        actual: usize,
    },
    /// Two images that must share a shape do not.
    ShapeMismatch {
        /// Shape `(width, height, depth)` of the first image.
        left: (usize, usize, usize),
        /// Shape `(width, height, depth)` of the second image.
        right: (usize, usize, usize),
    },
    /// An operation that is defined only for 2D images got a volume.
    NotTwoDimensional {
        /// Depth of the offending image.
        depth: usize,
    },
    /// A floating-point image contains NaN, which has no ordering.
    NonFinite,
}

impl fmt::Display for ImageError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::SizeMismatch { expected, actual } => {
                write!(f, "image has {actual} samples, expected {expected}")
            }
            Self::ShapeMismatch { left, right } => {
                write!(f, "image shapes differ: {left:?} vs {right:?}")
            }
            Self::NotTwoDimensional { depth } => {
                write!(f, "operation needs a 2D image, got depth {depth}")
            }
            Self::NonFinite => f.write_str("image contains NaN"),
        }
    }
}

impl std::error::Error for ImageError {}

/// A dense, row-major image with `x` fastest, then `y`, then `z`.
#[derive(Debug, Clone, PartialEq)]
pub struct Image<T> {
    width: usize,
    height: usize,
    depth: usize,
    data: Vec<T>,
}

impl<T> Image<T> {
    /// Wraps a 2D row-major buffer. The image has depth 1.
    ///
    /// # Errors
    /// Returns [`ImageError::SizeMismatch`] if `data.len() != width * height`.
    pub fn new_2d(width: usize, height: usize, data: Vec<T>) -> Result<Self, ImageError> {
        Self::new_3d(width, height, 1, data)
    }

    /// Wraps a 3D buffer laid out with `x` fastest, then `y`, then `z`.
    ///
    /// # Errors
    /// Returns [`ImageError::SizeMismatch`] if `data.len() != width * height * depth`.
    pub fn new_3d(
        width: usize,
        height: usize,
        depth: usize,
        data: Vec<T>,
    ) -> Result<Self, ImageError> {
        let expected = width * height * depth;
        if data.len() != expected {
            return Err(ImageError::SizeMismatch { expected, actual: data.len() });
        }
        Ok(Self { width, height, depth, data })
    }

    /// Image width (`x` extent).
    #[must_use]
    pub fn width(&self) -> usize {
        self.width
    }

    /// Image height (`y` extent).
    #[must_use]
    pub fn height(&self) -> usize {
        self.height
    }

    /// Image depth (`z` extent). Equal to 1 for 2D images.
    #[must_use]
    pub fn depth(&self) -> usize {
        self.depth
    }

    /// Shape as `(width, height, depth)`.
    #[must_use]
    pub fn shape(&self) -> (usize, usize, usize) {
        (self.width, self.height, self.depth)
    }

    /// Number of samples.
    #[must_use]
    pub fn len(&self) -> usize {
        self.data.len()
    }

    /// True when the image has no samples.
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.data.is_empty()
    }

    /// The samples in row-major order.
    #[must_use]
    pub fn data(&self) -> &[T] {
        &self.data
    }

    /// Mutable access to the samples in row-major order.
    #[must_use]
    pub fn data_mut(&mut self) -> &mut [T] {
        &mut self.data
    }

    /// Consumes the image and returns its sample buffer.
    #[must_use]
    pub fn into_data(self) -> Vec<T> {
        self.data
    }

    /// Linear index of `(x, y, z)`. Does not bounds-check the coordinates.
    #[must_use]
    pub fn index(&self, x: usize, y: usize, z: usize) -> usize {
        (z * self.height + y) * self.width + x
    }

    /// Sample at `(x, y, z)`, or `None` when the coordinates are out of range.
    #[must_use]
    pub fn get(&self, x: usize, y: usize, z: usize) -> Option<&T> {
        if x < self.width && y < self.height && z < self.depth {
            self.data.get(self.index(x, y, z))
        } else {
            None
        }
    }

    /// Converts a linear index back to `(x, y, z)`.
    #[must_use]
    pub fn coords(&self, idx: usize) -> (usize, usize, usize) {
        let plane = self.width * self.height;
        let z = idx / plane;
        let rem = idx % plane;
        (rem % self.width, rem / self.width, z)
    }

    /// Index of the neighbour at `offset` from `idx`, or `None` if it falls outside.
    #[must_use]
    pub(crate) fn offset(&self, idx: usize, offset: (isize, isize, isize)) -> Option<usize> {
        let (x, y, z) = self.coords(idx);
        let nx = x.checked_add_signed(offset.0)?;
        let ny = y.checked_add_signed(offset.1)?;
        let nz = z.checked_add_signed(offset.2)?;
        if nx < self.width && ny < self.height && nz < self.depth {
            Some(self.index(nx, ny, nz))
        } else {
            None
        }
    }

    /// Applies `f` to every sample, producing an image of the same shape.
    #[must_use]
    pub fn map<U>(&self, f: impl FnMut(&T) -> U) -> Image<U> {
        Image {
            width: self.width,
            height: self.height,
            depth: self.depth,
            data: self.data.iter().map(f).collect(),
        }
    }

    /// Checks that `other` has the same shape as `self`.
    ///
    /// # Errors
    /// Returns [`ImageError::ShapeMismatch`] when the shapes differ.
    pub fn check_same_shape<U>(&self, other: &Image<U>) -> Result<(), ImageError> {
        if self.shape() == other.shape() {
            Ok(())
        } else {
            Err(ImageError::ShapeMismatch { left: self.shape(), right: other.shape() })
        }
    }
}

/// Face (6) neighbourhood in 3D; also the 4-neighbourhood when `dz` is ignored.
pub(crate) const FACE_OFFSETS: [(isize, isize, isize); 6] =
    [(1, 0, 0), (-1, 0, 0), (0, 1, 0), (0, -1, 0), (0, 0, 1), (0, 0, -1)];

/// Full 26-neighbourhood in 3D; its `dz == 0` slice is the 8-neighbourhood in 2D.
pub(crate) const FULL_OFFSETS: [(isize, isize, isize); 26] = {
    let mut out = [(0isize, 0isize, 0isize); 26];
    let mut n = 0;
    let mut dz = -1isize;
    while dz <= 1 {
        let mut dy = -1isize;
        while dy <= 1 {
            let mut dx = -1isize;
            while dx <= 1 {
                if dx != 0 || dy != 0 || dz != 0 {
                    out[n] = (dx, dy, dz);
                    n += 1;
                }
                dx += 1;
            }
            dy += 1;
        }
        dz += 1;
    }
    out
};

/// The 2D 4-neighbourhood.
pub(crate) const OFFSETS_4: [(isize, isize, isize); 4] =
    [(1, 0, 0), (-1, 0, 0), (0, 1, 0), (0, -1, 0)];

/// The 2D 8-neighbourhood.
pub(crate) const OFFSETS_8: [(isize, isize, isize); 8] =
    [(1, 0, 0), (-1, 0, 0), (0, 1, 0), (0, -1, 0), (1, 1, 0), (-1, 1, 0), (1, -1, 0), (-1, -1, 0)];

#[allow(clippy::float_cmp)]
#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn indexing_round_trips() {
        let img = Image::new_3d(4, 3, 2, (0..24).collect::<Vec<u32>>()).unwrap();
        assert_eq!(img.shape(), (4, 3, 2));
        // (z * height + y) * width + x = (1 * 3 + 1) * 4 + 2
        assert_eq!(*img.get(2, 1, 1).unwrap(), 18);
        assert_eq!(img.get(4, 0, 0), None);
        for idx in 0..img.len() {
            let (x, y, z) = img.coords(idx);
            assert_eq!(img.index(x, y, z), idx);
        }
    }

    #[test]
    fn rejects_wrong_length() {
        assert_eq!(
            Image::new_2d(3, 2, vec![0u8; 5]),
            Err(ImageError::SizeMismatch { expected: 6, actual: 5 })
        );
    }

    #[test]
    fn offset_stops_at_borders() {
        let img = Image::new_2d(3, 3, vec![0u8; 9]).unwrap();
        assert_eq!(img.offset(0, (-1, 0, 0)), None);
        assert_eq!(img.offset(0, (1, 1, 0)), Some(4));
        assert_eq!(img.offset(8, (1, 0, 0)), None);
    }

    #[test]
    fn map_keeps_shape() {
        let img = Image::new_2d(2, 2, vec![1u8, 2, 3, 4]).unwrap();
        let doubled = img.map(|&v| f64::from(v) * 2.0);
        assert_eq!(doubled.data(), &[2.0, 4.0, 6.0, 8.0]);
        assert_eq!(doubled.shape(), (2, 2, 1));
    }
}
