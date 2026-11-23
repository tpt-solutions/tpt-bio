//! Multi-channel volumes with physical voxel spacing, and trilinear sampling.

use super::RenderError;
use crate::image::Image;

/// A multi-channel 3D volume with physical voxel spacing.
///
/// Every channel is an [`Image<f32>`] of the same shape. Integer microscope data is
/// converted losslessly with [`Volume::from_u16`]. The spacing is the physical size of one
/// voxel along `x`, `y` and `z` in world units. Voxel `(i, j, k)` is centred at world position
/// `((i + 0.5) * sx, (j + 0.5) * sy, (k + 0.5) * sz)`, so the volume occupies the box
/// `[0, extent]` in world space.
#[derive(Debug, Clone, PartialEq)]
pub struct Volume {
    channels: Vec<Image<f32>>,
    spacing: [f64; 3],
}

impl Volume {
    /// Builds a volume from `f32` channels.
    ///
    /// # Errors
    /// Returns [`RenderError::Empty`] when there are no channels or no samples,
    /// [`RenderError::Image`] when the channels differ in shape, and
    /// [`RenderError::InvalidParameter`] when a spacing component is not finite and positive.
    pub fn new(channels: Vec<Image<f32>>, spacing: [f64; 3]) -> Result<Self, RenderError> {
        let first = channels.first().ok_or(RenderError::Empty)?;
        if first.is_empty() {
            return Err(RenderError::Empty);
        }
        for ch in &channels[1..] {
            first.check_same_shape(ch)?;
        }
        if spacing.iter().any(|s| !s.is_finite() || *s <= 0.0) {
            return Err(RenderError::InvalidParameter("voxel spacing must be finite and positive"));
        }
        Ok(Self { channels, spacing })
    }

    /// Builds a volume from `u16` channels, converting each sample to `f32` exactly.
    ///
    /// # Errors
    /// Same as [`Volume::new`].
    pub fn from_u16(channels: &[Image<u16>], spacing: [f64; 3]) -> Result<Self, RenderError> {
        let converted = channels.iter().map(|c| c.map(|&v| f32::from(v))).collect();
        Self::new(converted, spacing)
    }

    /// Number of channels.
    #[must_use]
    pub fn channel_count(&self) -> usize {
        self.channels.len()
    }

    /// Shape `(width, height, depth)` shared by every channel.
    #[must_use]
    pub fn shape(&self) -> (usize, usize, usize) {
        self.channels[0].shape()
    }

    /// Physical voxel spacing `[sx, sy, sz]`.
    #[must_use]
    pub fn spacing(&self) -> [f64; 3] {
        self.spacing
    }

    /// Physical size `[shape * spacing]` of the volume's bounding box.
    #[must_use]
    pub fn extent(&self) -> [f64; 3] {
        let (w, h, d) = self.shape();
        [w as f64 * self.spacing[0], h as f64 * self.spacing[1], d as f64 * self.spacing[2]]
    }

    /// The channel at `index`, if it exists.
    #[must_use]
    pub fn channel(&self, index: usize) -> Option<&Image<f32>> {
        self.channels.get(index)
    }

    /// A new volume holding only the listed channels, in the given order.
    ///
    /// # Errors
    /// Returns [`RenderError::Empty`] for an empty selection and
    /// [`RenderError::ChannelIndex`] for an index that does not exist.
    pub fn select_channels(&self, indices: &[usize]) -> Result<Self, RenderError> {
        if indices.is_empty() {
            return Err(RenderError::Empty);
        }
        let mut picked = Vec::with_capacity(indices.len());
        for &i in indices {
            picked.push(self.channel(i).ok_or(RenderError::ChannelIndex(i))?.clone());
        }
        Self::new(picked, self.spacing)
    }

    /// Minimum and maximum finite sample of a channel, or `None` if the channel is missing
    /// or holds no finite samples.
    #[must_use]
    pub fn intensity_range(&self, index: usize) -> Option<(f32, f32)> {
        let ch = self.channel(index)?;
        ch.data().iter().filter(|v| v.is_finite()).fold(None, |acc, &v| match acc {
            None => Some((v, v)),
            Some((lo, hi)) => Some((lo.min(v), hi.max(v))),
        })
    }

    pub(super) fn channels(&self) -> &[Image<f32>] {
        &self.channels
    }

    /// Replaces the channels, keeping the spacing. Used by level-of-detail code.
    pub(super) fn with_channels(
        channels: Vec<Image<f32>>,
        spacing: [f64; 3],
    ) -> Result<Self, RenderError> {
        Self::new(channels, spacing)
    }
}

/// Trilinearly interpolates `img` at world position `p`.
///
/// Positions outside the volume are clamped to its outermost voxel centres.
pub(super) fn sample_trilinear(img: &Image<f32>, spacing: [f64; 3], p: [f64; 3]) -> f32 {
    let dims = [img.width(), img.height(), img.depth()];
    let mut lo = [0usize; 3];
    let mut hi = [0usize; 3];
    let mut frac = [0f64; 3];
    for a in 0..3 {
        let last = (dims[a] - 1) as f64;
        let u = (p[a] / spacing[a] - 0.5).clamp(0.0, last);
        let i = u.floor();
        lo[a] = i as usize;
        hi[a] = (lo[a] + 1).min(dims[a] - 1);
        frac[a] = u - i;
    }
    let data = img.data();
    let mut out = 0.0f64;
    for corner in 0..8usize {
        let pick = |bit: usize, a: usize| if (corner >> bit) & 1 == 0 { lo[a] } else { hi[a] };
        let weight = |bit: usize, a: usize| {
            if (corner >> bit) & 1 == 0 { 1.0 - frac[a] } else { frac[a] }
        };
        let (x, y, z) = (pick(0, 0), pick(1, 1), pick(2, 2));
        let w = weight(0, 0) * weight(1, 1) * weight(2, 2);
        out += w * f64::from(data[img.index(x, y, z)]);
    }
    out as f32
}

/// Linear index of the voxel containing world position `p`, clamped to the volume.
pub(super) fn nearest_voxel(shape: (usize, usize, usize), spacing: [f64; 3], p: [f64; 3]) -> usize {
    let dims = [shape.0, shape.1, shape.2];
    let mut idx = [0usize; 3];
    for a in 0..3 {
        let u = (p[a] / spacing[a]).floor().clamp(0.0, (dims[a] - 1) as f64);
        idx[a] = u as usize;
    }
    (idx[2] * shape.1 + idx[1]) * shape.0 + idx[0]
}

#[allow(clippy::float_cmp)]
#[cfg(test)]
mod tests {
    use super::*;

    fn ramp(w: usize, h: usize, d: usize) -> Image<f32> {
        let data = (0..w * h * d).map(|i| i as f32).collect();
        Image::new_3d(w, h, d, data).unwrap()
    }

    #[test]
    fn render3d_volume_rejects_bad_input() {
        assert_eq!(Volume::new(vec![], [1.0; 3]), Err(RenderError::Empty));
        assert!(Volume::new(vec![ramp(2, 2, 2)], [0.0, 1.0, 1.0]).is_err());
        let a = ramp(2, 2, 2);
        let b = ramp(3, 2, 2);
        assert!(Volume::new(vec![a, b], [1.0; 3]).is_err());
    }

    #[test]
    fn render3d_volume_extent_and_channels() {
        let v = Volume::new(vec![ramp(4, 3, 2), ramp(4, 3, 2)], [0.5, 1.0, 2.0]).unwrap();
        assert_eq!(v.shape(), (4, 3, 2));
        assert_eq!(v.extent(), [2.0, 3.0, 4.0]);
        let one = v.select_channels(&[1]).unwrap();
        assert_eq!(one.channel_count(), 1);
        assert_eq!(v.select_channels(&[5]), Err(RenderError::ChannelIndex(5)));
        assert_eq!(v.intensity_range(0), Some((0.0, 23.0)));
    }

    #[test]
    fn render3d_volume_from_u16_is_exact() {
        let img = Image::new_2d(2, 1, vec![0u16, 65535]).unwrap();
        let v = Volume::from_u16(&[img], [1.0; 3]).unwrap();
        assert_eq!(v.channel(0).unwrap().data(), &[0.0, 65535.0]);
    }

    #[test]
    fn render3d_trilinear_hits_voxel_centres_and_midpoints() {
        let img = ramp(2, 1, 1); // values 0 and 1 along x
        let s = [1.0; 3];
        assert_eq!(sample_trilinear(&img, s, [0.5, 0.5, 0.5]), 0.0);
        assert_eq!(sample_trilinear(&img, s, [1.5, 0.5, 0.5]), 1.0);
        assert!((sample_trilinear(&img, s, [1.0, 0.5, 0.5]) - 0.5).abs() < 1e-6);
        // Outside the volume the nearest edge voxel is used.
        assert_eq!(sample_trilinear(&img, s, [-10.0, 0.5, 0.5]), 0.0);
        assert_eq!(sample_trilinear(&img, s, [10.0, 0.5, 0.5]), 1.0);
    }
}
