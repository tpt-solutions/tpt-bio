//! Level of detail: block-max downsampling and screen-footprint level selection.
//!
//! [`downsample_max`] reduces a volume by an integer factor per axis, keeping the maximum
//! of each block. Blocks at the high edges may be partial. Their maximum uses only the voxels
//! that exist, and the coarse spacing is still `factor * spacing`, so the coarse extent can
//! exceed the original by less than one coarse voxel.
//!
//! [`choose_level`] compares candidate factors against a camera. A level's footprint is the
//! number of pixels one of its voxels covers at the camera's target depth. It returns the
//! coarsest level whose footprint is at most `max_pixels`, so voxels never draw larger
//! than the limit. Distant volumes use coarse data and close volumes use fine data.

use super::RenderError;
use super::camera::Camera;
use super::volume::Volume;
use crate::image::Image;

/// Block-max downsampling of one channel by an integer `factor` on each axis.
///
/// The output has shape `ceil(n / factor)` per axis.
///
/// # Errors
/// Returns [`RenderError::InvalidParameter`] when `factor` is zero.
pub fn downsample_max(img: &Image<f32>, factor: usize) -> Result<Image<f32>, RenderError> {
    if factor == 0 {
        return Err(RenderError::InvalidParameter("downsample factor must be at least 1"));
    }
    let (w, h, d) = img.shape();
    let (nw, nh, nd) = (w.div_ceil(factor), h.div_ceil(factor), d.div_ceil(factor));
    let mut out = vec![f32::NEG_INFINITY; nw * nh * nd];
    let data = img.data();
    for z in 0..d {
        for y in 0..h {
            for x in 0..w {
                let v = data[img.index(x, y, z)];
                let o = ((z / factor) * nh + y / factor) * nw + x / factor;
                if v > out[o] {
                    out[o] = v;
                }
            }
        }
    }
    Ok(Image::new_3d(nw, nh, nd, out)?)
}

impl Volume {
    /// Block-max downsampling of every channel. The spacing becomes `factor * spacing`.
    ///
    /// # Errors
    /// Returns [`RenderError::InvalidParameter`] when `factor` is zero.
    pub fn downsample_max(&self, factor: usize) -> Result<Volume, RenderError> {
        let channels = self
            .channels()
            .iter()
            .map(|ch| downsample_max(ch, factor))
            .collect::<Result<Vec<_>, _>>()?;
        let s = self.spacing();
        let f = factor as f64;
        Volume::with_channels(channels, [s[0] * f, s[1] * f, s[2] * f])
    }
}

/// Pixels covered by one voxel of size `spacing * factor` at the camera's target depth.
#[must_use]
pub fn footprint_pixels(camera: &Camera, spacing: f64, factor: usize) -> f64 {
    camera.pixels_per_unit(camera.distance) * spacing * factor as f64
}

/// Index into `factors` of the coarsest level whose voxels still cover at most `max_pixels`.
///
/// This is the usual level-of-detail rule: a voxel drawn larger than `max_pixels` on screen
/// is too coarse to look right, so the coarsest level that stays within the limit is used.
/// Distant volumes have small footprints, so they get the coarse levels; close volumes
/// get the fine levels. Footprints come from [`footprint_pixels`]. If no level meets the
/// limit, the finest level (smallest factor) is returned. Returns `0` for an empty list, so
/// callers should not pass one.
#[must_use]
pub fn choose_level(camera: &Camera, spacing: f64, factors: &[usize], max_pixels: f64) -> usize {
    let mut best: Option<(usize, usize)> = None; // (index, factor) of the best fit so far
    let mut finest = (0usize, usize::MAX);
    for (i, &f) in factors.iter().enumerate() {
        if f < finest.1 {
            finest = (i, f);
        }
        if footprint_pixels(camera, spacing, f) <= max_pixels && best.is_none_or(|(_, bf)| f > bf) {
            best = Some((i, f));
        }
    }
    best.map_or(finest.0, |(i, _)| i)
}

#[allow(clippy::float_cmp)]
#[cfg(test)]
mod tests {
    use super::*;
    use crate::render3d::camera::Projection;

    #[test]
    fn render3d_lod_constant_volume_stays_constant() {
        let img = Image::new_3d(7, 5, 3, vec![0.25f32; 105]).unwrap();
        for factor in [1, 2, 3, 4] {
            let down = downsample_max(&img, factor).unwrap();
            assert!(down.data().iter().all(|&v| v == 0.25), "factor {factor}");
            assert_eq!(
                down.shape(),
                (7usize.div_ceil(factor), 5usize.div_ceil(factor), 3usize.div_ceil(factor))
            );
        }
        let vol = Volume::new(vec![img], [1.0, 2.0, 3.0]).unwrap();
        let coarse = vol.downsample_max(2).unwrap();
        assert_eq!(coarse.spacing(), [2.0, 4.0, 6.0]);
        assert!(coarse.channel(0).unwrap().data().iter().all(|&v| v == 0.25));
    }

    #[test]
    fn render3d_lod_takes_block_maximum() {
        let img = Image::new_2d(4, 2, vec![1.0f32, 5.0, 0.0, 0.0, 2.0, 3.0, 0.0, 9.0]).unwrap();
        let down = downsample_max(&img, 2).unwrap();
        assert_eq!(down.shape(), (2, 1, 1));
        assert_eq!(down.data(), &[5.0, 9.0]);
        assert!(downsample_max(&img, 0).is_err());
    }

    fn camera_at(distance: f64) -> Camera {
        Camera {
            target: [0.0; 3],
            yaw: 0.3,
            pitch: 0.2,
            distance,
            projection: Projection::Perspective { fov_y: 1.0 },
            viewport: (64, 64),
        }
    }

    #[test]
    fn render3d_lod_level_changes_with_distance() {
        let factors = [1usize, 2, 4, 8, 16];
        // Limit of 4 px per voxel. Near (distance 10) voxels are about 5.9 px at factor 1,
        // so no level fits and the finest (factor 1) is used.
        let near = choose_level(&camera_at(10.0), 1.0, &factors, 4.0);
        // Far (distance 1000) voxels are about 0.06 px at factor 1, so the coarsest
        // level (factor 16) still fits within the limit.
        let far = choose_level(&camera_at(1000.0), 1.0, &factors, 4.0);
        assert_ne!(near, far);
        assert_eq!(factors[near], 1);
        assert_eq!(factors[far], 16);
        // Footprint at distance 10 is 64 / (2 * 10 * tan(0.5)), about 5.9 px per unit.
        let fp = footprint_pixels(&camera_at(10.0), 1.0, 1);
        assert!((fp - 64.0 / (2.0 * 10.0 * 0.5f64.tan())).abs() < 1e-9);
    }
}
