//! Streamed level of detail for volumes too large to hold at full resolution.
//!
//! A [`SlabSource`] yields a volume one z-slab at a time, so a reader backed by a file, a
//! memory map or a remote store never has to make the whole volume resident. The coarse
//! levels produced here are block-max reductions, identical to [`lod::downsample_max`] on
//! the whole volume. Only one slab of full-resolution data is held at a time, plus the
//! coarse output.
//!
//! Slabs start at multiples of the slab depth, and the slab depth must be a multiple of the
//! downsample factor. That keeps every block inside one slab, so the streamed result matches
//! the in-memory result exactly.
//!
//! [`lod::downsample_max`]: super::lod::downsample_max

use std::marker::PhantomData;

use super::RenderError;
use super::camera::Camera;
use super::lod::{choose_level, downsample_max};
use super::volume::Volume;
use crate::image::Image;
use crate::io::ome::OmeTiff;
use crate::io::sample::Sample;
use crate::source::ByteSource;

/// A volume that can be read one z-slab at a time.
pub trait SlabSource {
    /// Shape `(width, height, depth)` shared by every channel.
    fn shape(&self) -> (usize, usize, usize);

    /// Number of channels.
    fn channel_count(&self) -> usize;

    /// Physical voxel spacing `[sx, sy, sz]`.
    fn spacing(&self) -> [f64; 3];

    /// Reads slices `z0..z1` of `channel`. The returned image has depth `z1 - z0`.
    ///
    /// # Errors
    /// Implementations return a [`RenderError`] when the slab cannot be read.
    fn read_slab(&mut self, channel: usize, z0: usize, z1: usize) -> Result<Image<f32>, RenderError>;
}

impl SlabSource for Volume {
    fn shape(&self) -> (usize, usize, usize) {
        Volume::shape(self)
    }

    fn channel_count(&self) -> usize {
        Volume::channel_count(self)
    }

    fn spacing(&self) -> [f64; 3] {
        Volume::spacing(self)
    }

    fn read_slab(&mut self, channel: usize, z0: usize, z1: usize) -> Result<Image<f32>, RenderError> {
        let ch = self.channel(channel).ok_or(RenderError::ChannelIndex(channel))?;
        let (w, h, d) = ch.shape();
        if z0 > z1 || z1 > d {
            return Err(RenderError::InvalidParameter("slab range is outside the volume"));
        }
        let plane = w * h;
        Ok(Image::new_3d(w, h, z1 - z0, ch.data()[z0 * plane..z1 * plane].to_vec())?)
    }
}

/// Largest slab depth, a multiple of `factor`, whose slab of `f32` samples fits in
/// `budget_bytes`. Returns `factor` when even one block does not fit, so a streamed pass
/// always makes progress.
///
/// # Errors
/// Returns [`RenderError::InvalidParameter`] when `factor` is zero.
pub fn slab_depth_for_budget(
    shape: (usize, usize, usize),
    factor: usize,
    budget_bytes: usize,
) -> Result<usize, RenderError> {
    if factor == 0 {
        return Err(RenderError::InvalidParameter("downsample factor must be at least 1"));
    }
    let plane_bytes = shape.0.max(1) * shape.1.max(1) * std::mem::size_of::<f32>();
    let slices = budget_bytes / plane_bytes;
    Ok((slices / factor).max(1) * factor)
}

/// Block-max downsamples every channel of `source` by `factor`, reading `slab_depth` slices
/// at a time.
///
/// The result equals `Volume::downsample_max(factor)` on the fully loaded volume.
///
/// # Errors
/// Returns [`RenderError::InvalidParameter`] when `factor` is zero, when `slab_depth` is zero
/// or not a multiple of `factor`, or when the source returns a slab of the wrong shape. Source
/// read errors and [`RenderError::Empty`] for an empty volume are also returned.
pub fn downsample_streamed<S: SlabSource + ?Sized>(
    source: &mut S,
    factor: usize,
    slab_depth: usize,
) -> Result<Volume, RenderError> {
    if factor == 0 {
        return Err(RenderError::InvalidParameter("downsample factor must be at least 1"));
    }
    if slab_depth == 0 || !slab_depth.is_multiple_of(factor) {
        return Err(RenderError::InvalidParameter(
            "slab depth must be a positive multiple of the downsample factor",
        ));
    }
    let (w, h, d) = source.shape();
    if w == 0 || h == 0 || d == 0 || source.channel_count() == 0 {
        return Err(RenderError::Empty);
    }
    let (nw, nh, nd) = (w.div_ceil(factor), h.div_ceil(factor), d.div_ceil(factor));
    let mut channels = Vec::with_capacity(source.channel_count());
    for c in 0..source.channel_count() {
        let mut out = Vec::with_capacity(nw * nh * nd);
        let mut z0 = 0;
        while z0 < d {
            let z1 = (z0 + slab_depth).min(d);
            let slab = source.read_slab(c, z0, z1)?;
            if slab.shape() != (w, h, z1 - z0) {
                return Err(RenderError::InvalidParameter("source returned a slab of the wrong shape"));
            }
            out.extend_from_slice(downsample_max(&slab, factor)?.data());
            z0 = z1;
        }
        channels.push(Image::new_3d(nw, nh, nd, out)?);
    }
    let s = source.spacing();
    let f = factor as f64;
    Volume::new(channels, [s[0] * f, s[1] * f, s[2] * f])
}

/// Streams one coarse [`Volume`] per factor in `factors`, in the same order. Each level
/// reads the source again, so no full-resolution data is kept between levels.
///
/// # Errors
/// Returns the first error from [`downsample_streamed`].
pub fn streamed_levels<S: SlabSource + ?Sized>(
    source: &mut S,
    factors: &[usize],
    slab_depth_budget_bytes: usize,
) -> Result<Vec<Volume>, RenderError> {
    factors
        .iter()
        .map(|&f| {
            let depth = slab_depth_for_budget(source.shape(), f, slab_depth_budget_bytes)?;
            downsample_streamed(source, f, depth)
        })
        .collect()
}

/// One series of an OME-TIFF file, read as a z-stack one plane at a time.
///
/// Channel `c` of the volume is the file's channel `c` at time point `t`. Voxel spacing is
/// not stored in the plane data, so the caller gives it.
pub struct OmeSlabs<'a, S: ByteSource, T> {
    tiff: &'a OmeTiff<S>,
    series: usize,
    t: usize,
    shape: (usize, usize, usize),
    channels: usize,
    spacing: [f64; 3],
    _sample: PhantomData<T>,
}

impl<S: ByteSource, T> std::fmt::Debug for OmeSlabs<'_, S, T> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("OmeSlabs")
            .field("series", &self.series)
            .field("t", &self.t)
            .field("shape", &self.shape)
            .field("channels", &self.channels)
            .finish_non_exhaustive()
    }
}

impl<'a, S: ByteSource, T: Sample + Into<f32>> OmeSlabs<'a, S, T> {
    /// Wraps series `series` at time point `t`.
    ///
    /// # Errors
    /// Returns [`RenderError::InvalidParameter`] if the series does not exist, its pixel type
    /// is not `T`, or the time point is out of range, or if the spacing is not positive.
    pub fn new(
        tiff: &'a OmeTiff<S>,
        series: usize,
        t: usize,
        spacing: [f64; 3],
    ) -> Result<Self, RenderError> {
        let image = tiff
            .images()
            .get(series)
            .ok_or(RenderError::InvalidParameter("OME-TIFF series out of range"))?;
        let p = &image.pixels;
        if p.pixel_type != T::PIXEL_TYPE {
            return Err(RenderError::InvalidParameter("OME-TIFF pixel type does not match T"));
        }
        if t >= p.size_t {
            return Err(RenderError::InvalidParameter("OME-TIFF time point out of range"));
        }
        if spacing.iter().any(|s| !s.is_finite() || *s <= 0.0) {
            return Err(RenderError::InvalidParameter("voxel spacing must be finite and positive"));
        }
        Ok(Self {
            tiff,
            series,
            t,
            shape: (p.size_x, p.size_y, p.size_z),
            channels: p.size_c,
            spacing,
            _sample: PhantomData,
        })
    }
}

impl<S: ByteSource, T: Sample + Into<f32>> SlabSource for OmeSlabs<'_, S, T> {
    fn shape(&self) -> (usize, usize, usize) {
        self.shape
    }

    fn channel_count(&self) -> usize {
        self.channels
    }

    fn spacing(&self) -> [f64; 3] {
        self.spacing
    }

    fn read_slab(&mut self, channel: usize, z0: usize, z1: usize) -> Result<Image<f32>, RenderError> {
        if channel >= self.channels {
            return Err(RenderError::ChannelIndex(channel));
        }
        let (w, h, d) = self.shape;
        if z0 > z1 || z1 > d {
            return Err(RenderError::InvalidParameter("slab range is outside the volume"));
        }
        let mut data = Vec::with_capacity(w * h * (z1 - z0));
        for z in z0..z1 {
            let plane = self
                .tiff
                .read_plane::<T>(self.series, z, channel, self.t)
                .map_err(|e| RenderError::Read(e.to_string()))?;
            data.extend(plane.data().iter().map(|&v| Into::<f32>::into(v)));
        }
        Ok(Image::new_3d(w, h, z1 - z0, data)?)
    }
}

/// Streams only the level of detail that suits the camera, and returns it with its factor.
///
/// The factor is chosen by [`choose_level`]: the coarsest candidate whose voxels cover at most
/// `max_pixels` on screen. The largest axis spacing is used, so the estimate errs toward finer
/// levels. Only that level is built, in one pass over the source. The full-resolution volume
/// is never kept, unless it is itself the chosen level.
///
/// # Errors
/// Returns [`RenderError::Empty`] when `factors` is empty, and otherwise the errors of
/// [`downsample_streamed`] for the chosen level.
pub fn stream_level_for_camera<S: SlabSource + ?Sized>(
    source: &mut S,
    camera: &Camera,
    factors: &[usize],
    max_pixels: f64,
    budget_bytes: usize,
) -> Result<(usize, Volume), RenderError> {
    if factors.is_empty() {
        return Err(RenderError::Empty);
    }
    let spacing = source.spacing().into_iter().fold(f64::MIN, f64::max);
    let factor = factors[choose_level(camera, spacing, factors, max_pixels)];
    let depth = slab_depth_for_budget(source.shape(), factor, budget_bytes)?;
    Ok((factor, downsample_streamed(source, factor, depth)?))
}

#[allow(clippy::float_cmp)]
#[cfg(test)]
mod tests {
    use super::*;

    fn pattern(w: usize, h: usize, d: usize, seed: u32) -> Image<f32> {
        let data = (0..w * h * d)
            .map(|i| {
                let v = (i as u32).wrapping_mul(2_654_435_761).wrapping_add(seed);
                (v % 997) as f32
            })
            .collect();
        Image::new_3d(w, h, d, data).unwrap()
    }

    /// Wraps a volume and records every slab read, to check memory bounds.
    struct Counting {
        inner: Volume,
        max_slab_depth: usize,
        slices_read: Vec<usize>,
    }

    impl SlabSource for Counting {
        fn shape(&self) -> (usize, usize, usize) {
            self.inner.shape()
        }
        fn channel_count(&self) -> usize {
            self.inner.channel_count()
        }
        fn spacing(&self) -> [f64; 3] {
            self.inner.spacing()
        }
        fn read_slab(&mut self, c: usize, z0: usize, z1: usize) -> Result<Image<f32>, RenderError> {
            self.max_slab_depth = self.max_slab_depth.max(z1 - z0);
            if self.slices_read.len() <= c {
                self.slices_read.resize(c + 1, 0);
            }
            self.slices_read[c] += z1 - z0;
            self.inner.read_slab(c, z0, z1)
        }
    }

    fn two_channel(w: usize, h: usize, d: usize) -> Volume {
        Volume::new(vec![pattern(w, h, d, 1), pattern(w, h, d, 7)], [0.5, 1.0, 2.0]).unwrap()
    }

    #[test]
    fn render3d_stream_matches_in_memory_downsample() {
        // Odd sizes exercise partial blocks at the high edges in every axis.
        let vol = two_channel(7, 5, 11);
        for factor in [1usize, 2, 3, 4] {
            let reference = vol.downsample_max(factor).unwrap();
            let mut src = vol.clone();
            // Every slab depth that is a multiple of the factor, up to the whole volume.
            for k in 1..=(11usize.div_ceil(factor)) {
                let got = downsample_streamed(&mut src, factor, k * factor).unwrap();
                assert_eq!(got, reference, "factor {factor}, slab depth {}", k * factor);
            }
        }
    }

    #[test]
    fn render3d_stream_reads_bounded_slabs_and_covers_volume() {
        let vol = two_channel(4, 4, 10);
        let mut src = Counting { inner: vol, max_slab_depth: 0, slices_read: vec![] };
        downsample_streamed(&mut src, 2, 4).unwrap();
        assert!(src.max_slab_depth <= 4, "slab of {} slices", src.max_slab_depth);
        // Each channel is read exactly once through, in full.
        assert_eq!(src.slices_read, vec![10, 10]);
    }

    #[test]
    fn render3d_stream_rejects_bad_parameters() {
        let mut vol = two_channel(4, 4, 4);
        assert!(downsample_streamed(&mut vol, 0, 4).is_err());
        assert!(downsample_streamed(&mut vol, 2, 0).is_err());
        assert!(downsample_streamed(&mut vol, 2, 3).is_err(), "depth not a multiple of factor");
        assert!(slab_depth_for_budget((4, 4, 4), 0, 1024).is_err());
    }

    #[test]
    fn render3d_stream_budget_picks_block_aligned_depth() {
        // A 10x10 plane is 400 bytes, so a 2000-byte budget holds 5 slices.
        assert_eq!(slab_depth_for_budget((10, 10, 50), 1, 2000).unwrap(), 5);
        assert_eq!(slab_depth_for_budget((10, 10, 50), 2, 2000).unwrap(), 4);
        // A budget smaller than one block still gives one full block of progress.
        assert_eq!(slab_depth_for_budget((10, 10, 50), 4, 10).unwrap(), 4);
    }

    #[test]
    fn render3d_stream_detects_wrong_slab_shape() {
        struct Liar(Volume);
        impl SlabSource for Liar {
            fn shape(&self) -> (usize, usize, usize) {
                self.0.shape()
            }
            fn channel_count(&self) -> usize {
                1
            }
            fn spacing(&self) -> [f64; 3] {
                self.0.spacing()
            }
            fn read_slab(&mut self, _: usize, _: usize, _: usize) -> Result<Image<f32>, RenderError> {
                Ok(pattern(1, 1, 1, 0))
            }
        }
        let mut liar = Liar(two_channel(2, 2, 2));
        assert!(downsample_streamed(&mut liar, 1, 2).is_err());
    }

    #[test]
    fn render3d_stream_levels_match_per_factor_reference() {
        let vol = two_channel(9, 6, 8);
        let mut src = vol.clone();
        let levels = streamed_levels(&mut src, &[1, 2, 4], 600).unwrap();
        assert_eq!(levels.len(), 3);
        for (level, factor) in levels.iter().zip([1usize, 2, 4]) {
            assert_eq!(*level, vol.downsample_max(factor).unwrap(), "factor {factor}");
        }
    }

    fn camera_at(distance: f64) -> Camera {
        use super::super::camera::Projection;
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
    fn render3d_stream_level_matches_chosen_in_memory_level() {
        let vol = two_channel(12, 10, 9);
        let factors = [1usize, 2, 4];
        for distance in [5.0, 200.0, 5000.0] {
            let cam = camera_at(distance);
            let mut src = vol.clone();
            let (factor, level) =
                stream_level_for_camera(&mut src, &cam, &factors, 2.0, 4096).unwrap();
            let idx = choose_level(&cam, 2.0, &factors, 2.0);
            assert_eq!(factor, factors[idx]);
            assert_eq!(level, vol.downsample_max(factor).unwrap(), "distance {distance}");
        }
    }

    #[test]
    fn render3d_stream_level_is_coarser_when_far() {
        let vol = two_channel(16, 16, 16);
        let factors = [1usize, 2, 4, 8];
        let mut near_src = vol.clone();
        let mut far_src = vol.clone();
        let (near, _) = stream_level_for_camera(&mut near_src, &camera_at(3.0), &factors, 1.5, 8192).unwrap();
        let (far, _) = stream_level_for_camera(&mut far_src, &camera_at(3000.0), &factors, 1.5, 8192).unwrap();
        assert!(far > near, "far {far} should exceed near {near}");
    }

    #[test]
    fn render3d_stream_level_reads_source_once_and_rejects_empty_factors() {
        let vol = two_channel(6, 6, 10);
        let mut src = Counting { inner: vol, max_slab_depth: 0, slices_read: vec![] };
        stream_level_for_camera(&mut src, &camera_at(3000.0), &[2, 4], 1.0, 1 << 20).unwrap();
        assert_eq!(src.slices_read, vec![10, 10]);
        let mut empty = two_channel(4, 4, 4);
        assert_eq!(
            stream_level_for_camera(&mut empty, &camera_at(3.0), &[], 1.0, 1024),
            Err(RenderError::Empty)
        );
    }

    #[test]
    fn render3d_stream_ome_tiff_slabs_match_in_memory_volume() {
        use crate::io::ome::{OmeTiff, encode_ome_tiff};
        // Single-channel, single time point: plane index equals z in the file.
        let (w, h, d) = (7usize, 5usize, 6usize);
        let data: Vec<u16> = (0..w * h * d).map(|i| ((i * 37) % 997) as u16).collect();
        let img = Image::new_3d(w, h, d, data).unwrap();
        let bytes = encode_ome_tiff(&img, (d, 1, 1), Some("slabs"), false).unwrap();
        let tiff = OmeTiff::open(bytes).unwrap();
        let spacing = [0.5, 0.5, 2.0];
        let reference = Volume::from_u16(std::slice::from_ref(&img), spacing).unwrap();
        let mut slabs = OmeSlabs::<_, u16>::new(&tiff, 0, 0, spacing).unwrap();
        assert_eq!(slabs.shape(), (w, h, d));
        assert_eq!(slabs.channel_count(), 1);
        for depth in [2usize, 4] {
            let mut slabs = OmeSlabs::<_, u16>::new(&tiff, 0, 0, spacing).unwrap();
            let got = downsample_streamed(&mut slabs, 2, depth).unwrap();
            assert_eq!(got, reference.downsample_max(2).unwrap(), "slab depth {depth}");
        }
        // Reading a slab directly returns the planes in z order.
        let slab = slabs.read_slab(0, 1, 3).unwrap();
        assert_eq!(slab.shape(), (w, h, 2));
        assert_eq!(slab.data()[0], reference.channel(0).unwrap().data()[w * h]);
    }

    #[test]
    fn render3d_stream_ome_tiff_rejects_wrong_type_and_range() {
        use crate::io::ome::{OmeTiff, encode_ome_tiff};
        let img = Image::new_3d(3, 3, 2, vec![1u16; 18]).unwrap();
        let bytes = encode_ome_tiff(&img, (2, 1, 1), None, false).unwrap();
        let tiff = OmeTiff::open(bytes).unwrap();
        assert!(OmeSlabs::<_, f32>::new(&tiff, 0, 0, [1.0; 3]).is_err(), "pixel type");
        assert!(OmeSlabs::<_, u16>::new(&tiff, 1, 0, [1.0; 3]).is_err(), "series");
        assert!(OmeSlabs::<_, u16>::new(&tiff, 0, 3, [1.0; 3]).is_err(), "time point");
        assert!(OmeSlabs::<_, u16>::new(&tiff, 0, 0, [0.0, 1.0, 1.0]).is_err(), "spacing");
        let mut slabs = OmeSlabs::<_, u16>::new(&tiff, 0, 0, [1.0; 3]).unwrap();
        assert!(matches!(slabs.read_slab(1, 0, 1), Err(RenderError::ChannelIndex(1))));
        assert!(slabs.read_slab(0, 0, 5).is_err());
    }

    #[test]
    fn render3d_stream_ome_tiff_multichannel_planes_in_order() {
        use crate::io::ome::{OmeTiff, encode_ome_tiff};
        // Two channels of three planes. The file stores them channel-major, so channel c,
        // plane z is volume plane c * size_z + z.
        let (w, h, sz, sc) = (4usize, 3usize, 3usize, 2usize);
        let data: Vec<u16> = (0..w * h * sz * sc).map(|i| (i as u16) * 3 + 1).collect();
        let img = Image::new_3d(w, h, sz * sc, data).unwrap();
        let bytes = encode_ome_tiff(&img, (sz, sc, 1), None, false).unwrap();
        let tiff = OmeTiff::open(bytes).unwrap();
        let reference = Volume::from_u16(std::slice::from_ref(&img), [1.0; 3]).unwrap();
        let mut slabs = OmeSlabs::<_, u16>::new(&tiff, 0, 0, [1.0; 3]).unwrap();
        assert_eq!(slabs.channel_count(), sc);
        for c in 0..sc {
            let slab = slabs.read_slab(c, 0, sz).unwrap();
            let want = &reference.channel(0).unwrap().data()[c * sz * w * h..(c + 1) * sz * w * h];
            assert_eq!(slab.data(), want, "channel {c}");
        }
    }
}
