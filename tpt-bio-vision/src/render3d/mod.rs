//! CPU reference renderer for 3D microscopy volumes.
//!
//! This is the foundation for the Phase 4 3D viewer. Everything here runs on the CPU and
//! renders to plain RGBA buffers. No GPU, windowing or interactive code exists yet.
//!
//! - [`backend`]: the [`backend::RenderBackend`] trait and the CPU implementation.
//! - [`volume`]: multi-channel volumes with physical voxel spacing and trilinear sampling.
//! - [`camera`]: orbit camera with perspective or orthographic projection.
//! - [`transfer`]: per-channel windowed transfer functions and a clipping plane.
//! - [`raycast`]: front-to-back alpha compositing along rays through a volume.
//! - [`overlay`]: label overlay with first-hit surface blending and a deterministic palette.
//! - [`isosurface`]: marching-tetrahedra isosurfaces producing a triangle [`isosurface::Mesh`].
//! - [`lod`]: block-max downsampling and screen-footprint level selection.
//! - `webgpu` (feature `webgpu`, wasm32): the same marching through the browser's WebGPU.
//! - [`stream`]: block-max levels streamed from a [`stream::SlabSource`], one z-slab at a time.
//! - [`rasterise`]: flat-shaded z-buffer rasteriser for isosurface meshes.

pub mod backend;
pub mod camera;
pub mod isosurface;
pub mod lod;
pub mod overlay;
pub mod rasterise;
pub mod raycast;
pub mod stream;
pub mod transfer;
pub mod volume;
#[cfg(all(feature = "webgpu", target_arch = "wasm32"))]
pub mod webgpu;

use std::fmt;

use crate::image::ImageError;

/// An 8-bit RGBA pixel.
pub type Rgba8 = [u8; 4];

/// Errors from building or rendering a 3D scene.
#[derive(Debug, Clone, PartialEq)]
pub enum RenderError {
    /// An image operation failed, for example because shapes differ.
    Image(ImageError),
    /// A volume had no channels or no samples.
    Empty,
    /// A numeric parameter was out of its valid range.
    InvalidParameter(&'static str),
    /// A transfer function does not have one entry per volume channel.
    ChannelMismatch {
        /// Number of channels in the volume.
        expected: usize,
        /// Number of channel entries supplied.
        actual: usize,
    },
    /// A requested channel index does not exist.
    ChannelIndex(usize),
    /// A slab could not be read from its source. The message says why.
    Read(String),
}

impl fmt::Display for RenderError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Image(e) => write!(f, "{e}"),
            Self::Empty => f.write_str("volume has no channels or no samples"),
            Self::InvalidParameter(msg) => f.write_str(msg),
            Self::ChannelMismatch { expected, actual } => {
                write!(f, "transfer function has {actual} channels, volume has {expected}")
            }
            Self::ChannelIndex(i) => write!(f, "channel index {i} is out of range"),
            Self::Read(m) => write!(f, "could not read slab: {m}"),
        }
    }
}

impl std::error::Error for RenderError {}

impl From<ImageError> for RenderError {
    fn from(e: ImageError) -> Self {
        Self::Image(e)
    }
}

type Vec3 = [f64; 3];

fn add(a: Vec3, b: Vec3) -> Vec3 {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}

fn sub(a: Vec3, b: Vec3) -> Vec3 {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

fn scale(a: Vec3, s: f64) -> Vec3 {
    [a[0] * s, a[1] * s, a[2] * s]
}

fn dot(a: Vec3, b: Vec3) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

fn cross(a: Vec3, b: Vec3) -> Vec3 {
    [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]]
}

fn length(a: Vec3) -> f64 {
    dot(a, a).sqrt()
}

/// Returns `a / |a|`, or `a` unchanged when it has zero length.
fn normalise(a: Vec3) -> Vec3 {
    let l = length(a);
    if l > 0.0 { scale(a, 1.0 / l) } else { a }
}

/// Converts a colour in `[0, 1]` to an 8-bit channel, clamping out-of-range values.
fn to_u8(v: f32) -> u8 {
    (v.clamp(0.0, 1.0) * 255.0).round() as u8
}
