//! A seam for choosing how a volume is rendered.
//!
//! [`RenderBackend`] is the synchronous entry point that a native renderer implements.
//! [`CpuBackend`] is the only implementation today. It is the reference that any other
//! backend is checked against, and it calls the free functions in [`super::raycast`]
//! unchanged.
//!
//! A native GPU backend is planned (see `docs/dependency-audit.md`, "Native GPU
//! candidates"). It would implement this trait without changing callers. The browser
//! WebGPU path in the `webgpu` module is asynchronous and wasm-only, so it does not
//! implement this trait.

use super::camera::Camera;
use super::raycast::{RenderSettings, render_volume};
use super::transfer::TransferFunction;
use super::volume::Volume;
use super::{RenderError, Rgba8};
use crate::image::Image;

/// Renders a [`Volume`] to an RGBA image of the camera's viewport.
pub trait RenderBackend {
    /// A short, stable name for this backend, for logs and diagnostics.
    fn name(&self) -> &'static str;

    /// Renders the intensity volume through `tf`.
    ///
    /// # Errors
    /// Same as [`render_volume`].
    fn render_volume(
        &self,
        volume: &Volume,
        camera: &Camera,
        tf: &TransferFunction,
        settings: &RenderSettings,
    ) -> Result<Image<Rgba8>, RenderError>;
}

/// The CPU reference renderer. Rows run in parallel on native targets.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct CpuBackend;

impl RenderBackend for CpuBackend {
    fn name(&self) -> &'static str {
        "cpu"
    }

    fn render_volume(
        &self,
        volume: &Volume,
        camera: &Camera,
        tf: &TransferFunction,
        settings: &RenderSettings,
    ) -> Result<Image<Rgba8>, RenderError> {
        render_volume(volume, camera, tf, settings)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::render3d::camera::Projection;
    use crate::render3d::transfer::ChannelTransfer;

    fn scene() -> (Volume, Camera, TransferFunction, RenderSettings) {
        let img = Image::new_3d(8, 8, 8, vec![1.0f32; 512]).unwrap();
        let vol = Volume::new(vec![img], [1.0; 3]).unwrap();
        let mut ch = ChannelTransfer::new([1.0, 0.5, 0.0], (0.0, 1.0));
        ch.opacity = vec![(0.0, 1.0), (1.0, 1.0)];
        let tf = TransferFunction { channels: vec![ch], clip: None };
        let cam = Camera {
            target: [4.0; 3],
            yaw: 0.3,
            pitch: 0.2,
            distance: 100.0,
            projection: Projection::Orthographic { height: 16.0 },
            viewport: (16, 16),
        };
        (vol, cam, tf, RenderSettings::new(0.5))
    }

    #[test]
    fn render3d_backend_cpu_matches_free_function() {
        let (vol, cam, tf, settings) = scene();
        let direct = render_volume(&vol, &cam, &tf, &settings).unwrap();
        let via_trait = CpuBackend.render_volume(&vol, &cam, &tf, &settings).unwrap();
        assert_eq!(direct, via_trait);
        assert_eq!(CpuBackend.name(), "cpu");
    }

    #[test]
    fn render3d_backend_is_usable_as_a_trait_object() {
        let (vol, cam, tf, settings) = scene();
        let backend: &dyn RenderBackend = &CpuBackend;
        assert!(backend.render_volume(&vol, &cam, &tf, &settings).is_ok());
    }
}
