//! Front-to-back ray casting of a [`Volume`] into an RGBA image.
//!
//! Each pixel shoots one ray from the camera. The ray is clipped to the volume's bounding
//! box and marched with a fixed world-space step, taking the midpoint of each step. Samples
//! are trilinear, and the transfer function gives colour and opacity per sample. Samples are
//! composited front to back with the over operator. A ray stops once its accumulated opacity
//! reaches [`OPACITY_CUTOFF`]. Rays that miss the box show the background colour.
//!
//! Opacity is per step, not per unit length: changing `step` changes the apparent density.
//! Output pixels are opaque (alpha 255).

use super::camera::{Camera, Ray};
use super::transfer::TransferFunction;
use super::volume::{Volume, nearest_voxel, sample_trilinear};
use super::{RenderError, Rgba8, add, length, scale, sub, to_u8};
use crate::image::Image;

/// Accumulated opacity at which a ray stops marching.
pub const OPACITY_CUTOFF: f32 = 0.99;

/// Settings for a render.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct RenderSettings {
    /// March step in world units. Must be finite and positive.
    pub step: f64,
    /// Background colour shown behind the volume and for rays that miss it.
    pub background: [u8; 3],
}

impl RenderSettings {
    /// Settings with the given step and a black background.
    #[must_use]
    pub fn new(step: f64) -> Self {
        Self { step, background: [0; 3] }
    }
}

/// Renders the intensity volume through `tf` to an RGBA image of the camera's viewport.
///
/// # Errors
/// Returns [`RenderError::InvalidParameter`] for a bad camera, transfer function or step,
/// and [`RenderError::ChannelMismatch`] when the transfer function does not have one entry
/// per volume channel.
pub fn render_volume(
    volume: &Volume,
    camera: &Camera,
    tf: &TransferFunction,
    settings: &RenderSettings,
) -> Result<Image<Rgba8>, RenderError> {
    check_inputs(volume, camera, tf, settings)?;
    let (w, h) = camera.viewport;
    let back = ortho_back(volume, camera);
    let mut data = Vec::with_capacity(w * h);
    for y in 0..h {
        for x in 0..w {
            let ray = camera.ray(x as f64 + 0.5, y as f64 + 0.5, back);
            let m = march(volume, None, &ray, tf, settings.step);
            data.push(finish(m.colour, m.alpha, settings.background));
        }
    }
    Ok(Image::new_2d(w, h, data)?)
}

/// Validates the inputs shared by every render entry point.
pub(super) fn check_inputs(
    volume: &Volume,
    camera: &Camera,
    tf: &TransferFunction,
    settings: &RenderSettings,
) -> Result<(), RenderError> {
    camera.validate()?;
    tf.validate()?;
    if !(settings.step.is_finite() && settings.step > 0.0) {
        return Err(RenderError::InvalidParameter("march step must be finite and positive"));
    }
    if tf.channels.len() != volume.channel_count() {
        return Err(RenderError::ChannelMismatch {
            expected: volume.channel_count(),
            actual: tf.channels.len(),
        });
    }
    Ok(())
}

/// The orthographic back offset that puts every point of the volume in front of the ray origin.
pub(super) fn ortho_back(volume: &Volume, camera: &Camera) -> f64 {
    let ext = volume.extent();
    let centre = [ext[0] / 2.0, ext[1] / 2.0, ext[2] / 2.0];
    let diag = length(ext);
    camera.distance + length(sub(camera.target, centre)) + diag
}

/// Result of marching one ray.
#[derive(Debug, Clone, PartialEq)]
pub(super) struct March {
    /// Premultiplied accumulated colour.
    pub colour: [f32; 3],
    /// Accumulated opacity.
    pub alpha: f32,
    /// Label of the first labelled voxel met along the ray, if labels were supplied.
    pub first_label: Option<u32>,
}

/// Marches `ray` through the volume. When `labels` is given, also records the first
/// nonzero label, which keeps being searched after opacity has saturated.
pub(super) fn march(
    volume: &Volume,
    labels: Option<&Image<u32>>,
    ray: &Ray,
    tf: &TransferFunction,
    step: f64,
) -> March {
    let mut out = March { colour: [0.0; 3], alpha: 0.0, first_label: None };
    let Some((t0, t1)) = intersect_box(ray, volume.extent()) else {
        return out;
    };
    let spacing = volume.spacing();
    let shape = volume.shape();
    let channels = volume.channels();
    let mut values = vec![0.0f32; channels.len()];
    let mut t = t0.max(0.0) + 0.5 * step;
    while t < t1 {
        let p = add(ray.origin, scale(ray.dir, t));
        let clipped = tf.clip.is_some_and(|c| !c.keeps(p));
        if !clipped {
            if out.alpha < OPACITY_CUTOFF {
                for (v, ch) in values.iter_mut().zip(channels) {
                    *v = sample_trilinear(ch, spacing, p);
                }
                let (c, a) = tf.combine(&values);
                let k = (1.0 - out.alpha) * a;
                for (acc, cv) in out.colour.iter_mut().zip(c) {
                    *acc += k * cv;
                }
                out.alpha += k;
            }
            if let (Some(lab), None) = (labels, out.first_label) {
                let id = lab.data()[nearest_voxel(shape, spacing, p)];
                if id != 0 {
                    out.first_label = Some(id);
                }
            }
        }
        let labels_done = labels.is_none() || out.first_label.is_some();
        if out.alpha >= OPACITY_CUTOFF && labels_done {
            break;
        }
        t += step;
    }
    out
}

/// Composites premultiplied `colour` with opacity `alpha` over `background` and quantises.
pub(super) fn finish(colour: [f32; 3], alpha: f32, background: [u8; 3]) -> Rgba8 {
    let rest = 1.0 - alpha;
    let mut out = [255u8; 4];
    for (o, (&c, &b)) in out.iter_mut().zip(colour.iter().zip(&background)) {
        *o = to_u8(c + rest * f32::from(b) / 255.0);
    }
    out
}

/// Parametric range `[t0, t1]` over which `ray` is inside the box `[0, extent]`.
fn intersect_box(ray: &Ray, extent: [f64; 3]) -> Option<(f64, f64)> {
    let mut t0 = f64::NEG_INFINITY;
    let mut t1 = f64::INFINITY;
    for ((&o, &d), &hi) in ray.origin.iter().zip(&ray.dir).zip(&extent) {
        if d.abs() < 1e-12 {
            if o < 0.0 || o > hi {
                return None;
            }
        } else {
            let inv = 1.0 / d;
            let (mut ta, mut tb) = ((0.0 - o) * inv, (hi - o) * inv);
            if ta > tb {
                std::mem::swap(&mut ta, &mut tb);
            }
            t0 = t0.max(ta);
            t1 = t1.min(tb);
            if t0 > t1 {
                return None;
            }
        }
    }
    Some((t0, t1))
}

#[allow(clippy::float_cmp)]
#[cfg(test)]
mod tests {
    use super::*;
    use crate::render3d::camera::Projection;
    use crate::render3d::transfer::ChannelTransfer;

    /// A uniform box of `n^3` voxels, spacing 1, all with value `value`.
    fn uniform_box(n: usize, value: f32) -> Volume {
        let img = Image::new_3d(n, n, n, vec![value; n * n * n]).unwrap();
        Volume::new(vec![img], [1.0; 3]).unwrap()
    }

    fn opaque_tf(colour: [f32; 3]) -> TransferFunction {
        let mut ch = ChannelTransfer::new(colour, (0.0, 1.0));
        ch.opacity = vec![(0.0, 1.0), (1.0, 1.0)];
        TransferFunction { channels: vec![ch], clip: None }
    }

    fn ortho(yaw: f64, pitch: f64, height: f64, viewport: (usize, usize)) -> Camera {
        Camera {
            target: [10.0, 10.0, 10.0],
            yaw,
            pitch,
            distance: 200.0,
            projection: Projection::Orthographic { height },
            viewport,
        }
    }

    #[test]
    fn render3d_raycast_head_on_box_is_transfer_colour() {
        let colour = [0.8, 0.25, 0.5];
        let vol = uniform_box(20, 1.0);
        let cam = ortho(0.0, 0.0, 40.0, (80, 80));
        let settings = RenderSettings { step: 0.5, background: [0, 0, 255] };
        let img = render_volume(&vol, &cam, &opaque_tf(colour), &settings).unwrap();
        // Viewport 80 px for 40 units: 2 px per unit, so the box covers pixels 20..60.
        let want = [to_u8(colour[0]), to_u8(colour[1]), to_u8(colour[2]), 255];
        for y in 20..60 {
            for x in 20..60 {
                assert_eq!(*img.get(x, y, 0).unwrap(), want, "pixel ({x}, {y})");
            }
        }
        // A corner pixel is outside the box and shows the background.
        assert_eq!(*img.get(0, 0, 0).unwrap(), [0, 0, 255, 255]);
    }

    #[test]
    fn render3d_raycast_rotated_box_covers_projected_area() {
        // Orthographic silhouette of an axis-aligned cube of side L has area
        // L^2 * (|dx| + |dy| + |dz|), where d is the unit view direction.
        let n = 20usize;
        let (yaw, pitch) = (0.5f64, 0.3f64);
        let vol = uniform_box(n, 1.0);
        let cam = ortho(yaw, pitch, 40.0, (120, 120));
        let settings = RenderSettings::new(0.25);
        let img = render_volume(&vol, &cam, &opaque_tf([1.0; 3]), &settings).unwrap();
        let covered = img.data().iter().filter(|p| p[0] != 0).count() as f64;
        let d = [pitch.cos() * yaw.cos(), pitch.cos() * yaw.sin(), pitch.sin()];
        let scale = 120.0 / 40.0;
        let side = n as f64;
        let expected = side * side * (d[0].abs() + d[1].abs() + d[2].abs()) * scale * scale;
        let rel = (covered - expected).abs() / expected;
        assert!(rel < 0.03, "covered {covered}, expected {expected}, rel {rel}");
    }

    #[test]
    fn render3d_raycast_clip_plane_discards_samples() {
        use crate::render3d::transfer::ClipPlane;
        let vol = uniform_box(20, 1.0);
        let cam = ortho(0.0, 0.0, 40.0, (80, 80));
        let settings = RenderSettings::new(0.5);
        // The camera looks along -x from +x. Keeping x >= 10 leaves the box visible.
        let mut tf = opaque_tf([1.0, 1.0, 1.0]);
        tf.clip = Some(ClipPlane { normal: [1.0, 0.0, 0.0], offset: -10.0 });
        let kept = render_volume(&vol, &cam, &tf, &settings).unwrap();
        assert_eq!(kept.get(40, 40, 0).unwrap(), &[255, 255, 255, 255]);
        // Keeping x >= 25 discards the whole box (which spans x in [0, 20]).
        tf.clip = Some(ClipPlane { normal: [1.0, 0.0, 0.0], offset: -25.0 });
        let none = render_volume(&vol, &cam, &tf, &settings).unwrap();
        assert_eq!(none.get(40, 40, 0).unwrap(), &[0, 0, 0, 255]);
    }

    #[test]
    fn render3d_raycast_rejects_channel_mismatch() {
        let vol = uniform_box(4, 1.0);
        let mut tf = opaque_tf([1.0; 3]);
        tf.channels.push(ChannelTransfer::new([1.0; 3], (0.0, 1.0)));
        let cam = ortho(0.0, 0.0, 8.0, (8, 8));
        let err = render_volume(&vol, &cam, &tf, &RenderSettings::new(1.0)).unwrap_err();
        assert_eq!(err, RenderError::ChannelMismatch { expected: 1, actual: 2 });
    }
}
