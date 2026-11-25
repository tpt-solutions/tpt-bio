//! Label overlay drawn over the intensity rendering.
//!
//! A label volume (`Image<u32>`, `0` = background) with the same shape as the intensity
//! volume is blended over the ray-cast image. Each ray records the first nonzero label it
//! meets while marching. That is the first-hit surface: the ray's first labelled voxel,
//! found whether or not intensity opacity has already saturated. The pixel becomes
//! `(1 - opacity) * intensity + opacity * label_colour(id)`. Pixels whose ray meets no
//! label keep the plain intensity colour.
//!
//! Label colours come from [`label_colour`], a fixed hash of the label id, so the same id
//! always gets the same colour in every render.

use super::camera::Camera;
use super::raycast::{RenderSettings, check_inputs, finish, march, ortho_back, par_pixels};
use super::transfer::TransferFunction;
use super::volume::Volume;
use super::{RenderError, Rgba8};
use crate::image::Image;

/// How strongly labels are blended over the intensity image.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct LabelOverlay {
    /// Blend weight in `[0, 1]`: `0` hides labels, `1` shows only label colours.
    pub opacity: f32,
}

/// Deterministic RGB colour for a label id.
///
/// The id is passed through the `SplitMix64` finaliser and each byte is lifted into
/// `[64, 255]`, so colours are never dark. Background (`0`) gets a colour too, but the
/// overlay never draws it.
#[must_use]
pub fn label_colour(id: u32) -> [u8; 3] {
    let mut z = u64::from(id).wrapping_add(0x9E37_79B9_7F4A_7C15);
    z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
    z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
    z ^= z >> 31;
    let lift = |b: u64| (64 + (b & 0xFF) * 191 / 255) as u8;
    [lift(z), lift(z >> 8), lift(z >> 16)]
}

/// Renders the intensity volume with the label overlay drawn on top.
///
/// # Errors
/// Same as [`super::raycast::render_volume`], plus [`RenderError::Image`] when `labels`
/// does not have the volume's shape and [`RenderError::InvalidParameter`] when the
/// opacity is outside `[0, 1]`.
pub fn render_labelled(
    volume: &Volume,
    labels: &Image<u32>,
    camera: &Camera,
    tf: &TransferFunction,
    settings: &RenderSettings,
    overlay: &LabelOverlay,
) -> Result<Image<Rgba8>, RenderError> {
    check_inputs(volume, camera, tf, settings)?;
    let channel = volume.channel(0).ok_or(RenderError::Empty)?;
    channel.check_same_shape(labels)?;
    if !(0.0..=1.0).contains(&overlay.opacity) {
        return Err(RenderError::InvalidParameter("label opacity must be in [0, 1]"));
    }
    let (w, h) = camera.viewport;
    let back = ortho_back(volume, camera);
    let data = par_pixels(w, h, |x, y| {
        let ray = camera.ray(x as f64 + 0.5, y as f64 + 0.5, back);
        let m = march(volume, Some(labels), &ray, tf, settings.step);
        let mut rgb = [0.0f32; 3];
        let base = finish(m.colour, m.alpha, settings.background);
        for (c, b) in rgb.iter_mut().zip(base) {
            *c = f32::from(b) / 255.0;
        }
        if let Some(id) = m.first_label {
            rgb = blend(rgb, label_colour(id), overlay.opacity);
        }
        quantise(rgb)
    });
    Ok(Image::new_2d(w, h, data)?)
}

/// Blends `label` over `base` with weight `opacity`.
#[must_use]
pub fn blend(base: [f32; 3], label: [u8; 3], opacity: f32) -> [f32; 3] {
    let mut out = base;
    for (o, l) in out.iter_mut().zip(label) {
        *o = *o * (1.0 - opacity) + f32::from(l) / 255.0 * opacity;
    }
    out
}

fn quantise(rgb: [f32; 3]) -> Rgba8 {
    let q = |v: f32| (v.clamp(0.0, 1.0) * 255.0).round() as u8;
    [q(rgb[0]), q(rgb[1]), q(rgb[2]), 255]
}

#[allow(clippy::float_cmp)]
#[cfg(test)]
mod tests {
    use super::*;
    use crate::render3d::camera::Projection;
    use crate::render3d::transfer::ChannelTransfer;

    #[test]
    fn render3d_overlay_palette_is_deterministic_and_distinct() {
        assert_eq!(label_colour(7), label_colour(7));
        let ids = [1u32, 2, 3, 4, 5, 100];
        for (i, a) in ids.iter().enumerate() {
            for b in &ids[i + 1..] {
                assert_ne!(label_colour(*a), label_colour(*b), "ids {a} and {b}");
            }
        }
    }

    #[test]
    fn render3d_overlay_two_labels_render_two_colours() {
        // Two label blocks stacked along z (the camera's up axis at yaw 0). Seen head-on from
        // +x, the first labelled voxel is on the near face of whichever block is in front.
        let n = 12usize;
        let mut labels = vec![0u32; n * n * n];
        let vals = vec![0.5f32; n * n * n];
        for z in 0..n {
            for y in 0..n {
                for x in 0..n {
                    let i = (z * n + y) * n + x;
                    labels[i] = if z < n / 2 { 1 } else { 2 };
                }
            }
        }
        let vol = Volume::new(vec![Image::new_3d(n, n, n, vals).unwrap()], [1.0; 3]).unwrap();
        let labels = Image::new_3d(n, n, n, labels).unwrap();
        let mut ch = ChannelTransfer::new([0.5, 0.5, 0.5], (0.0, 1.0));
        ch.opacity = vec![(0.0, 0.0), (1.0, 0.0)];
        let tf = TransferFunction { channels: vec![ch], clip: None };
        let cam = Camera {
            target: [6.0, 6.0, 6.0],
            yaw: 0.0,
            pitch: 0.0,
            distance: 100.0,
            projection: Projection::Orthographic { height: 12.0 },
            viewport: (48, 48),
        };
        let opts = LabelOverlay { opacity: 1.0 };
        let img =
            render_labelled(&vol, &labels, &cam, &tf, &RenderSettings::new(0.5), &opts).unwrap();
        // Top rows look at high z (label 2); bottom rows at low z (label 1).
        let top = *img.get(24, 4, 0).unwrap();
        let bottom = *img.get(24, 43, 0).unwrap();
        assert_eq!(top, [label_colour(2)[0], label_colour(2)[1], label_colour(2)[2], 255]);
        assert_eq!(bottom, [label_colour(1)[0], label_colour(1)[1], label_colour(1)[2], 255]);
        assert_ne!(top, bottom);
    }

    #[test]
    fn render3d_overlay_zero_opacity_is_plain_intensity() {
        let n = 4usize;
        let vol =
            Volume::new(vec![Image::new_3d(n, n, n, vec![1.0; n * n * n]).unwrap()], [1.0; 3])
                .unwrap();
        let labels = Image::new_3d(n, n, n, vec![3u32; n * n * n]).unwrap();
        let mut ch = ChannelTransfer::new([0.2, 0.4, 0.6], (0.0, 1.0));
        ch.opacity = vec![(0.0, 1.0), (1.0, 1.0)];
        let tf = TransferFunction { channels: vec![ch], clip: None };
        let cam = Camera {
            target: [2.0; 3],
            yaw: 0.0,
            pitch: 0.0,
            distance: 50.0,
            projection: Projection::Orthographic { height: 4.0 },
            viewport: (8, 8),
        };
        let settings = RenderSettings::new(0.5);
        let plain = crate::render3d::raycast::render_volume(&vol, &cam, &tf, &settings).unwrap();
        let off =
            render_labelled(&vol, &labels, &cam, &tf, &settings, &LabelOverlay { opacity: 0.0 })
                .unwrap();
        assert_eq!(plain, off);
    }

    #[test]
    fn render3d_overlay_blend_endpoints() {
        assert_eq!(blend([0.5; 3], [255, 0, 0], 0.0), [0.5; 3]);
        let full = blend([0.5; 3], [255, 0, 0], 1.0);
        assert_eq!(full, [1.0, 0.0, 0.0]);
    }
}
