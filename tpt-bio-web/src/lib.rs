//! Browser demo for tpt-bio. It renders one scene with the CPU ray caster and with WebGPU,
//! then reports the largest per-channel difference between the two images.
//!
//! The crate builds only for wasm32. Build it with `wasm-bindgen` as described in `README.md`.
#![cfg(target_arch = "wasm32")]
// Grid coordinates (x, y, z) and the scene's w, h, d read most clearly as single letters.
#![allow(clippy::many_single_char_names)]

use js_sys::{Object, Reflect};
use tpt_bio_vision::image::Image;
use tpt_bio_vision::render3d::camera::{Camera, Projection};
use tpt_bio_vision::render3d::raycast::{RenderSettings, render_volume};
use tpt_bio_vision::render3d::transfer::{ChannelTransfer, TransferFunction};
use tpt_bio_vision::render3d::volume::Volume;
use tpt_bio_vision::render3d::webgpu::WebGpuRenderer;
use wasm_bindgen::prelude::*;

/// Largest per-channel difference, out of 255, that the check accepts. The GPU uses a
/// lookup table for transfer functions, so exact equality is not expected.
const TOLERANCE: i32 = 3;

/// A two-channel scene: a bright ball and a dim shell on an odd-sized grid.
fn scene(width: usize, height: usize) -> (Volume, Camera, TransferFunction, RenderSettings) {
    let (w, h, d) = (13usize, 11usize, 9usize);
    let mut a = Vec::with_capacity(w * h * d);
    let mut b = Vec::with_capacity(w * h * d);
    for z in 0..d {
        for y in 0..h {
            for x in 0..w {
                let r2 =
                    (x as f32 - 6.0).powi(2) + (y as f32 - 5.0).powi(2) + (z as f32 - 4.0).powi(2);
                a.push(200.0 * (-r2 / 12.0).exp());
                b.push(90.0 * (-(r2 - 20.0).powi(2) / 30.0).exp());
            }
        }
    }
    let volume = Volume::new(
        vec![
            Image::new_3d(w, h, d, a).expect("shape matches"),
            Image::new_3d(w, h, d, b).expect("shape matches"),
        ],
        [1.0, 1.2, 0.8],
    )
    .expect("valid volume");
    let mut red = ChannelTransfer::new([1.0, 0.2, 0.1], (0.0, 200.0));
    red.opacity = vec![(0.0, 0.0), (0.5, 0.4), (1.0, 0.9)];
    let mut green = ChannelTransfer::new([0.1, 1.0, 0.3], (0.0, 90.0));
    green.opacity = vec![(0.0, 0.0), (1.0, 0.3)];
    green.gain = 0.8;
    let tf = TransferFunction { channels: vec![red, green], clip: None };
    let camera = Camera {
        target: [6.0, 5.0, 4.0],
        yaw: 0.6,
        pitch: 0.3,
        distance: 30.0,
        projection: Projection::Perspective { fov_y: 0.9 },
        viewport: (width, height),
    };
    (volume, camera, tf, RenderSettings { step: 0.5, background: [10, 20, 30] })
}

fn js_error(e: impl std::fmt::Display) -> JsValue {
    JsValue::from_str(&e.to_string())
}

/// Renders the demo scene on the CPU and with WebGPU, and compares the two.
///
/// Resolves to an object with `ok` (true when the largest difference is within tolerance),
/// `maxDiff`, `pixels`, `width` and `height`.
///
/// # Errors
/// Rejects with a message when WebGPU is unavailable or a render step fails.
#[wasm_bindgen]
pub async fn run_comparison(width: u32, height: u32) -> Result<JsValue, JsValue> {
    let (width, height) = (width.max(1) as usize, height.max(1) as usize);
    let (volume, camera, tf, settings) = scene(width, height);
    let cpu = render_volume(&volume, &camera, &tf, &settings).map_err(js_error)?;
    let gpu = WebGpuRenderer::new().await.map_err(js_error)?;
    let got = gpu.render(&volume, &camera, &tf, &settings).await.map_err(js_error)?;

    let max_diff = cpu
        .data()
        .iter()
        .zip(got.data())
        .flat_map(|(a, b)| (0..4).map(move |c| (i32::from(a[c]) - i32::from(b[c])).abs()))
        .max()
        .unwrap_or(0);

    let out = Object::new();
    let set = |k: &str, v: JsValue| Reflect::set(&out, &JsValue::from_str(k), &v).map(|_| ());
    set("ok", JsValue::from_bool(max_diff <= TOLERANCE))?;
    set("maxDiff", JsValue::from_f64(f64::from(max_diff)))?;
    set("pixels", JsValue::from_f64((width * height) as f64))?;
    set("width", JsValue::from_f64(width as f64))?;
    set("height", JsValue::from_f64(height as f64))?;
    Ok(out.into())
}
