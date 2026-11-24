//! Browser ray casting with WebGPU, behind the `webgpu` feature (wasm32 only).
//!
//! [`WebGpuRenderer::render`] produces the same image as [`super::raycast::render_volume`].
//! It uses the browser's WebGPU through `web-sys`, so no native GPU library is linked. Ray
//! setup and box clipping run on the CPU with the same code as the CPU renderer. Marching runs
//! in the compute shader `raycast.wgsl`, which is standard WGSL that the browser compiles.
//!
//! Each channel's transfer function is evaluated on the CPU into a lookup table of
//! [`LUT_LEN`] entries across its window, and the shader interpolates between entries. The
//! result therefore matches the CPU renderer to within the table spacing.
//!
//! Volumes are uploaded as one storage buffer. WebGPU's default maximum storage binding is
//! 128 MiB, so larger volumes need a streamed level from [`super::stream`].
//!
//! The `web-sys` WebGPU bindings are gated behind `--cfg=web_sys_unstable_apis`, which the
//! workspace sets for wasm32 in `.cargo/config.toml`.

use std::fmt;

use js_sys::Float32Array;
use wasm_bindgen::{JsCast, JsValue};
use wasm_bindgen_futures::JsFuture;
use web_sys::{
    GpuAdapter, GpuAutoLayoutMode, GpuBindGroupDescriptor, GpuBindGroupEntry, GpuBufferBinding,
    GpuBufferDescriptor, GpuComputePipeline, GpuComputePipelineDescriptor, GpuDevice,
    GpuProgrammableStage, GpuQueue, GpuShaderModuleDescriptor,
};

use super::RenderError;
use super::Rgba8;
use super::camera::Camera;
use super::raycast::{
    OPACITY_CUTOFF, RenderSettings, check_inputs, finish, intersect_box, ortho_back,
};
use super::transfer::TransferFunction;
use super::volume::Volume;
use crate::image::Image;

/// Entries in each channel's transfer lookup table.
pub const LUT_LEN: usize = 4096;

const SHADER: &str = include_str!("raycast.wgsl");
const WORKGROUP: u32 = 8;
/// WebGPU's default `maxStorageBufferBindingSize`.
const DEFAULT_MAX_BINDING: u64 = 128 * 1024 * 1024;

// WebGPU flag values from the specification (GPUBufferUsage, GPUShaderStage, GPUMapMode).
const USAGE_MAP_READ: u32 = 0x0001;
const USAGE_COPY_SRC: u32 = 0x0004;
const USAGE_COPY_DST: u32 = 0x0008;
const USAGE_UNIFORM: u32 = 0x0040;
const USAGE_STORAGE: u32 = 0x0080;
const MAP_MODE_READ: u32 = 0x0001;

/// Errors from browser ray casting.
#[derive(Debug, Clone, PartialEq)]
pub enum WebGpuError {
    /// The browser has no WebGPU, or no adapter was returned.
    NoAdapter(String),
    /// The device or pipeline could not be created.
    Device(String),
    /// A buffer, dispatch or readback failed. The message says which step.
    Gpu(String),
    /// The scene is invalid, or too large for one storage binding.
    Render(RenderError),
}

impl fmt::Display for WebGpuError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::NoAdapter(m) => write!(f, "no WebGPU adapter: {m}"),
            Self::Device(m) => write!(f, "WebGPU device: {m}"),
            Self::Gpu(m) => write!(f, "WebGPU: {m}"),
            Self::Render(e) => write!(f, "{e}"),
        }
    }
}

impl std::error::Error for WebGpuError {}

impl From<RenderError> for WebGpuError {
    fn from(e: RenderError) -> Self {
        Self::Render(e)
    }
}

fn js_err(step: &str, e: &JsValue) -> WebGpuError {
    WebGpuError::Gpu(format!("{step}: {}", e.as_string().unwrap_or_else(|| format!("{e:?}"))))
}

/// A WebGPU device with the ray-marching pipeline ready.
pub struct WebGpuRenderer {
    device: GpuDevice,
    queue: GpuQueue,
    pipeline: GpuComputePipeline,
}

impl fmt::Debug for WebGpuRenderer {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("WebGpuRenderer").finish_non_exhaustive()
    }
}

impl WebGpuRenderer {
    /// Requests an adapter and device from the browser and builds the pipeline.
    ///
    /// # Errors
    /// Returns [`WebGpuError::NoAdapter`] when the browser has no WebGPU or no adapter, and
    /// [`WebGpuError::Device`] or [`WebGpuError::Gpu`] when device or pipeline creation fails.
    pub async fn new() -> Result<Self, WebGpuError> {
        let window = web_sys::window().ok_or_else(|| WebGpuError::NoAdapter("no window".into()))?;
        let gpu = window.navigator().gpu();
        let adapter_value = JsFuture::from(gpu.request_adapter())
            .await
            .map_err(|e| js_err("requestAdapter", &e))?;
        if adapter_value.is_null() || adapter_value.is_undefined() {
            return Err(WebGpuError::NoAdapter("requestAdapter returned null".into()));
        }
        let adapter: GpuAdapter = adapter_value
            .dyn_into()
            .map_err(|_| WebGpuError::NoAdapter("requestAdapter returned a non-adapter".into()))?;
        let device: GpuDevice = JsFuture::from(adapter.request_device())
            .await
            .map_err(|e| WebGpuError::Device(format!("{e:?}")))?
            .dyn_into()
            .map_err(|_| WebGpuError::Device("requestDevice returned a non-device".into()))?;
        let queue = device.queue();

        let module = device.create_shader_module(&GpuShaderModuleDescriptor::new(SHADER));
        let stage = GpuProgrammableStage::new(&module);
        stage.set_entry_point("main");
        let desc = GpuComputePipelineDescriptor::new_with_gpu_auto_layout_mode(
            GpuAutoLayoutMode::Auto,
            &stage,
        );
        let pipeline = device.create_compute_pipeline(&desc);
        Ok(Self { device, queue, pipeline })
    }

    /// Renders `volume` through `tf` from `camera`, like [`super::raycast::render_volume`].
    ///
    /// # Errors
    /// Returns [`WebGpuError::Render`] for invalid inputs or a volume over the storage binding
    /// limit, and [`WebGpuError::Gpu`] when a buffer, dispatch or readback fails.
    pub async fn render(
        &self,
        volume: &Volume,
        camera: &Camera,
        tf: &TransferFunction,
        settings: &RenderSettings,
    ) -> Result<Image<Rgba8>, WebGpuError> {
        check_inputs(volume, camera, tf, settings)?;
        let (w, h) = camera.viewport;
        let channels = volume.channel_count();
        let (dw, dh, dd) = volume.shape();
        let voxels = dw * dh * dd;
        if (voxels * channels * 4) as u64 > DEFAULT_MAX_BINDING {
            return Err(WebGpuError::Render(RenderError::InvalidParameter(
                "volume exceeds the WebGPU storage binding limit; use a streamed level",
            )));
        }

        let back = ortho_back(volume, camera);
        let mut rays = Vec::with_capacity(w * h * 2);
        for y in 0..h {
            for x in 0..w {
                let ray = camera.ray(x as f64 + 0.5, y as f64 + 0.5, back);
                let (t0, t1) = intersect_box(&ray, volume.extent()).unwrap_or((0.0, 0.0));
                let o = ray.origin;
                let d = ray.dir;
                rays.push([o[0] as f32, o[1] as f32, o[2] as f32, t0.max(0.0) as f32]);
                rays.push([d[0] as f32, d[1] as f32, d[2] as f32, t1 as f32]);
            }
        }

        let mut voxel_data = Vec::with_capacity(voxels * channels);
        for c in 0..channels {
            voxel_data.extend_from_slice(volume.channel(c).ok_or(RenderError::ChannelIndex(c))?.data());
        }

        let mut lut = Vec::with_capacity(channels * LUT_LEN);
        let mut windows = Vec::with_capacity(channels * 2);
        for ch in &tf.channels {
            let (lo, hi) = ch.window;
            windows.extend_from_slice(&[lo, hi]);
            for i in 0..LUT_LEN {
                let v = lo + (hi - lo) * i as f32 / (LUT_LEN - 1) as f32;
                let (c, a) = ch.evaluate(v);
                lut.extend_from_slice(&[c[0], c[1], c[2], a]);
            }
        }

        let spacing = volume.spacing();
        let (clip_on, clip) = match tf.clip {
            Some(p) => (
                1u32,
                [p.normal[0] as f32, p.normal[1] as f32, p.normal[2] as f32, p.offset as f32],
            ),
            None => (0, [0.0; 4]),
        };
        let params = pack_params(
            [dw as u32, dh as u32, dd as u32, channels as u32],
            [
                spacing[0] as f32,
                spacing[1] as f32,
                spacing[2] as f32,
                settings.step as f32,
            ],
            [LUT_LEN as u32, clip_on, w as u32, h as u32],
            clip,
            OPACITY_CUTOFF,
        );

        let uniform = self.upload(&params, USAGE_UNIFORM)?;
        let vol_buf = self.upload(&f32_bytes(&voxel_data), USAGE_STORAGE)?;
        let ray_values: Vec<f32> = rays.iter().flatten().copied().collect();
        let rays_buf = self.upload(&f32_bytes(&ray_values), USAGE_STORAGE)?;
        let lut_buf = self.upload(&f32_bytes(&lut), USAGE_STORAGE)?;
        let win_buf = self.upload(&f32_bytes(&windows), USAGE_STORAGE)?;

        let out_size = (w * h * 16) as u32;
        let out_buf = self.create_buffer(out_size, USAGE_STORAGE | USAGE_COPY_SRC)?;
        let staging = self.create_buffer(out_size, USAGE_MAP_READ | USAGE_COPY_DST)?;

        let layout = self.pipeline.get_bind_group_layout(0);
        let entries = [
            bind(0, &uniform),
            bind(1, &vol_buf),
            bind(2, &rays_buf),
            bind(3, &lut_buf),
            bind(4, &win_buf),
            bind(5, &out_buf),
        ];
        let bind_group = self
            .device
            .create_bind_group(&GpuBindGroupDescriptor::new(&entries, &layout));

        let encoder = self.device.create_command_encoder();
        let pass = encoder.begin_compute_pass();
        pass.set_pipeline(&self.pipeline);
        pass.set_bind_group(0, Some(&bind_group));
        pass.dispatch_workgroups_with_workgroup_count_y(
            (w as u32).div_ceil(WORKGROUP),
            (h as u32).div_ceil(WORKGROUP),
        );
        pass.end();
        encoder
            .copy_buffer_to_buffer_with_u32_and_u32_and_u32(&out_buf, 0, &staging, 0, out_size)
            .map_err(|e| js_err("copyBufferToBuffer", &e))?;
        let command = encoder.finish();
        self.queue.submit(&[command]);

        let floats = read_back(&staging, out_size).await?;
        let pixels: Vec<Rgba8> = floats
            .chunks_exact(4)
            .map(|px| finish([px[0], px[1], px[2]], px[3], settings.background))
            .collect();
        Image::new_2d(w, h, pixels).map_err(|e| WebGpuError::Render(e.into()))
    }

    /// Creates a storage or uniform buffer holding `bytes`.
    fn upload(&self, bytes: &[u8], usage: u32) -> Result<web_sys::GpuBuffer, WebGpuError> {
        let buf = self.create_buffer(bytes.len() as u32, usage | USAGE_COPY_DST)?;
        self.queue
            .write_buffer_with_u32_and_u8_slice(&buf, 0, bytes)
            .map_err(|e| js_err("writeBuffer", &e))?;
        Ok(buf)
    }

    fn create_buffer(&self, size: u32, usage: u32) -> Result<web_sys::GpuBuffer, WebGpuError> {
        self.device
            .create_buffer(&GpuBufferDescriptor::new(size, usage))
            .map_err(|e| js_err("createBuffer", &e))
    }
}

/// Maps a staging buffer after the copy has run, and returns its `f32` contents.
async fn read_back(staging: &web_sys::GpuBuffer, size: u32) -> Result<Vec<f32>, WebGpuError> {
    JsFuture::from(staging.map_async_with_u32_and_u32(MAP_MODE_READ, 0, size))
        .await
        .map_err(|e| js_err("mapAsync", &e))?;
    let mapped = staging.get_mapped_range().map_err(|e| js_err("getMappedRange", &e))?;
    let floats = Float32Array::new(&mapped).to_vec();
    staging.unmap();
    Ok(floats)
}

fn bind(binding: u32, buffer: &web_sys::GpuBuffer) -> GpuBindGroupEntry {
    GpuBindGroupEntry::new_with_gpu_buffer_binding(binding, &GpuBufferBinding::new(buffer))
}

/// Packs the five `vec4` fields of the shader's `Params` struct, 80 bytes in total.
fn pack_params(
    dims: [u32; 4],
    spacing_step: [f32; 4],
    counts: [u32; 4],
    clip: [f32; 4],
    cutoff: f32,
) -> Vec<u8> {
    let mut out = Vec::with_capacity(80);
    for v in dims {
        out.extend_from_slice(&v.to_le_bytes());
    }
    for v in spacing_step {
        out.extend_from_slice(&v.to_le_bytes());
    }
    for v in counts {
        out.extend_from_slice(&v.to_le_bytes());
    }
    for v in clip {
        out.extend_from_slice(&v.to_le_bytes());
    }
    for v in [cutoff, 0.0, 0.0, 0.0] {
        out.extend_from_slice(&v.to_le_bytes());
    }
    out
}

fn f32_bytes(values: &[f32]) -> Vec<u8> {
    values.iter().flat_map(|v| v.to_le_bytes()).collect()
}
