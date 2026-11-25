# 3D viewer architecture: sharing a paradigm with tpt-vertex

Status: open. This note records what is known on the tpt-bio side and lists the
questions that need answers. It does not decide the design.

## What could not be verified

The todo item asks for a note on sharing the rendering paradigm with `tpt-vertex`
(spec ADR-003). `tpt-vertex` is not in this workspace, and its source is not available
in this environment. Its design, its crate boundaries, its GPU abstraction and its
license therefore could not be checked. Nothing below describes `tpt-vertex`. Every
statement about it is an open question.

## The tpt-bio side

ADR-003 says 3D rendering will use wgpu (WebGPU) and share the paradigm of `tpt-vertex`.
wgpu has since been rejected by the licence policy (`docs/dependency-audit.md`), and the
ADR carries a note saying so. Today the repository has the following, and nothing else:

- **No native GPU dependency.** `wgpu` is not in any `Cargo.toml` and will not be. The
  browser path uses `web-sys` WebGPU (feature `webgpu`, wasm32 only). The planned native
  backend is Vulkan through `ash` with an offline-compiled SPIR-V shader, which passed the
  licence audit. It is not started.
- **A CPU reference renderer** in `tpt-bio-vision::render3d` (see `docs/viewer.md`).
  Its public types are the boundary that any GPU path would have to meet:

  | Type | Shape of the data |
  | --- | --- |
  | `Volume` | f32 channels, one `Image<f32>` each, with `[f64; 3]` spacing |
  | `Camera` | target `[f64; 3]`, `yaw`, `pitch`, `distance`, `Projection`, viewport; Z-up world |
  | `TransferFunction` | per channel: colour and opacity knots, a window and a gain; optional clip plane |
  | `Mesh` | `vertices: Vec<[f32; 3]>`, `triangles: Vec<[u32; 3]>` (counter-clockwise from outside) |
  | output | `Image<Rgba8>`, one `[u8; 4]` per pixel |

- **Behaviour the GPU version would have to reproduce**, if it is meant to match the
  reference: trilinear sampling, front-to-back compositing with a step in world units,
  first-hit label blending, flat shading in the rasteriser, and the orbit camera's
  convention (zero yaw looks along `-x`, pitch measured up from the `xy` plane).

The CPU renderer runs in `f64` for camera geometry and `f32` for samples, and writes
8-bit colour. A GPU path would likely use `f32` throughout, so comparisons need a
tolerance.

## Questions to answer when tpt-vertex is available

1. **Licence and availability.** Is it MIT or Apache-2.0 (or dual), and can it be
   depended on from this workspace, or only copied? The `deny.toml` policy applies to
   its dependency tree too.
2. **Crate boundary.** Is it a library crate, an application, or a set of modules? Which
   parts would tpt-bio depend on, and which would be duplicated?
3. **Retained or immediate.** Does it keep a scene graph between frames, or rebuild
   draw calls each frame? A volume renderer with a large texture needs the first kind of
   resource lifetime.
4. **GPU abstraction.** Does it use wgpu directly, and which version? Are pipeline and
   bind-group layouts shared, or built per application?
5. **Coordinate conventions.** Z-up or Y-up, right-handed or left-handed, and what depth
   range is used (`0..1` or `-1..1`)? The tpt-bio camera is Z-up and right-handed.
6. **Camera model.** Does it take a matrix, or orbit parameters like `Camera`? Which
   projections does it support? Is there an orthographic mode for scale bars?
7. **Transfer functions and colour maps.** How are they represented on the GPU: a 1D
   lookup texture, shader code, or per-draw uniforms? A 1D lookup would match the
   knot-based `ChannelTransfer`.
8. **Volume textures.** How does it upload 3D data, which formats (f32, u16, R8), and does
   it stream or chunk large volumes? The tpt-bio side has no chunked path yet.
9. **Shaders.** Is it WGSL only? Are shaders shared with the browser build?
10. **Targets.** Does it target WebAssembly and WebGPU in the browser, and does it run
    under a Tauri desktop shell? The spec names both.
11. **Input and events.** Does it own the event loop, and how do camera controls reach
    the application?
12. **Testing.** Does it have a headless or reference rendering path that tpt-bio could
    use to check its GPU output against `render3d`?
13. **Precision and units.** `f32` or `f64` in the public API, and how it names
    physical units.

Until these are answered, the tpt-bio work does not depend on tpt-vertex. The CPU
renderer and its public types stay as they are. A native GPU path is added behind the
`RenderBackend` trait in `render3d::backend`, whose only implementation today is
`CpuBackend`.
