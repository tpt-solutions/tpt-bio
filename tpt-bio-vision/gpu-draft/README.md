# GPU ray casting (retired native path, kept as reference)

`raycast.wgsl` is the compute shader. It is standard WGSL, so a browser can run it unchanged
through the WebGPU API. It is not compiled by this crate.

`wgpu_host.rs.txt` is the native wgpu 26 host code that drove the shader. It is not compiled.
It was verified on an Intel UHD 770 (Windows, DX12): perspective two-channel and orthographic
clip-plane renders matched the CPU renderer to within 3 of 255 per channel.

Retired, because wgpu cannot meet the licence policy. Its tree includes `codespan-reporting`,
which is an Apache-2.0-only non-optional dependency of naga, so no feature selection removes
it. `gl_generator` and `glutin_wgl_sys` are also Apache-2.0 only, `hexf-parse` is CC0-1.0,
and `paste` has advisory RUSTSEC-2024-0436.

The replacement for the browser is WebGPU through `web-sys` and `wasm-bindgen` (both
`MIT OR Apache-2.0`), which call the browser's own WebGPU implementation and need no native
GPU crates. Native GPU rendering is not planned. The CPU renderer and streamed levels cover
large volumes.
