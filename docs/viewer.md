# 3D volume rendering

`tpt-bio-vision::render3d` is a CPU reference renderer for multi-channel 3D microscopy
volumes. It turns a volume into an RGBA image for one camera position, and it can
extract and rasterise an isosurface mesh. It is the reference that a GPU renderer would
be checked against.

**What this is not.** There is no GPU code, no WebGPU or wgpu, no browser build and no
desktop application. There is no window and no interactive control loop. The `tptbio`
command line does not render either: `tptbio view` prints a data summary of an image.
The GPU viewer described in spec ADR-003 is not started (todo.md, Phase 4).

## What is in the module

| Module | What it does |
| --- | --- |
| `volume` | `Volume`: f32 channels with physical voxel spacing, trilinear sampling, `from_u16` |
| `camera` | `Camera` with a perspective or orthographic `Projection`, orbit angles, `project` and `ray` |
| `transfer` | `ChannelTransfer` (colour knots, opacity knots, window, gain) and `TransferFunction` with an optional `ClipPlane` |
| `raycast` | `render_volume`: front-to-back alpha compositing along each pixel's ray |
| `overlay` | `render_labelled`: label colours blended over the volume, first-hit surface, deterministic `label_colour` |
| `isosurface` | `isosurface`: marching tetrahedra to a triangle `Mesh`; `label_indicator` for one label |
| `rasterise` | `rasterise_mesh`: flat-shaded z-buffer rasteriser for a `Mesh` |
| `lod` | `downsample_max` (block maximum) and `choose_level` from the screen footprint |

All output is `Image<Rgba8>` (`[u8; 4]` per pixel). Errors are `RenderError`.

## Sample data

The example needs a volume. `make_sample_volume` writes a synthetic one: 96 x 96 x 32
voxels with three 16-bit channels (nuclei, membrane and a diffuse marker) and six cells.
It is saved as OME-TIFF, the format `OmeTiff` reads.

```sh
cargo run --release -p tpt-bio-vision --example make_sample_volume -- sample-volume.ome.tif
```

## Render the sample

```sh
cargo run --release -p tpt-bio-vision --example render3d_viewer -- sample-volume.ome.tif viewer-out
```

The run prints the volume shape, the number of nuclei found and the size of the mesh,
then writes three binary PPM images (`P6`, 320 x 240) to `viewer-out/`:

- `volume.ppm`: the three channels composited with per-channel colours and windows.
- `labels.ppm`: the same view with the 3D-labelled nuclei overlaid at 60 % opacity.
- `mesh.ppm`: the nuclei isosurface at intensity 1500, flat-shaded.

Expected output for the sample: `nuclei objects: 6`, and a mesh of about 9,500
vertices and 19,000 triangles. Use `--release`; the ray caster is slow in debug builds.

PPM is a plain format that most image viewers and converters read. It is not a viewer.

## Code

The core of the example is the camera, the transfer function and one render call. The
excerpt below is from `tpt-bio-vision/examples/render3d_viewer.rs`, unchanged:

```rust
let volume = Volume::from_u16(&channels, [1.0, 1.0, 2.0])?;
let tf = TransferFunction { channels: transfers, clip: None };
let camera = Camera {
    target: centre,
    yaw: 0.6,
    pitch: 0.4,
    distance: 2.6 * extent.iter().copied().fold(0.0, f64::max),
    projection: Projection::Perspective { fov_y: 0.8 },
    viewport: (320, 240),
};
let settings = RenderSettings::new(1.0);
let image = render_volume(&volume, &camera, &tf, &settings)?;
```

Parameters that matter:

- `Volume::from_u16` takes the spacing in each axis. Anisotropic microscopy stacks
  have a larger z spacing. The renderer uses it for the volume's extent and for
  sampling; the march step is in world units, not voxels.
- `ChannelTransfer::new(colour, (lo, hi))` maps intensity linearly to opacity between
  `lo` and `hi`. `volume.intensity_range(c)` gives the data range.
- `RenderSettings::step` is the march step in world units. Smaller steps look better
  and cost more time.
- `TransferFunction::clip` discards samples on one side of a plane, for example to
  cut a cell open.
- `Camera::validate` checks the camera; the render functions call it for you.

For a label overlay, call `label_3d` on a mask, then `render_labelled`. For a surface,
call `isosurface` on an `Image<f32>` field with an iso value, then `rasterise_mesh`.

## Limits

- **CPU only, and slow for large volumes.** The renderer marches every pixel's ray
  through the whole volume. Use `downsample_max` and `choose_level` to render a
  coarser level, but nothing does this automatically.
- **No chunked or streamed volumes.** A `Volume` is held in memory as f32 per channel.
  The todo item for chunked streaming is not done.
- **Flat shading only.** The mesh rasteriser colours each triangle with one colour; there
  are no normals, lights or textures.
- **Labels are blended by first hit.** The overlay shows the label at the first
  surface a ray meets, not per-voxel label colour along the ray.
- **No interaction.** The camera is set in code. A caller who wants orbiting must
  build a new `Camera` for each frame.
- **Sampling is trilinear, and the isosurface is marching tetrahedra.** Both are
  deterministic, and the mesh is not decimated.
