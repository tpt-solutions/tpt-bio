# tpt-bio-vision

Biological image I/O, segmentation and feature extraction for tpt-bio.

`tpt-bio-vision` reads microscopy and volumetric image formats, runs
segmentation and filtering on 2D and 3D samples, measures per-object features,
and includes a CPU reference renderer for 3D volumes. Filters run on
`tpt-simd`, and remote object storage access is signed with AWS SigV4.

## Capabilities

| Area | Modules | What it does |
| --- | --- | --- |
| Containers | `image`, `source` | `Image` type for 2D and 3D samples; `ByteSource` trait with in-memory, local-file and HTTP sources |
| Filtering | `filter`, `spatial`, `conv3d`, `edges` | Separable 2D convolution and Gaussian blur, median/mean/Gaussian filters, background subtraction, 3D convolution, Sobel, Scharr and Canny edges |
| Thresholding | `threshold` | Global and tile-based Otsu |
| Segmentation | `label`, `watershed`, `seeds`, `graphcut`, `distance` | 2D and 3D connected components, relabelling, marker-controlled watershed, LoG and local-maxima seeds, graph cut, exact Euclidean distance transform |
| Morphology | `morph` | Binary and grayscale morphology and top-hat in 2D and 3D |
| Objects | `objects`, `features`, `shape3d`, `zernike`, `texture`, `spatial_features` | Filtering by size and border contact, per-object shape and intensity features, 3D shape descriptors, Zernike moments, texture and spatial features |
| Evaluation | `metrics` | Object-level IoU matching with precision, recall and F1 |
| 3D rendering | `render3d` | Camera, transfer functions, ray casting, label overlay, marching-tetrahedra isosurfaces, level-of-detail and mesh rasterising |
| Formats | `io` | PNG, TIFF, BigTIFF, OME-TIFF, memory-mapped TIFF pages, Zarr v2 and v3, OME-NGFF plate metadata, ndarray conversion |
| Caching | `cache` | `BlockCache` wraps a `ByteSource` and keeps aligned fixed-size blocks, aimed at remote reads |
| Remote access | `http`, `sigv4` | HTTP byte-range source and AWS Signature Version 4 request signing |

## Installation

```toml
[dependencies]
tpt-bio-vision = "0.1"
```

## Example

The crate-level documentation describes the module layout. Run the API docs
locally to see the full signatures:

```sh
cargo doc -p tpt-bio-vision --open
```

## Dependencies on the workspace

- `tpt-bio-core` for shared types and statistics.
- `tpt-bio-seq` for sequence-side data used in feature extraction.
- `tpt-crypto-hash` for SHA-256 and HMAC used by SigV4.
- `tpt-simd` for vectorised filter kernels.

External formats are handled by `image`, `tiff`, `ndarray`, `quick-xml`,
`serde_json` and `flate2`.

## Testing

```sh
cargo test -p tpt-bio-vision
```

## Minimum supported Rust version

Rust 1.95, set in the workspace manifest.

## Status

Version 0.1.0, pre-release. The capabilities in the table above are listed as implemented in the crate documentation. Open items are tracked in `todo.md` in the repository root.

## Licence

Dual-licensed under MIT OR Apache-2.0, at your option.
