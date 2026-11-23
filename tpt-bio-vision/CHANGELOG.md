# Changelog

All notable changes to `tpt-bio-vision` are recorded here. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and the crate uses
[Semantic Versioning](https://semver.org/).

## [Unreleased]

### Added

- Image container for 2D and 3D samples, and the `ByteSource` abstraction with
  in-memory, local-file and HTTP sources.
- Separable convolution and Gaussian filters, median/mean/Gaussian filters,
  background subtraction, 3D convolution, Sobel/Scharr/Canny edges.
- Global and tile-based Otsu thresholding.
- 2D and 3D connected-component labelling, watershed, seed detection, graph cut,
  and exact Euclidean distance transform.
- 2D and 3D binary and grayscale morphology.
- Object filtering, per-object features, 3D shape descriptors, Zernike moments,
  texture and spatial features.
- IoU-based object matching with precision, recall and F1.
- I/O for PNG, TIFF, BigTIFF, OME-TIFF, memory-mapped TIFF pages, Zarr v2 and v3,
  OME-NGFF plate metadata, and ndarray conversion.
- `BlockCache` for remote byte sources.
- AWS Signature Version 4 request signing.
- CPU reference 3D renderer: camera, transfer functions, ray casting, label
  overlay, marching-tetrahedra isosurfaces, level of detail and mesh rasterising.

## [0.1.0]

Initial workspace version. Not yet published.
