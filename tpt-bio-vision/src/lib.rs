//! Biological image analysis for tpt-bio.
//!
//! Implemented so far:
//! - [`filter`]: separable 2D convolution and Gaussian blur on `tpt-simd`.
//! - [`sigv4`]: AWS Signature Version 4 request signing for remote object storage,
//!   using `tpt-crypto`'s SHA-256 and HMAC.
//! - [`source`]: the [`source::ByteSource`] trait, with in-memory and local-file sources.
//! - [`image`]: the [`image::Image`] container for 2D and 3D samples.
//! - [`threshold`]: global and tile-based Otsu thresholding.
//! - [`label`]: 2D and 3D connected-component labelling, and relabelling.
//! - [`distance`]: exact Euclidean distance transform.
//! - [`watershed`]: marker-controlled watershed.
//! - [`seeds`]: local maxima and Laplacian-of-Gaussian seed detection.
//! - [`morph`]: binary and grayscale morphology, and top-hat, in 2D and 3D.
//! - [`spatial`]: median, mean and Gaussian filters, and background subtraction.
//! - [`edges`]: Sobel, Scharr and Canny edge detection.
//! - [`objects`]: object filtering by size and border contact.
//! - [`features`]: per-object shape and intensity features.
//! - [`metrics`]: object-level `IoU` matching with precision, recall and `F1`.
//! - [`render3d`]: CPU reference renderer for 3D volumes: camera, transfer functions,
//!   ray casting, label overlay, marching-tetrahedra isosurfaces, LOD and mesh rasterising.
//!
//! - [`io`]: PNG, TIFF, `BigTIFF`, OME-TIFF, memory-mapped TIFF pages, Zarr v2 and v3,
//!   OME-NGFF plate metadata, and `ndarray` conversion.

#![forbid(unsafe_code)]
#![warn(missing_docs)]
// Numeric code indexes pixels by x, y and z; single-letter names are intentional.
#![allow(clippy::many_single_char_names)]

pub mod cache;
pub mod conv3d;
pub mod distance;
pub mod edges;
pub mod features;
pub mod filter;
pub mod graphcut;
pub mod http;
pub mod image;
pub mod io;
pub mod label;
pub mod metrics;
pub mod morph;
pub mod objects;
pub mod render3d;
pub mod seeds;
pub mod shape3d;
pub mod sigv4;
pub mod source;
pub mod spatial;
pub mod spatial_features;
pub mod texture;
pub mod threshold;
pub mod watershed;
pub mod zernike;
