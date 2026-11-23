//! Biological image analysis for tpt-bio.
//!
//! Implemented so far:
//! - [`filter`]: separable 2D convolution and Gaussian blur on `tpt-simd`.
//! - [`sigv4`]: AWS Signature Version 4 request signing for remote object storage,
//!   using `tpt-crypto`'s SHA-256 and HMAC.
//! - [`source`]: the [`source::ByteSource`] trait, with in-memory and local-file sources.
//!
//! Image formats, segmentation and feature extraction are not implemented yet.

#![forbid(unsafe_code)]
#![warn(missing_docs)]

pub mod filter;
pub mod sigv4;
pub mod source;
