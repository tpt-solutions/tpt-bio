//! Shared biological types, statistics and bio-math for tpt-bio.
//!
//! This crate is the `no_std` + `alloc` leaf of the workspace. It has no I/O and
//! no `unsafe` code. Enable the default `std` feature for the rayon-based helpers
//! in [`par`].
//!
//! ```
//! use tpt_bio_core::nucleotide::reverse_complement;
//!
//! assert_eq!(reverse_complement(b"ACGTN").unwrap(), b"NACGT".to_vec());
//! ```
#![no_std]
#![forbid(unsafe_code)]
#![warn(missing_docs)]

extern crate alloc;

#[cfg(feature = "std")]
extern crate std;

pub mod amino;
pub mod coords;
pub mod error;
pub mod kmer;
pub mod nucleotide;
pub mod quality;
pub mod stats;

#[cfg(feature = "std")]
pub mod par;

pub use amino::{AminoAcid, CodonTable};
pub use coords::{Interval, Region};
pub use error::{InvalidSymbol, RegionParseError};
pub use nucleotide::{Base, Iupac};
pub use quality::Phred;
