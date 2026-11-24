//! CRAM 3.0 reading and writing (partial; see the module list below).
//!
//! Implemented: ITF8/LTF8 integers and CRC32 ([`bytes`]), block framing with
//! raw, gzip and rANS methods, compression and slice headers, data-series value
//! decoding, read-feature record decoding with reference reconstruction, the
//! [`reader::CramReader`] and a [`writer::CramWriter`] for the subset described in their docs.

pub mod block;
pub mod bytes;
pub mod codec;
pub mod compression;
pub mod container;
pub mod crai;
pub mod decode;
pub mod huffman;
pub mod md5;
pub mod rans;
pub mod rans_nx16;
pub mod reader;
pub mod reference;
pub mod slice;
pub mod writer;

#[cfg(test)]
mod roundtrip;
