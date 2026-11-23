//! CRAM block framing: compression method, content type, sizes, payload and CRC32.
//!
//! Block layout: method (u8), content type (u8), content id (ITF8), compressed
//! size (ITF8), raw size (ITF8), payload, then CRC32 (u32 LE) over every byte of
//! the block before the checksum.
//!
//! Methods: 0 raw, 1 gzip, 4 rANS 4x8. Methods 2 (bzip2) and 3 (LZMA) are
//! reported as [`SeqError::Unsupported`]: no audited pure-Rust, permissively
//! licensed implementation is in the workspace.

use std::io::{Read, Write};

use flate2::{Compression, read::GzDecoder, write::GzEncoder};

use crate::cram::bytes::{Cursor, crc32, put_itf8};
use crate::cram::rans;
use crate::error::{Result, SeqError};

const FORMAT: &str = "CRAM block";

/// Block compression method codes.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Method {
    /// Stored without compression.
    Raw,
    /// gzip (RFC 1952).
    Gzip,
    /// bzip2. Not supported.
    Bzip2,
    /// LZMA. Not supported.
    Lzma,
    /// rANS 4x8.
    Rans4x8,
}

impl Method {
    fn from_code(c: u8) -> Result<Self> {
        Ok(match c {
            0 => Self::Raw,
            1 => Self::Gzip,
            2 => Self::Bzip2,
            3 => Self::Lzma,
            4 => Self::Rans4x8,
            _ => return Err(SeqError::malformed(FORMAT, format!("unknown method {c}"))),
        })
    }

    fn code(self) -> u8 {
        match self {
            Self::Raw => 0,
            Self::Gzip => 1,
            Self::Bzip2 => 2,
            Self::Lzma => 3,
            Self::Rans4x8 => 4,
        }
    }
}

/// Content type codes.
pub mod content {
    /// File header block.
    pub const FILE_HEADER: u8 = 0;
    /// Compression header block.
    pub const COMPRESSION_HEADER: u8 = 1;
    /// Slice header block.
    pub const SLICE_HEADER: u8 = 2;
    /// External data block (one per external content id).
    pub const EXTERNAL_DATA: u8 = 4;
    /// Core data block (bit-coded record fields).
    pub const CORE_DATA: u8 = 5;
}

/// A decoded block.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Block {
    /// Content type code (see [`content`]).
    pub content_type: u8,
    /// Content id (zero for core data).
    pub content_id: i32,
    /// Uncompressed payload.
    pub data: Vec<u8>,
}

/// Decompresses a payload with the given method.
///
/// # Errors
/// Returns [`SeqError::Unsupported`] for bzip2 and LZMA, and
/// [`SeqError::Malformed`] on corrupt data or a size mismatch.
pub fn decompress(method: Method, payload: &[u8], raw_size: usize) -> Result<Vec<u8>> {
    let out = match method {
        Method::Raw => payload.to_vec(),
        Method::Gzip => {
            let mut out = Vec::with_capacity(raw_size);
            GzDecoder::new(payload)
                .read_to_end(&mut out)
                .map_err(|e| SeqError::malformed(FORMAT, format!("gzip: {e}")))?;
            out
        }
        Method::Rans4x8 => rans::decode(payload)?,
        Method::Bzip2 => {
            return Err(SeqError::Unsupported("CRAM block method bzip2 (method 2)".into()));
        }
        Method::Lzma => {
            return Err(SeqError::Unsupported("CRAM block method LZMA (method 3)".into()));
        }
    };
    if out.len() != raw_size {
        return Err(SeqError::malformed(
            FORMAT,
            format!("decoded {} bytes, block declares {raw_size}", out.len()),
        ));
    }
    Ok(out)
}

/// Parses one block from the cursor, verifying its CRC32.
///
/// # Errors
/// Returns [`SeqError::Malformed`] on truncation, checksum mismatch or corrupt data.
pub fn read_block(c: &mut Cursor<'_>) -> Result<Block> {
    let start = c.position();
    let method = Method::from_code(c.u8()?)?;
    let content_type = c.u8()?;
    let content_id = c.itf8()?;
    let comp_len = usize::try_from(c.itf8()?)
        .map_err(|_| SeqError::malformed(FORMAT, "negative compressed size"))?;
    let raw_len =
        usize::try_from(c.itf8()?).map_err(|_| SeqError::malformed(FORMAT, "negative raw size"))?;
    let payload = c.bytes(comp_len)?;
    let end = c.position();
    let stored = c.u32_le()?;
    // The checksum covers the header and payload bytes just consumed.
    let computed = crc32(c.slice_between(start, end)?);
    if stored != computed {
        return Err(SeqError::malformed(FORMAT, "block CRC32 mismatch"));
    }
    let data = decompress(method, payload, raw_len)?;
    Ok(Block { content_type, content_id, data })
}

/// Encodes a block with the given method. The stored CRC32 covers the header and payload.
///
/// # Errors
/// Returns an error only if the gzip encoder fails, or for unsupported methods.
pub fn write_block(
    method: Method,
    content_type: u8,
    content_id: i32,
    raw: &[u8],
    out: &mut Vec<u8>,
) -> Result<()> {
    let payload = match method {
        Method::Raw => raw.to_vec(),
        Method::Gzip => {
            let mut enc = GzEncoder::new(Vec::new(), Compression::default());
            enc.write_all(raw)?;
            enc.finish()?
        }
        Method::Rans4x8 => rans::encode_o0(raw),
        Method::Bzip2 | Method::Lzma => {
            return Err(SeqError::Unsupported("writing bzip2 or LZMA blocks".into()));
        }
    };
    let start = out.len();
    out.push(method.code());
    out.push(content_type);
    put_itf8(out, content_id);
    put_itf8(out, payload.len() as i32);
    put_itf8(out, raw.len() as i32);
    out.extend_from_slice(&payload);
    let sum = crc32(&out[start..]);
    out.extend_from_slice(&sum.to_le_bytes());
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn raw_and_gzip_roundtrip() {
        for method in [Method::Raw, Method::Gzip, Method::Rans4x8] {
            let mut buf = Vec::new();
            write_block(method, content::CORE_DATA, 0, b"hello cram block", &mut buf).unwrap();
            let mut c = Cursor::new(&buf);
            let b = read_block(&mut c).unwrap();
            assert_eq!(b.data, b"hello cram block");
            assert_eq!(b.content_type, content::CORE_DATA);
            assert!(c.is_empty());
        }
    }

    #[test]
    fn corrupt_payload_fails_crc() {
        let mut buf = Vec::new();
        write_block(Method::Raw, content::EXTERNAL_DATA, 7, b"abcdef", &mut buf).unwrap();
        let n = buf.len();
        buf[n - 6] ^= 0x01;
        let err = read_block(&mut Cursor::new(&buf)).unwrap_err();
        assert!(matches!(err, SeqError::Malformed { .. }));
    }

    #[test]
    fn bzip2_block_is_unsupported() {
        let mut buf = Vec::new();
        write_block(Method::Raw, 0, 0, b"x", &mut buf).unwrap();
        buf[0] = 2;
        // Fix the CRC so the failure comes from the method, not the checksum.
        let n = buf.len();
        let sum = crc32(&buf[..n - 4]);
        buf[n - 4..].copy_from_slice(&sum.to_le_bytes());
        let err = read_block(&mut Cursor::new(&buf)).unwrap_err();
        assert!(matches!(err, SeqError::Unsupported(_)));
    }
}
