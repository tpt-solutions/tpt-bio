//! CRAM block framing: compression method, content type, sizes, payload and CRC32.
//!
//! Block layout: method (u8), content type (u8), content id (ITF8), compressed
//! size (ITF8), raw size (ITF8), payload, then CRC32 (u32 LE) over every byte of
//! the block before the checksum.
//!
//! Methods: 0 raw, 1 gzip, 4 rANS 4x8 (orders 0 and 1), 5 rANS Nx16, 6 arithmetic,
//! 8 TOK3 (read names, decode only). Method 2 (bzip2) is decoded by the pure-Rust
//! decoder in `bzip2.rs`. Method 3 (LZMA) is reported as
//! [`SeqError::Unsupported`]: no audited pure-Rust, permissively licensed
//! implementation is in the workspace.

use std::io::{Read, Write};

use flate2::{Compression, read::GzDecoder, write::GzEncoder};

use crate::cram::bytes::{Cursor, crc32, put_itf8};
use crate::cram::{arith, bzip2, rans, rans_nx16, tok3};
use crate::error::{Result, SeqError};

const FORMAT: &str = "CRAM block";

/// Block compression method codes.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Method {
    /// Stored without compression.
    Raw,
    /// gzip (RFC 1952).
    Gzip,
    /// bzip2. Decoded by `bzip2.rs`; not written.
    Bzip2,
    /// LZMA. Not supported.
    Lzma,
    /// rANS 4x8 (orders 0 and 1).
    Rans4x8,
    /// rANS Nx16 (CRAM 3.1).
    RansNx16,
    /// TOK3 tokenised read names (CRAM 3.1).
    Tok3,
    /// Adaptive arithmetic coder (CRAM 3.1). The stream's flag byte selects the order.
    Arith,
}

impl Method {
    fn from_code(c: u8) -> Result<Self> {
        Ok(match c {
            0 => Self::Raw,
            1 => Self::Gzip,
            2 => Self::Bzip2,
            3 => Self::Lzma,
            4 => Self::Rans4x8,
            5 => Self::RansNx16,
            6 => Self::Arith,
            8 => Self::Tok3,
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
            Self::RansNx16 => 5,
            Self::Arith => 6,
            Self::Tok3 => 8,
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
/// Returns [`SeqError::Unsupported`] for LZMA, and
/// [`SeqError::Malformed`] on corrupt data or a size mismatch.
pub fn decompress(method: Method, payload: &[u8], raw_size: usize) -> Result<Vec<u8>> {
    // The size is untrusted: refuse impossible sizes before any buffer is reserved.
    if raw_size > super::rans::MAX_OUTPUT {
        return Err(SeqError::malformed("CRAM", format!("block size {raw_size} is too large")));
    }
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
        Method::RansNx16 => rans_nx16::decode(payload, raw_size)?,
        Method::Tok3 => tok3::decode(payload, raw_size)?,
        Method::Arith => arith::decode(payload, raw_size)?,
        Method::Bzip2 => bzip2::decompress(payload, raw_size)?,
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
        Method::RansNx16 => rans_nx16::encode(raw),
        Method::Bzip2 | Method::Lzma | Method::Tok3 | Method::Arith => {
            return Err(SeqError::Unsupported("writing bzip2, LZMA or TOK3 blocks".into()));
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

    /// Frames `payload` as a block with the given method code, for tests that
    /// need a payload the writer does not produce.
    fn frame_block(code: u8, payload: &[u8], raw_len: usize) -> Vec<u8> {
        let mut buf = vec![code, content::CORE_DATA];
        put_itf8(&mut buf, 0);
        put_itf8(&mut buf, payload.len() as i32);
        put_itf8(&mut buf, raw_len as i32);
        buf.extend_from_slice(payload);
        let sum = crc32(&buf);
        buf.extend_from_slice(&sum.to_le_bytes());
        buf
    }

    #[test]
    fn rans_order1_and_nx16_blocks_decode() {
        let data = b"ACGTTGCAACGGTTAACCGGTTAAACCCGGGTTT";
        let o1 = frame_block(4, &rans::encode_o1(data), data.len());
        let block = read_block(&mut Cursor::new(&o1)).unwrap();
        assert_eq!(block.data, data);
        let nx = frame_block(5, &rans_nx16::encode(data), data.len());
        assert_eq!(read_block(&mut Cursor::new(&nx)).unwrap().data, data);
        let mut buf = Vec::new();
        write_block(Method::RansNx16, content::CORE_DATA, 0, data, &mut buf).unwrap();
        assert_eq!(read_block(&mut Cursor::new(&buf)).unwrap().data, data);
    }

    #[test]
    fn bzip2_payload_decodes_to_its_declared_size() {
        // one.bz2 is libbzip2's encoding of b"x".
        let payload = include_bytes!("../../tests/fixtures/bzip2/one.bz2");
        assert_eq!(decompress(Method::Bzip2, payload, 1).unwrap(), b"x");
        assert!(matches!(decompress(Method::Bzip2, payload, 2), Err(SeqError::Malformed { .. })));
    }
}
