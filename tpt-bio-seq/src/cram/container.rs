//! CRAM file definition and container header.
//!
//! The file starts with the magic `CRAM`, a major and minor version byte and a
//! 20-byte file id. Each container starts with a header: block length (i32 LE),
//! reference id, alignment start, span, record count, record counter, base count,
//! block count, landmark array, then a CRC32 (u32 LE) over the header bytes.
//! Blocks follow the header.

use crate::cram::bytes::{Cursor, crc32, put_itf8, put_ltf8};
use crate::error::{Result, SeqError};

const FORMAT: &str = "CRAM";
/// File magic.
pub const MAGIC: &[u8; 4] = b"CRAM";
/// Length of the file definition in bytes (magic, version, file id).
pub const FILE_DEF_LEN: usize = 26;

/// Parsed container header.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ContainerHeader {
    /// Total length of the blocks that follow the header, in bytes.
    pub length: i32,
    /// Reference sequence id, -1 for unmapped, -2 for multi-reference.
    pub ref_id: i32,
    /// Alignment start (1-based, or 0 when none).
    pub start: i32,
    /// Alignment span on the reference.
    pub span: i32,
    /// Number of records in the container.
    pub n_records: i32,
    /// Global record counter.
    pub record_counter: i64,
    /// Number of bases in the container.
    pub n_bases: i64,
    /// Number of blocks.
    pub n_blocks: i32,
    /// Byte offsets of the slice headers, relative to the first block.
    pub landmarks: Vec<i32>,
}

/// CRAM format version.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Version {
    /// Major version (3 for CRAM 3.0).
    pub major: u8,
    /// Minor version.
    pub minor: u8,
}

/// Writes the file definition for CRAM 3.0 with the given 20-byte file id.
#[must_use]
pub fn write_file_definition(file_id: &[u8; 20]) -> Vec<u8> {
    let mut out = Vec::with_capacity(FILE_DEF_LEN);
    out.extend_from_slice(MAGIC);
    out.push(3);
    out.push(0);
    out.extend_from_slice(file_id);
    out
}

/// Reads the file definition, returning the version and file id.
///
/// # Errors
/// Returns [`SeqError::Malformed`] for a bad magic number or truncation, and
/// [`SeqError::Unsupported`] for a major version other than 2 or 3.
pub fn read_file_definition(data: &[u8]) -> Result<(Version, [u8; 20])> {
    if data.len() < FILE_DEF_LEN || &data[..4] != MAGIC {
        return Err(SeqError::malformed(FORMAT, "missing CRAM magic"));
    }
    let version = Version { major: data[4], minor: data[5] };
    if !(version.major == 2 || version.major == 3) {
        return Err(SeqError::Unsupported(format!(
            "CRAM version {}.{}",
            version.major, version.minor
        )));
    }
    let mut id = [0u8; 20];
    id.copy_from_slice(&data[6..26]);
    Ok((version, id))
}

/// Encodes a container header, including its CRC32.
#[must_use]
pub fn write_container_header(h: &ContainerHeader) -> Vec<u8> {
    let mut out = Vec::new();
    out.extend_from_slice(&h.length.to_le_bytes());
    put_itf8(&mut out, h.ref_id);
    put_itf8(&mut out, h.start);
    put_itf8(&mut out, h.span);
    put_itf8(&mut out, h.n_records);
    put_ltf8(&mut out, h.record_counter);
    put_ltf8(&mut out, h.n_bases);
    put_itf8(&mut out, h.n_blocks);
    put_itf8(&mut out, h.landmarks.len() as i32);
    for &l in &h.landmarks {
        put_itf8(&mut out, l);
    }
    let sum = crc32(&out);
    out.extend_from_slice(&sum.to_le_bytes());
    out
}

/// Reads a container header from the cursor, verifying its CRC32.
///
/// # Errors
/// Returns [`SeqError::Malformed`] on truncation, a negative landmark count or a
/// checksum mismatch.
pub fn read_container_header(c: &mut Cursor<'_>) -> Result<ContainerHeader> {
    let start = c.position();
    let length = c.i32_le()?;
    let ref_id = c.itf8()?;
    let astart = c.itf8()?;
    let span = c.itf8()?;
    let n_records = c.itf8()?;
    let record_counter = c.ltf8()?;
    let n_bases = c.ltf8()?;
    let n_blocks = c.itf8()?;
    let n_landmarks = c.itf8()?;
    let n_landmarks = usize::try_from(n_landmarks)
        .map_err(|_| SeqError::malformed(FORMAT, "negative landmark count"))?;
    let mut landmarks = Vec::with_capacity(n_landmarks.min(1024));
    for _ in 0..n_landmarks {
        landmarks.push(c.itf8()?);
    }
    let end = c.position();
    let stored = c.u32_le()?;
    if stored != crc32(c.slice_between(start, end)?) {
        return Err(SeqError::malformed_at(
            FORMAT,
            start as u64,
            "container header CRC32 mismatch",
        ));
    }
    Ok(ContainerHeader {
        length,
        ref_id,
        start: astart,
        span,
        n_records,
        record_counter,
        n_bases,
        n_blocks,
        landmarks,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn file_definition_roundtrip() {
        let id = *b"tpt-test-file-id....";
        let def = write_file_definition(&id);
        assert_eq!(def.len(), FILE_DEF_LEN);
        let (v, got) = read_file_definition(&def).unwrap();
        assert_eq!(v, Version { major: 3, minor: 0 });
        assert_eq!(got, id);
        assert!(read_file_definition(b"BAM\x01\x03\x00").is_err());
    }

    #[test]
    fn container_header_roundtrip_and_crc() {
        let h = ContainerHeader {
            length: 1234,
            ref_id: 2,
            start: 100,
            span: 50,
            n_records: 3,
            record_counter: 1 << 40,
            n_bases: 300,
            n_blocks: 5,
            landmarks: vec![0, 99],
        };
        let mut bytes = write_container_header(&h);
        assert_eq!(read_container_header(&mut Cursor::new(&bytes)).unwrap(), h);
        let i = bytes.len() - 6;
        bytes[i] ^= 0x40;
        assert!(read_container_header(&mut Cursor::new(&bytes)).is_err());
    }
}
