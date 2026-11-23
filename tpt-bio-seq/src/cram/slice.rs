//! CRAM slice header: the block that describes one slice of records.
//!
//! Fields, in order: reference id (ITF8), alignment start, alignment span,
//! record count, record counter (LTF8), number of blocks, content id count and
//! ids, embedded reference content id (ITF8, -1 when absent), reference MD5
//! (16 bytes, all zero when absent) and optional tag bytes. The slice header is
//! itself a block of content type 2.

use crate::cram::bytes::{Cursor, put_itf8, put_ltf8};
use crate::error::{Result, SeqError};

const FORMAT: &str = "CRAM slice header";

/// A parsed slice header.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SliceHeader {
    /// Reference id: an index, -1 for unmapped, -2 for multi-reference.
    pub ref_id: i32,
    /// 1-based alignment start, or 0 when none.
    pub start: i32,
    /// Alignment span on the reference.
    pub span: i32,
    /// Number of records in the slice.
    pub n_records: i32,
    /// Global index of the slice's first record.
    pub record_counter: i64,
    /// Number of data blocks in the slice (core block plus external blocks).
    pub n_blocks: i32,
    /// Content ids of the external blocks.
    pub content_ids: Vec<i32>,
    /// Content id of an embedded reference block, if any.
    pub embedded_ref_id: Option<i32>,
    /// MD5 of the reference region, if recorded.
    pub md5: Option<[u8; 16]>,
    /// Optional tag bytes (kept verbatim).
    pub tags: Vec<u8>,
}

/// Parses a slice header block payload.
///
/// # Errors
/// Returns [`SeqError::Malformed`] for truncated data or negative counts.
pub fn parse(data: &[u8]) -> Result<SliceHeader> {
    let mut c = Cursor::new(data);
    let ref_id = c.itf8()?;
    let start = c.itf8()?;
    let span = c.itf8()?;
    let n_records = c.itf8()?;
    let record_counter = c.ltf8()?;
    let n_blocks = c.itf8()?;
    let n_ids = usize::try_from(c.itf8()?)
        .map_err(|_| SeqError::malformed(FORMAT, "negative content id count"))?;
    let mut content_ids = Vec::with_capacity(n_ids.min(4096));
    for _ in 0..n_ids {
        content_ids.push(c.itf8()?);
    }
    let embedded = c.itf8()?;
    let md5_bytes = c.bytes(16)?;
    let md5 = if md5_bytes.iter().all(|&b| b == 0) {
        None
    } else {
        let mut m = [0u8; 16];
        m.copy_from_slice(md5_bytes);
        Some(m)
    };
    let tags = c.bytes(c.remaining())?.to_vec();
    Ok(SliceHeader {
        ref_id,
        start,
        span,
        n_records,
        record_counter,
        n_blocks,
        content_ids,
        embedded_ref_id: (embedded >= 0).then_some(embedded),
        md5,
        tags,
    })
}

/// Serializes a slice header block payload. [`parse`] reads the result back.
#[must_use]
pub fn write(h: &SliceHeader) -> Vec<u8> {
    let mut out = Vec::new();
    put_itf8(&mut out, h.ref_id);
    put_itf8(&mut out, h.start);
    put_itf8(&mut out, h.span);
    put_itf8(&mut out, h.n_records);
    put_ltf8(&mut out, h.record_counter);
    put_itf8(&mut out, h.n_blocks);
    put_itf8(&mut out, h.content_ids.len() as i32);
    for &id in &h.content_ids {
        put_itf8(&mut out, id);
    }
    put_itf8(&mut out, h.embedded_ref_id.unwrap_or(-1));
    out.extend_from_slice(&h.md5.unwrap_or([0; 16]));
    out.extend_from_slice(&h.tags);
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn slice_header_roundtrip() {
        let h = SliceHeader {
            ref_id: 1,
            start: 500,
            span: 120,
            n_records: 7,
            record_counter: 1 << 33,
            n_blocks: 3,
            content_ids: vec![1, 2],
            embedded_ref_id: None,
            md5: Some(*b"0123456789abcdef"),
            tags: Vec::new(),
        };
        assert_eq!(parse(&write(&h)).unwrap(), h);
    }

    #[test]
    fn hand_derived_bytes_and_absent_md5() {
        // ref 0, start 1, span 2, 1 record, counter 0, 1 block, no content ids,
        // embedded -1 (ITF8 0xFF 0xFF 0xFF 0xFF 0x0F), zero MD5.
        let mut bytes = vec![0x00, 0x01, 0x02, 0x01, 0x00, 0x01, 0x00];
        bytes.extend_from_slice(&[0xFF, 0xFF, 0xFF, 0xFF, 0x0F]);
        bytes.extend_from_slice(&[0u8; 16]);
        let h = parse(&bytes).unwrap();
        assert_eq!(h.ref_id, 0);
        assert_eq!(h.span, 2);
        assert_eq!(h.embedded_ref_id, None);
        assert_eq!(h.md5, None);
        assert_eq!(write(&h), bytes);
    }
}
