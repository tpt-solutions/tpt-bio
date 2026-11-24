//! Pure Rust BGZF (blocked gzip) reading and writing.
//!
//! BGZF is a series of independent gzip members, each at most 64 KiB, with a
//! `BC` extra subfield that records the block size. Readers and writers here use
//! `flate2` with its `miniz_oxide` backend, so no C code is linked.
//!
//! A *virtual offset* packs the compressed offset of a block (upper 48 bits) and
//! an offset inside its decompressed payload (lower 16 bits). It is the position
//! format used by BAI, CSI and TBI indexes.

use std::collections::VecDeque;
use std::io::{self, Read, Seek, SeekFrom, Write};

/// Reads exactly `len` bytes into a new vector.
///
/// Length fields in untrusted input can claim gigabytes. The vector grows only as
/// bytes actually arrive, so a short input allocates little, and a truncated input
/// returns [`io::ErrorKind::UnexpectedEof`].
pub(crate) fn read_bounded<R: Read>(r: &mut R, len: usize) -> io::Result<Vec<u8>> {
    let mut out = Vec::new();
    let got = r.by_ref().take(len as u64).read_to_end(&mut out)?;
    if got != len {
        return Err(io::Error::new(io::ErrorKind::UnexpectedEof, "truncated input"));
    }
    Ok(out)
}

use flate2::{Compression, Crc, read::DeflateDecoder, write::DeflateEncoder};
use rayon::prelude::*;

use crate::error::{Result, SeqError};

/// Largest uncompressed payload stored in one block.
pub const MAX_BLOCK_PAYLOAD: usize = 0xff00;
/// Largest allowed size of a whole BGZF block.
pub const MAX_BLOCK_SIZE: usize = 0x1_0000;
/// The 28-byte empty block that terminates a well-formed BGZF file.
pub const EOF_MARKER: [u8; 28] = [
    0x1f, 0x8b, 0x08, 0x04, 0x00, 0x00, 0x00, 0x00, 0x00, 0xff, 0x06, 0x00, 0x42, 0x43, 0x02, 0x00,
    0x1b, 0x00, 0x03, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
];

const FIXED_HEADER_LEN: usize = 12;
const FOOTER_LEN: usize = 8;

/// A position in a BGZF file: compressed block offset and in-block offset.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct VirtualOffset(u64);

impl VirtualOffset {
    /// Builds an offset. Returns `None` if `coffset` needs more than 48 bits.
    #[must_use]
    pub const fn new(coffset: u64, uoffset: u16) -> Option<Self> {
        if coffset >= 1 << 48 { None } else { Some(Self((coffset << 16) | uoffset as u64)) }
    }

    /// Interprets a raw 64-bit value as stored in an index file.
    #[must_use]
    pub const fn from_raw(raw: u64) -> Self {
        Self(raw)
    }

    /// The raw 64-bit value as stored in an index file.
    #[must_use]
    pub const fn raw(self) -> u64 {
        self.0
    }

    /// Byte offset of the containing block within the compressed file.
    #[must_use]
    pub const fn compressed_offset(self) -> u64 {
        self.0 >> 16
    }

    /// Offset within the decompressed block payload.
    #[must_use]
    pub const fn uncompressed_offset(self) -> u16 {
        (self.0 & 0xffff) as u16
    }
}

/// Location and sizes of one block inside a complete BGZF byte slice.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct BlockSpan {
    /// Offset of the block's first byte.
    pub start: usize,
    /// Total block size in bytes, header and footer included.
    pub size: usize,
}

/// Parses the fixed 12-byte gzip header and the extra field, returning the total
/// block size. `header` must hold at least 12 bytes; `extra` the XLEN bytes.
fn block_size_from_extra(extra: &[u8]) -> Option<usize> {
    let mut i = 0;
    while i + 4 <= extra.len() {
        let slen = u16::from_le_bytes([extra[i + 2], extra[i + 3]]) as usize;
        let body = i + 4;
        if body + slen > extra.len() {
            return None;
        }
        if extra[i] == b'B' && extra[i + 1] == b'C' && slen == 2 {
            let bsize = u16::from_le_bytes([extra[body], extra[body + 1]]);
            return Some(usize::from(bsize) + 1);
        }
        i = body + slen;
    }
    None
}

fn check_fixed_header(h: &[u8]) -> Result<usize> {
    if h[0] != 0x1f || h[1] != 0x8b || h[2] != 8 || h[3] & 0x04 == 0 {
        return Err(SeqError::malformed("BGZF", "block does not start with a gzip FEXTRA header"));
    }
    Ok(u16::from_le_bytes([h[10], h[11]]) as usize)
}

fn inflate_block(cdata: &[u8], crc: u32, isize: u32) -> Result<Vec<u8>> {
    // The footer is untrusted. Reject an impossible size before reserving memory for it.
    if isize as usize > MAX_BLOCK_SIZE {
        return Err(SeqError::malformed(
            "BGZF",
            format!("ISIZE {isize} exceeds the 64 KiB block limit"),
        ));
    }
    let mut out = Vec::with_capacity(isize as usize);
    DeflateDecoder::new(cdata)
        .read_to_end(&mut out)
        .map_err(|e| SeqError::malformed("BGZF", format!("deflate error: {e}")))?;
    if out.len() != isize as usize {
        return Err(SeqError::malformed(
            "BGZF",
            format!("ISIZE mismatch: footer says {isize}, decoded {}", out.len()),
        ));
    }
    let mut hasher = Crc::new();
    hasher.update(&out);
    if hasher.sum() != crc {
        return Err(SeqError::malformed("BGZF", "CRC32 mismatch"));
    }
    Ok(out)
}

/// Streaming BGZF reader. Implements [`Read`].
///
/// Use [`BgzfReader::seek_virtual`] (requires `R: Seek`) for random access.
#[derive(Debug)]
pub struct BgzfReader<R> {
    inner: R,
    buf: Vec<u8>,
    pos: usize,
    block_start: u64,
    next_start: u64,
    last_block_was_eof: bool,
    /// Compressed input read ahead but not yet split into whole blocks.
    raw: Vec<u8>,
    /// File offset of `raw[0]`.
    raw_start: u64,
    /// Blocks inflated ahead of the reader, in file order.
    ready: VecDeque<ReadyBlock>,
    /// Compressed bytes to read on the next refill. Starts small after a seek, so a
    /// short region query does not inflate a large batch, and doubles up to the limit.
    ahead: usize,
}

/// A block that has been inflated ahead of use.
#[derive(Debug)]
struct ReadyBlock {
    /// File offset of the block.
    start: u64,
    /// Compressed size of the block.
    size: u64,
    /// The inflated bytes.
    data: Vec<u8>,
}

/// Compressed bytes read ahead per batch, at most. Blocks in a batch inflate in parallel.
const READ_AHEAD_BYTES: usize = 4 << 20;
/// Compressed bytes read ahead first, after the reader starts or seeks.
const READ_AHEAD_MIN: usize = 1 << 16;

impl<R: Read> BgzfReader<R> {
    /// Wraps a reader positioned at the start of a BGZF stream.
    pub fn new(inner: R) -> Self {
        Self {
            inner,
            buf: Vec::new(),
            pos: 0,
            block_start: 0,
            next_start: 0,
            last_block_was_eof: false,
            raw: Vec::new(),
            raw_start: 0,
            ready: VecDeque::new(),
            ahead: READ_AHEAD_MIN,
        }
    }

    /// Reads ahead a batch of compressed blocks and inflates them in parallel.
    /// Returns `false` when the input has ended and nothing is left.
    fn refill(&mut self) -> Result<bool> {
        let mut eof = false;
        let target = self.ahead;
        self.ahead = (self.ahead * 2).min(READ_AHEAD_BYTES);
        while self.raw.len() < target && !eof {
            let old = self.raw.len();
            self.raw.resize(old + (1 << 16), 0);
            let n = self.inner.read(&mut self.raw[old..])?;
            self.raw.truncate(old + n);
            eof = n == 0;
        }
        let spans = complete_block_spans(&self.raw)?;
        let complete: usize = spans.last().map_or(0, |s| s.start + s.size);
        if complete == 0 {
            if self.raw.is_empty() {
                return Ok(false);
            }
            return Err(SeqError::malformed_at("BGZF", self.raw_start, "truncated block"));
        }
        let raw = &self.raw;
        let inflated: Vec<Result<Vec<u8>>> = spans
            .par_iter()
            .map(|span| inflate_span(&raw[span.start..span.start + span.size]))
            .collect();
        let mut addr = self.raw_start;
        for (span, data) in spans.iter().zip(inflated) {
            self.ready.push_back(ReadyBlock { start: addr, size: span.size as u64, data: data? });
            addr += span.size as u64;
        }
        self.raw.drain(..complete);
        self.raw_start += complete as u64;
        if eof && !self.raw.is_empty() {
            return Err(SeqError::malformed_at("BGZF", self.raw_start, "truncated block"));
        }
        Ok(true)
    }

    /// Reads and decodes the next block. Returns `false` at a clean end of stream.
    fn load_next_block(&mut self) -> Result<bool> {
        if self.ready.is_empty() && !self.refill()? {
            return Ok(false);
        }
        let Some(block) = self.ready.pop_front() else { return Ok(false) };
        self.buf = block.data;
        self.pos = 0;
        self.block_start = block.start;
        self.next_start = block.start + block.size;
        self.last_block_was_eof = self.buf.is_empty();
        Ok(true)
    }

    /// The virtual offset of the next byte to be read.
    #[must_use]
    pub fn virtual_position(&self) -> VirtualOffset {
        if self.pos >= self.buf.len() {
            // Positioned at a block boundary: point at the start of the next block.
            VirtualOffset::new(self.next_start, 0).unwrap_or(VirtualOffset(0))
        } else {
            VirtualOffset::new(self.block_start, self.pos as u16).unwrap_or(VirtualOffset(0))
        }
    }

    /// Whether the stream ended with the standard empty EOF block.
    ///
    /// Only meaningful after reading to the end. A missing marker usually means a
    /// truncated transfer, so callers may want to warn.
    #[must_use]
    pub fn has_eof_marker(&self) -> bool {
        self.last_block_was_eof
    }

    /// Reads all remaining data into a vector.
    pub fn read_to_vec(&mut self) -> Result<Vec<u8>> {
        let mut out = Vec::new();
        self.read_to_end(&mut out)?;
        Ok(out)
    }
}

impl<R: Read + Seek> BgzfReader<R> {
    /// Moves to a virtual offset previously returned by [`BgzfReader::virtual_position`].
    pub fn seek_virtual(&mut self, voffset: VirtualOffset) -> Result<()> {
        let coffset = voffset.compressed_offset();
        self.inner.seek(SeekFrom::Start(coffset))?;
        self.buf.clear();
        self.pos = 0;
        self.next_start = coffset;
        self.raw.clear();
        self.raw_start = coffset;
        self.ready.clear();
        self.ahead = READ_AHEAD_MIN;
        if voffset.uncompressed_offset() == 0 {
            return Ok(());
        }
        if !self.load_next_block()? {
            return Err(SeqError::NotFound(format!("BGZF block at offset {coffset}")));
        }
        let u = usize::from(voffset.uncompressed_offset());
        if u > self.buf.len() {
            return Err(SeqError::malformed_at("BGZF", coffset, "in-block offset past block end"));
        }
        self.pos = u;
        Ok(())
    }
}

impl<R: Read> Read for BgzfReader<R> {
    fn read(&mut self, out: &mut [u8]) -> io::Result<usize> {
        if out.is_empty() {
            return Ok(0);
        }
        while self.pos >= self.buf.len() {
            match self.load_next_block() {
                Ok(true) => {}
                Ok(false) => return Ok(0),
                Err(SeqError::Io(e)) => return Err(e),
                Err(e) => return Err(io::Error::new(io::ErrorKind::InvalidData, e.to_string())),
            }
        }
        let n = out.len().min(self.buf.len() - self.pos);
        out[..n].copy_from_slice(&self.buf[self.pos..self.pos + n]);
        self.pos += n;
        Ok(n)
    }
}

/// Reads up to `buf.len()` bytes, returning how many were read. Unlike
/// `read_exact`, a short read is not an error.
/// Finds the complete blocks at the start of `data`. A partial block at the end is not
/// an error here, since more input may follow.
fn complete_block_spans(data: &[u8]) -> Result<Vec<BlockSpan>> {
    let mut spans = Vec::new();
    let mut pos = 0usize;
    while data.len() - pos >= FIXED_HEADER_LEN {
        let xlen = check_fixed_header(&data[pos..])?;
        let extra_end = pos + FIXED_HEADER_LEN + xlen;
        if extra_end > data.len() {
            break;
        }
        let size = block_size_from_extra(&data[pos + FIXED_HEADER_LEN..extra_end])
            .ok_or_else(|| SeqError::malformed_at("BGZF", pos as u64, "missing BC subfield"))?;
        if size < FIXED_HEADER_LEN + xlen + FOOTER_LEN || size > MAX_BLOCK_SIZE {
            return Err(SeqError::malformed_at("BGZF", pos as u64, "invalid block size"));
        }
        if pos + size > data.len() {
            break;
        }
        spans.push(BlockSpan { start: pos, size });
        pos += size;
    }
    Ok(spans)
}

/// Inflates one whole block, checking its CRC and size.
fn inflate_span(block: &[u8]) -> Result<Vec<u8>> {
    let xlen = check_fixed_header(block)?;
    let rest = &block[FIXED_HEADER_LEN + xlen..];
    let (cdata, footer) = rest.split_at(rest.len() - FOOTER_LEN);
    let crc = u32::from_le_bytes([footer[0], footer[1], footer[2], footer[3]]);
    let isize = u32::from_le_bytes([footer[4], footer[5], footer[6], footer[7]]);
    inflate_block(cdata, crc, isize)
}

/// Compresses one payload into a complete BGZF block.
fn compress_block(payload: &[u8], level: Compression) -> Result<Vec<u8>> {
    debug_assert!(payload.len() <= MAX_BLOCK_PAYLOAD);
    let mut cdata = deflate(payload, level)?;
    if cdata.len() + 26 > MAX_BLOCK_SIZE {
        // Incompressible data can expand past the block limit. Store it instead.
        cdata = deflate(payload, Compression::none())?;
    }
    let total = 18 + cdata.len() + FOOTER_LEN;
    let mut block = Vec::with_capacity(total);
    block.extend_from_slice(&[0x1f, 0x8b, 0x08, 0x04, 0, 0, 0, 0, 0, 0xff, 6, 0]);
    block.extend_from_slice(&[b'B', b'C', 2, 0]);
    block.extend_from_slice(&((total - 1) as u16).to_le_bytes());
    block.extend_from_slice(&cdata);
    let mut hasher = Crc::new();
    hasher.update(payload);
    block.extend_from_slice(&hasher.sum().to_le_bytes());
    block.extend_from_slice(&(payload.len() as u32).to_le_bytes());
    Ok(block)
}

fn deflate(payload: &[u8], level: Compression) -> Result<Vec<u8>> {
    let mut enc = DeflateEncoder::new(Vec::with_capacity(payload.len() / 2), level);
    enc.write_all(payload)?;
    Ok(enc.finish()?)
}

/// Streaming BGZF writer. Implements [`Write`].
///
/// Call [`BgzfWriter::finish`] to flush the last block and write the EOF marker.
/// Dropping the writer without calling `finish` leaves the file without an EOF
/// marker.
#[derive(Debug)]
pub struct BgzfWriter<W: Write> {
    inner: W,
    buf: Vec<u8>,
    level: Compression,
    written: u64,
}

impl<W: Write> BgzfWriter<W> {
    /// Creates a writer with the default compression level.
    pub fn new(inner: W) -> Self {
        Self::with_level(inner, Compression::default())
    }

    /// Creates a writer with an explicit compression level.
    pub fn with_level(inner: W, level: Compression) -> Self {
        Self { inner, buf: Vec::with_capacity(MAX_BLOCK_PAYLOAD), level, written: 0 }
    }

    /// The virtual offset at which the next byte will be written.
    #[must_use]
    pub fn virtual_position(&self) -> VirtualOffset {
        VirtualOffset::new(self.written, self.buf.len() as u16).unwrap_or(VirtualOffset(0))
    }

    fn flush_block(&mut self) -> Result<()> {
        if self.buf.is_empty() {
            return Ok(());
        }
        let block = compress_block(&self.buf, self.level)?;
        self.inner.write_all(&block)?;
        self.written += block.len() as u64;
        self.buf.clear();
        Ok(())
    }

    /// Flushes the final block, writes the EOF marker and returns the inner writer.
    pub fn finish(mut self) -> Result<W> {
        self.flush_block()?;
        self.inner.write_all(&EOF_MARKER)?;
        self.inner.flush()?;
        Ok(self.inner)
    }
}

impl<W: Write> Write for BgzfWriter<W> {
    fn write(&mut self, mut data: &[u8]) -> io::Result<usize> {
        let total = data.len();
        while !data.is_empty() {
            let room = MAX_BLOCK_PAYLOAD - self.buf.len();
            let n = room.min(data.len());
            self.buf.extend_from_slice(&data[..n]);
            data = &data[n..];
            if self.buf.len() == MAX_BLOCK_PAYLOAD {
                self.flush_block().map_err(|e| io::Error::other(e.to_string()))?;
            }
        }
        Ok(total)
    }

    fn flush(&mut self) -> io::Result<()> {
        // Data is only written in whole blocks; `finish` writes the remainder.
        self.inner.flush()
    }
}

/// Splits a complete BGZF byte slice into block spans without inflating them.
///
/// # Errors
/// Returns [`SeqError::Malformed`] if a block header is invalid or the input ends mid-block.
pub fn scan_blocks(data: &[u8]) -> Result<Vec<BlockSpan>> {
    let mut spans = Vec::new();
    let mut pos = 0usize;
    while pos < data.len() {
        if data.len() - pos < FIXED_HEADER_LEN {
            return Err(SeqError::malformed_at("BGZF", pos as u64, "truncated block header"));
        }
        let xlen = check_fixed_header(&data[pos..])?;
        let extra_end = pos + FIXED_HEADER_LEN + xlen;
        if extra_end > data.len() {
            return Err(SeqError::malformed_at("BGZF", pos as u64, "truncated extra field"));
        }
        let size = block_size_from_extra(&data[pos + FIXED_HEADER_LEN..extra_end])
            .ok_or_else(|| SeqError::malformed_at("BGZF", pos as u64, "missing BC subfield"))?;
        if pos + size > data.len() {
            return Err(SeqError::malformed_at("BGZF", pos as u64, "truncated block"));
        }
        spans.push(BlockSpan { start: pos, size });
        pos += size;
    }
    Ok(spans)
}

/// Decompresses a complete BGZF byte slice, inflating blocks in parallel.
///
/// Block boundaries are found from the headers first, then every block is
/// inflated and CRC-checked on the rayon thread pool.
///
/// # Errors
/// Returns [`SeqError::Malformed`] for any corrupt block.
pub fn decompress_parallel(data: &[u8]) -> Result<Vec<u8>> {
    let spans = scan_blocks(data)?;
    let chunks: Vec<Result<Vec<u8>>> = spans
        .par_iter()
        .map(|span| {
            let block = &data[span.start..span.start + span.size];
            let xlen = check_fixed_header(block)?;
            let consumed = FIXED_HEADER_LEN + xlen;
            let rest = &block[consumed..];
            let (cdata, footer) = rest.split_at(rest.len() - FOOTER_LEN);
            let crc = u32::from_le_bytes([footer[0], footer[1], footer[2], footer[3]]);
            let isize = u32::from_le_bytes([footer[4], footer[5], footer[6], footer[7]]);
            inflate_block(cdata, crc, isize)
        })
        .collect();
    let total: usize = chunks.iter().map(|c| c.as_ref().map_or(0, Vec::len)).sum();
    let mut out = Vec::with_capacity(total);
    for chunk in chunks {
        out.extend_from_slice(&chunk?);
    }
    Ok(out)
}

/// Compresses a byte slice into a complete BGZF stream (with EOF marker),
/// compressing blocks in parallel.
///
/// # Errors
/// Returns [`SeqError::Io`] only if the underlying deflate encoder fails.
pub fn compress_parallel(data: &[u8], level: Compression) -> Result<Vec<u8>> {
    let blocks: Vec<Result<Vec<u8>>> =
        data.par_chunks(MAX_BLOCK_PAYLOAD).map(|c| compress_block(c, level)).collect();
    let mut out = Vec::with_capacity(data.len() / 2 + EOF_MARKER.len());
    for block in blocks {
        out.extend_from_slice(&block?);
    }
    out.extend_from_slice(&EOF_MARKER);
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Cursor;

    fn sample(n: usize) -> Vec<u8> {
        (0..n).map(|i| (i.wrapping_mul(31) ^ (i >> 7)) as u8).collect()
    }

    #[test]
    fn virtual_offset_packing() {
        let v = VirtualOffset::new(123_456, 789).unwrap();
        assert_eq!(v.compressed_offset(), 123_456);
        assert_eq!(v.uncompressed_offset(), 789);
        assert_eq!(VirtualOffset::from_raw(v.raw()), v);
        assert_eq!(VirtualOffset::new(1 << 48, 0), None);
    }

    #[test]
    fn writer_reader_roundtrip_multi_block() {
        let data = sample(3 * MAX_BLOCK_PAYLOAD + 1234);
        let mut w = BgzfWriter::new(Vec::new());
        w.write_all(&data).unwrap();
        let file = w.finish().unwrap();
        assert!(file.ends_with(&EOF_MARKER));

        let mut r = BgzfReader::new(Cursor::new(&file));
        let back = r.read_to_vec().unwrap();
        assert_eq!(back, data);
        assert!(r.has_eof_marker());
    }

    #[test]
    fn parallel_roundtrip_matches_streaming() {
        let data = sample(200_000);
        let file = compress_parallel(&data, Compression::default()).unwrap();
        assert_eq!(decompress_parallel(&file).unwrap(), data);

        let mut r = BgzfReader::new(Cursor::new(&file));
        assert_eq!(r.read_to_vec().unwrap(), data);
    }

    #[test]
    fn empty_stream_is_just_eof_marker() {
        let file = BgzfWriter::new(Vec::new()).finish().unwrap();
        assert_eq!(file, EOF_MARKER);
        assert!(decompress_parallel(&file).unwrap().is_empty());
    }

    #[test]
    fn virtual_seek_returns_same_bytes() {
        let data = sample(2 * MAX_BLOCK_PAYLOAD + 500);
        let mut w = BgzfWriter::new(Vec::new());
        let mut marks = Vec::new();
        for chunk in data.chunks(10_000) {
            marks.push((w.virtual_position(), chunk.len()));
            w.write_all(chunk).unwrap();
        }
        let file = w.finish().unwrap();

        let mut r = BgzfReader::new(Cursor::new(&file));
        let (voff, _) = marks[7];
        let mut skip = vec![0u8; 70_000];
        r.read_exact(&mut skip).unwrap();
        let mut expect_from_zero = data[..70_000].to_vec();
        assert_eq!(skip, expect_from_zero.split_off(0));

        let mut seek = BgzfReader::new(Cursor::new(&file));
        seek.seek_virtual(voff).unwrap();
        let mut got = vec![0u8; 10_000];
        seek.read_exact(&mut got).unwrap();
        let start = 7 * 10_000;
        assert_eq!(got, data[start..start + 10_000]);
    }

    #[test]
    fn detects_corruption_and_truncation() {
        let data = sample(10_000);
        let mut file = compress_parallel(&data, Compression::default()).unwrap();

        let mut corrupt = file.clone();
        let mid = corrupt.len() / 2;
        corrupt[mid] ^= 0xff;
        assert!(decompress_parallel(&corrupt).is_err());

        file.truncate(file.len() - 40);
        assert!(decompress_parallel(&file).is_err());
        let mut r = BgzfReader::new(Cursor::new(&file));
        assert!(r.read_to_vec().is_err());
    }

    #[test]
    fn missing_eof_marker_is_reported_not_fatal() {
        let data = sample(100);
        let mut w = BgzfWriter::new(Vec::new());
        w.write_all(&data).unwrap();
        let mut file = w.finish().unwrap();
        file.truncate(file.len() - EOF_MARKER.len());
        let mut r = BgzfReader::new(Cursor::new(&file));
        assert_eq!(r.read_to_vec().unwrap(), data);
        assert!(!r.has_eof_marker());
    }
}
