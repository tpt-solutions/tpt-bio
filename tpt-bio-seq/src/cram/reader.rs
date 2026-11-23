//! Streaming CRAM 3.0 reader.
//!
//! [`CramReader`] reads the file definition and the SAM header container, then
//! decodes one container at a time and yields [`Record`]s in file order. Without
//! a reference, records that need reference bases fail with
//! [`SeqError::Unsupported`]; unmapped records always decode.
//!
//! [`CramReader::build_index`] scans the file and returns the CRAI entries (one per
//! slice) that [`Crai::query`] answers region queries from.

use std::collections::VecDeque;
use std::io::{ErrorKind, Read, Seek, SeekFrom};

use crate::cram::block::{self, Block, content};
use crate::cram::bytes::Cursor;
use crate::cram::compression::{self, CompressionHeader};
use crate::cram::container::{self, ContainerHeader, FILE_DEF_LEN};
use crate::cram::crai::{Crai, CraiEntry};
use crate::cram::decode::{self, Series};
use crate::cram::reference::CramReference;
use crate::cram::slice::{self, SliceHeader};
use crate::error::{Result, SeqError};
use crate::header::Header;
use crate::record::Record;

const FORMAT: &str = "CRAM";
const CHUNK: usize = 16 * 1024;
const MAX_HEADER: usize = 1 << 20;

/// A byte source with a look-ahead buffer and an absolute file offset.
#[derive(Debug)]
struct Stream<R> {
    inner: R,
    buf: Vec<u8>,
    start: usize,
    abs: u64,
    eof: bool,
}

impl<R: Read> Stream<R> {
    fn new(inner: R) -> Self {
        Self { inner, buf: Vec::new(), start: 0, abs: 0, eof: false }
    }

    fn reset(&mut self) {
        self.buf.clear();
        self.start = 0;
        self.abs = 0;
        self.eof = false;
    }

    fn available(&self) -> &[u8] {
        &self.buf[self.start..]
    }

    /// Reads another chunk. Returns `false` at end of input.
    fn fill_more(&mut self) -> Result<bool> {
        if self.eof {
            return Ok(false);
        }
        if self.start > 0 {
            self.buf.drain(..self.start);
            self.start = 0;
        }
        let old = self.buf.len();
        self.buf.resize(old + CHUNK, 0);
        let n = loop {
            match self.inner.read(&mut self.buf[old..]) {
                Ok(n) => break n,
                Err(e) if e.kind() == ErrorKind::Interrupted => {}
                Err(e) => return Err(e.into()),
            }
        };
        self.buf.truncate(old + n);
        if n == 0 {
            self.eof = true;
            return Ok(false);
        }
        Ok(true)
    }

    /// Ensures at least `n` bytes are buffered. Returns `false` if input ends first.
    fn ensure(&mut self, n: usize) -> Result<bool> {
        while self.available().len() < n {
            if !self.fill_more()? {
                return Ok(false);
            }
        }
        Ok(true)
    }

    fn consume(&mut self, n: usize) {
        self.start += n;
        self.abs += n as u64;
    }

    fn take(&mut self, n: usize) -> Result<Vec<u8>> {
        if !self.ensure(n)? {
            return Err(SeqError::malformed_at(FORMAT, self.abs, "unexpected end of file"));
        }
        let v = self.available()[..n].to_vec();
        self.consume(n);
        Ok(v)
    }
}

/// A container read from the stream: its header, its file offset and its blocks
/// with their byte ranges relative to the first block.
#[derive(Debug)]
struct RawContainer {
    offset: u64,
    header: ContainerHeader,
    blocks: Vec<(u64, u64, Block)>,
}

/// Reads the next container, or `None` at end of input.
fn read_container<R: Read>(s: &mut Stream<R>) -> Result<Option<RawContainer>> {
    if !s.ensure(1)? {
        return Ok(None);
    }
    let offset = s.abs;
    let (header, used) = loop {
        let avail = s.available();
        let mut c = Cursor::new(avail);
        match container::read_container_header(&mut c) {
            Ok(h) => break (h, c.position()),
            Err(e) => {
                if avail.len() > MAX_HEADER || !s.fill_more()? {
                    return Err(e);
                }
            }
        }
    };
    s.consume(used);
    let length = usize::try_from(header.length)
        .map_err(|_| SeqError::malformed_at(FORMAT, offset, "negative container length"))?;
    if !s.ensure(length)? {
        return Err(SeqError::malformed_at(FORMAT, offset, "container extends past end of file"));
    }
    let body = s.available()[..length].to_vec();
    s.consume(length);
    let mut c = Cursor::new(&body);
    let mut blocks = Vec::new();
    while !c.is_empty() {
        let start = c.position() as u64;
        let b = block::read_block(&mut c)?;
        blocks.push((start, c.position() as u64, b));
    }
    Ok(Some(RawContainer { offset, header, blocks }))
}

/// One slice of a container: its header, its data blocks and its byte range.
#[derive(Debug)]
struct SliceUnit {
    header: SliceHeader,
    blocks: Vec<Block>,
    rel_start: u64,
    rel_end: u64,
}

/// Splits a container into its compression header and slices.
fn split_container(c: RawContainer) -> Result<(CompressionHeader, Vec<SliceUnit>)> {
    let mut q: VecDeque<(u64, u64, Block)> = c.blocks.into();
    let (_, _, first) = q
        .pop_front()
        .ok_or_else(|| SeqError::malformed_at(FORMAT, c.offset, "container has no blocks"))?;
    if first.content_type != content::COMPRESSION_HEADER {
        return Err(SeqError::malformed_at(
            FORMAT,
            c.offset,
            "first block is not a compression header",
        ));
    }
    let comp = compression::parse(&first.data)?;
    let mut units = Vec::new();
    while let Some((start, mut end, hb)) = q.pop_front() {
        if hb.content_type != content::SLICE_HEADER {
            return Err(SeqError::malformed_at(FORMAT, c.offset, "expected a slice header block"));
        }
        let sh = slice::parse(&hb.data)?;
        let nb = usize::try_from(sh.n_blocks)
            .map_err(|_| SeqError::malformed_at(FORMAT, c.offset, "negative block count"))?;
        let mut blocks = Vec::with_capacity(nb.min(1024));
        for _ in 0..nb {
            let (_, e, b) = q
                .pop_front()
                .ok_or_else(|| SeqError::malformed_at(FORMAT, c.offset, "slice blocks missing"))?;
            end = e;
            blocks.push(b);
        }
        units.push(SliceUnit { header: sh, blocks, rel_start: start, rel_end: end });
    }
    Ok((comp, units))
}

/// Reads CRAM 3.0 (and 2.x layout) files as a stream of records.
#[derive(Debug)]
pub struct CramReader<R: Read> {
    stream: Stream<R>,
    header: Header,
    reference: Option<CramReference>,
    pending: VecDeque<Record>,
    done: bool,
}

impl<R: Read> CramReader<R> {
    /// Opens a CRAM stream without a reference. Records that need reference bases
    /// fail with [`SeqError::Unsupported`].
    ///
    /// # Errors
    /// Returns an error if the file definition or the SAM header container is bad.
    pub fn new(inner: R) -> Result<Self> {
        Self::open(inner, None)
    }

    /// Opens a CRAM stream and decodes mapped records against `reference`. Slice
    /// MD5 values are checked against the reference when present.
    ///
    /// # Errors
    /// As [`CramReader::new`].
    pub fn with_reference(inner: R, reference: CramReference) -> Result<Self> {
        Self::open(inner, Some(reference))
    }

    fn open(inner: R, reference: Option<CramReference>) -> Result<Self> {
        let mut stream = Stream::new(inner);
        let def = stream.take(FILE_DEF_LEN)?;
        container::read_file_definition(&def)?;
        let first = read_container(&mut stream)?
            .ok_or_else(|| SeqError::malformed(FORMAT, "missing file header container"))?;
        let block = first
            .blocks
            .into_iter()
            .map(|(_, _, b)| b)
            .find(|b| b.content_type == content::FILE_HEADER)
            .ok_or_else(|| SeqError::malformed(FORMAT, "no file header block"))?;
        let mut c = Cursor::new(&block.data);
        let text_len = usize::try_from(c.i32_le()?)
            .map_err(|_| SeqError::malformed(FORMAT, "negative header length"))?;
        let text = c.bytes(text_len)?;
        let text = std::str::from_utf8(text)
            .map_err(|_| SeqError::malformed(FORMAT, "SAM header is not UTF-8"))?;
        let header = Header::parse(text)?;
        Ok(Self { stream, header, reference, pending: VecDeque::new(), done: false })
    }

    /// The SAM header from the file.
    #[must_use]
    pub fn header(&self) -> &Header {
        &self.header
    }

    /// Returns the next record, or `None` at end of file.
    ///
    /// # Errors
    /// Returns [`SeqError::Malformed`] for corrupt data, [`SeqError::Unsupported`]
    /// for features this reader does not implement, and [`SeqError::Io`] for read failures.
    pub fn next_record(&mut self) -> Result<Option<Record>> {
        loop {
            if let Some(r) = self.pending.pop_front() {
                return Ok(Some(r));
            }
            if self.done {
                return Ok(None);
            }
            match read_container(&mut self.stream)? {
                None => self.done = true,
                Some(raw) => {
                    let recs = self.decode_container(raw)?;
                    self.pending = recs.into();
                }
            }
        }
    }

    fn decode_container(&self, raw: RawContainer) -> Result<Vec<Record>> {
        let expected = raw.header.n_records;
        let (comp, units) = split_container(raw)?;
        let series = Series::new(&comp)?;
        let mut out = Vec::new();
        for u in units {
            let recs = decode::decode_slice(
                &series,
                &u.header,
                u.blocks,
                self.header.references(),
                self.reference.as_ref(),
            )?;
            out.extend(recs);
        }
        if i64::try_from(out.len()).ok() != Some(i64::from(expected)) {
            return Err(SeqError::malformed(FORMAT, "record count differs from container header"));
        }
        Ok(out)
    }
}

impl<R: Read + Seek> CramReader<R> {
    /// Scans the whole file and returns one CRAI entry per slice.
    ///
    /// The input is rewound and read to its end, so the reader is exhausted
    /// afterwards: build the index on a reader used only for that purpose.
    ///
    /// # Errors
    /// Returns an error for corrupt containers or I/O failures.
    pub fn build_index(&mut self) -> Result<Crai> {
        self.stream.inner.seek(SeekFrom::Start(0))?;
        self.stream.reset();
        let mut entries = Vec::new();
        // Skip the file definition and the header container.
        self.stream.take(FILE_DEF_LEN)?;
        read_container(&mut self.stream)?;
        while let Some(raw) = read_container(&mut self.stream)? {
            let offset = raw.offset;
            let (_, units) = split_container(raw)?;
            for u in units {
                entries.push(CraiEntry {
                    ref_id: u.header.ref_id,
                    start: u.header.start,
                    span: u.header.span,
                    container_offset: offset,
                    slice_offset: u.rel_start,
                    slice_len: u.rel_end - u.rel_start,
                });
            }
        }
        self.pending.clear();
        self.done = true;
        Ok(Crai::from_entries(entries))
    }
}
