//! Tabix (TBI) and CSI region indexes for BGZF text (VCF, BED, generic) and BCF.
//!
//! Both formats map each reference to binned chunks of BGZF virtual offsets. A TBI
//! index also keeps a 16 KiB linear index, and its bins use the UCSC scheme shared
//! with BAI ([`crate::bai`]). A CSI index has configurable `min_shift` and `depth`,
//! no linear index, and a per-bin `loffset`. Index files are BGZF-compressed, and
//! [`TabixIndex::read_from`] also accepts an uncompressed index.
//!
//! Query pruning uses only the linear index. The `loffset` of CSI bins is stored
//! and read back, but it is not used to skip chunks.
//!
//! The pseudo-bin 37450 (for TBI, and `depth + 1` levels for CSI) holds each
//! reference's span and mapped/unmapped counts.
//!
//! Coordinates are 0-based and half-open internally. Public query functions take
//! 1-based inclusive regions, as in [`crate::bai::fetch`]'s callers expect for text.

use std::collections::{BTreeMap, HashMap};
use std::io::{Read, Seek, Write};

use crate::bai::Chunk;
use crate::bgzf::{BgzfReader, BgzfWriter, VirtualOffset};
use crate::error::{Result, SeqError};

const FORMAT: &str = "tabix";
const TBI_MAGIC: &[u8; 4] = b"TBI\x01";
const CSI_MAGIC: &[u8; 4] = b"CSI\x01";
const LINEAR_SHIFT: u32 = 14;
const TBI_MIN_SHIFT: u32 = 14;
const TBI_DEPTH: u32 = 5;
/// Format flag meaning zero-based, half-open coordinates (UCSC-style BED).
const ZERO_BASED_FLAG: u32 = 0x1_0000;

/// The bin geometry and flavour of an index.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum IndexKind {
    /// A TBI index: 14-bit leaf bins, five levels, and a linear index.
    Tbi,
    /// A CSI index with leaf bins of `2^min_shift` bases and `depth` levels below the root.
    Csi {
        /// Log2 of the leaf bin size.
        min_shift: u32,
        /// Number of levels below the root.
        depth: u32,
    },
}

impl IndexKind {
    fn geometry(self) -> (u32, u32) {
        match self {
            Self::Tbi => (TBI_MIN_SHIFT, TBI_DEPTH),
            Self::Csi { min_shift, depth } => (min_shift, depth),
        }
    }

    fn magic(self) -> &'static [u8; 4] {
        match self {
            Self::Tbi => TBI_MAGIC,
            Self::Csi { .. } => CSI_MAGIC,
        }
    }

    /// The pseudo-bin that holds reference span and read counts.
    fn meta_bin(self) -> u32 {
        let (_, depth) = self.geometry();
        (((1u64 << (3 * (depth + 1))) - 1) / 7 + 1) as u32
    }

    /// One past the largest coordinate the index can hold.
    fn max_position(self) -> u64 {
        let (min_shift, depth) = self.geometry();
        1u64 << (min_shift + 3 * depth)
    }

    fn validate(self) -> Result<()> {
        if let Self::Csi { min_shift, depth } = self
            && (min_shift == 0 || depth > 10 || min_shift + 3 * depth > 60)
        {
            return Err(SeqError::Unsupported(format!(
                "CSI geometry min_shift={min_shift} depth={depth}"
            )));
        }
        Ok(())
    }
}

/// The text layout of a tabix index (its `format` field and column numbers).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct TabixConfig {
    /// Column preset.
    pub preset: Preset,
    /// Whether begin positions are 0-based (BED-style) rather than 1-based.
    pub zero_based: bool,
    /// 1-based column of the sequence name.
    pub seq_col: u32,
    /// 1-based column of the begin position.
    pub beg_col: u32,
    /// 1-based column of the end position, or 0 for none.
    pub end_col: u32,
    /// Lines starting with this byte are comments. 0 means none.
    pub meta: u8,
    /// Number of leading lines to skip.
    pub skip: u32,
}

/// Column preset stored in the index header.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Preset {
    /// Generic columns.
    Generic,
    /// SAM columns.
    Sam,
    /// VCF columns: the end comes from `END=` or the length of `REF`.
    Vcf,
}

impl TabixConfig {
    /// The VCF preset: sequence in column 1, 1-based `POS` in column 2, `#` comments.
    #[must_use]
    pub const fn vcf() -> Self {
        Self {
            preset: Preset::Vcf,
            zero_based: false,
            seq_col: 1,
            beg_col: 2,
            end_col: 0,
            meta: b'#',
            skip: 0,
        }
    }

    /// The BED preset: 0-based start in column 2 and end in column 3.
    #[must_use]
    pub const fn bed() -> Self {
        Self {
            preset: Preset::Generic,
            zero_based: true,
            seq_col: 1,
            beg_col: 2,
            end_col: 3,
            meta: b'#',
            skip: 0,
        }
    }

    fn format_bits(self) -> u32 {
        let base = match self.preset {
            Preset::Generic => 0,
            Preset::Sam => 1,
            Preset::Vcf => 2,
        };
        base | if self.zero_based { ZERO_BASED_FLAG } else { 0 }
    }

    fn from_fields(
        format: u32,
        seq: u32,
        beg: u32,
        end: u32,
        meta: u32,
        skip: u32,
    ) -> Result<Self> {
        let preset = match format & 0xffff {
            0 => Preset::Generic,
            1 => Preset::Sam,
            2 => Preset::Vcf,
            other => return Err(SeqError::Unsupported(format!("tabix preset {other}"))),
        };
        Ok(Self {
            preset,
            zero_based: format & ZERO_BASED_FLAG != 0,
            seq_col: seq,
            beg_col: beg,
            end_col: end,
            meta: u8::try_from(meta)
                .map_err(|_| SeqError::Unsupported(format!("meta character {meta}")))?,
            skip,
        })
    }

    /// The sequence name, 0-based begin, and exclusive end of a text line.
    fn text_extent<'a>(&self, line: &'a str) -> Result<(&'a str, u64, u64)> {
        let cols: Vec<&str> = line.split('\t').collect();
        let name = column(&cols, self.seq_col)?;
        let raw_beg = parse_u64(column(&cols, self.beg_col)?)?;
        let beg0 = if self.zero_based {
            raw_beg
        } else {
            raw_beg
                .checked_sub(1)
                .ok_or_else(|| SeqError::malformed(FORMAT, "position 0 in 1-based column"))?
        };
        let end0 = if self.preset == Preset::Vcf {
            match info_end(cols.get(7).copied()) {
                Some(e) => e,
                None => beg0 + column(&cols, 4)?.len() as u64,
            }
        } else if self.end_col > 0 {
            parse_u64(column(&cols, self.end_col)?)?
        } else {
            beg0 + 1
        };
        Ok((name, beg0, end0.max(beg0 + 1)))
    }
}

fn column<'a>(cols: &[&'a str], n: u32) -> Result<&'a str> {
    (n as usize)
        .checked_sub(1)
        .and_then(|i| cols.get(i).copied())
        .ok_or_else(|| SeqError::malformed(FORMAT, format!("line has no column {n}")))
}

fn parse_u64(s: &str) -> Result<u64> {
    s.parse().map_err(|_| SeqError::malformed(FORMAT, format!("not a position: {s}")))
}

fn info_end(info: Option<&str>) -> Option<u64> {
    info?.split(';').find_map(|kv| kv.strip_prefix("END=")?.parse().ok())
}

/// The bin for `[beg, end)` (0-based half-open) under the given geometry.
fn reg2bin(beg: u64, end: u64, min_shift: u32, depth: u32) -> u32 {
    let last = end - 1;
    for level in (0..=depth).rev() {
        let shift = min_shift + 3 * (depth - level);
        if beg >> shift == last >> shift {
            return (first_bin(level) + (beg >> shift)) as u32;
        }
    }
    0
}

/// Every bin that may hold a record overlapping `[beg, end)` (0-based half-open).
pub(crate) fn reg2bins(beg: u64, end: u64, min_shift: u32, depth: u32) -> Vec<u32> {
    let last = end - 1;
    let mut bins = Vec::new();
    for level in 0..=depth {
        let shift = min_shift + 3 * (depth - level);
        for b in (beg >> shift)..=(last >> shift) {
            bins.push((first_bin(level) + b) as u32);
        }
    }
    bins
}

/// The first bin id on a level. Level 0 is the root.
fn first_bin(level: u32) -> u64 {
    ((1u64 << (3 * level)) - 1) / 7
}

/// Index data for one bin.
#[derive(Debug, Clone, PartialEq, Eq)]
struct Bin {
    loffset: VirtualOffset,
    chunks: Vec<Chunk>,
}

/// Reference span and read counts from the pseudo-bin.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct RefMeta {
    beg: VirtualOffset,
    end: VirtualOffset,
    n_mapped: u64,
    n_unmapped: u64,
}

#[derive(Debug, Clone, Default, PartialEq, Eq)]
struct RefIndex {
    bins: BTreeMap<u32, Bin>,
    linear: Vec<VirtualOffset>,
    meta: Option<RefMeta>,
}

/// A tabix (TBI) or CSI index.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct TabixIndex {
    kind: IndexKind,
    config: Option<TabixConfig>,
    names: Vec<String>,
    refs: Vec<RefIndex>,
    /// Records with no coordinate, written after the reference data. Zero for an index
    /// built from placed records; read when present.
    n_no_coor: u64,
}

impl TabixIndex {
    /// The index flavour.
    #[must_use]
    pub fn kind(&self) -> IndexKind {
        self.kind
    }

    /// The text layout, or `None` for a BCF index.
    #[must_use]
    pub fn config(&self) -> Option<TabixConfig> {
        self.config
    }

    /// Sequence names in index order. Empty for a BCF index, which uses contig indices.
    #[must_use]
    pub fn names(&self) -> &[String] {
        &self.names
    }

    /// Number of references.
    #[must_use]
    pub fn n_refs(&self) -> usize {
        self.refs.len()
    }

    /// The index of a sequence name.
    #[must_use]
    pub fn ref_id(&self, name: &str) -> Option<usize> {
        self.names.iter().position(|n| n == name)
    }

    /// Chunks that may hold records overlapping a 1-based inclusive region.
    ///
    /// # Errors
    /// Returns [`SeqError::NotFound`] for an unknown sequence and
    /// [`SeqError::Malformed`] for an invalid region.
    pub fn query_chunks(&self, chrom: &str, start: u64, end: u64) -> Result<Vec<Chunk>> {
        check_region(start, end)?;
        let id =
            self.ref_id(chrom).ok_or_else(|| SeqError::NotFound(format!("sequence {chrom}")))?;
        Ok(self.query_by_id(id, start, end))
    }

    /// Chunks that may hold records overlapping a 1-based inclusive region on a
    /// reference index. Returns an empty list for an unknown reference or region.
    #[must_use]
    pub fn query_by_id(&self, ref_id: usize, start: u64, end: u64) -> Vec<Chunk> {
        let Some(r) = self.refs.get(ref_id) else {
            return Vec::new();
        };
        let beg0 = start.saturating_sub(1);
        let end0 = end.min(self.kind.max_position());
        if end0 <= beg0 {
            return Vec::new();
        }
        let min_off = match self.kind {
            IndexKind::Tbi => match r.linear.get((beg0 >> LINEAR_SHIFT) as usize) {
                Some(v) => *v,
                None => return Vec::new(),
            },
            IndexKind::Csi { .. } => VirtualOffset::from_raw(0),
        };
        let (min_shift, depth) = self.kind.geometry();
        let mut hits = Vec::new();
        for bin in reg2bins(beg0, end0, min_shift, depth) {
            if let Some(b) = r.bins.get(&bin) {
                hits.extend(b.chunks.iter().filter(|c| c.end > min_off).copied());
            }
        }
        hits.sort_by_key(|c| c.beg);
        merge_chunks(hits)
    }

    /// Builds an index by scanning BGZF-compressed text, such as a `.vcf.gz` or
    /// `.bed.gz`, from `input`. The input must be sorted by sequence (each name in one
    /// contiguous block) and by position within each sequence.
    ///
    /// # Errors
    /// Returns [`SeqError::Malformed`] for unsorted input or unparseable lines, or an
    /// I/O error.
    pub fn build_text<R: Read>(input: R, config: TabixConfig, kind: IndexKind) -> Result<Self> {
        kind.validate()?;
        let mut bgzf = BgzfReader::new(input);
        let mut builder = IndexBuilder::new(kind, Some(config));
        let mut line = Vec::new();
        let mut line_no = 0u64;
        loop {
            let begin = bgzf.virtual_position();
            line.clear();
            if read_line(&mut bgzf, &mut line)? == 0 {
                break;
            }
            let end = bgzf.virtual_position();
            line_no += 1;
            if line.last() == Some(&b'\r') {
                line.pop();
            }
            if line.is_empty()
                || line_no <= u64::from(config.skip)
                || (config.meta != 0 && line[0] == config.meta)
            {
                continue;
            }
            let text = std::str::from_utf8(&line)
                .map_err(|_| SeqError::malformed(FORMAT, format!("line {line_no} is not UTF-8")))?;
            let (name, beg0, end0) = config
                .text_extent(text)
                .map_err(|e| SeqError::malformed(FORMAT, format!("line {line_no}: {e}")))?;
            let id = builder.resolve_name(name);
            builder.add(id, beg0, end0, begin, end)?;
        }
        Ok(builder.finish())
    }

    /// Writes the index as BGZF-compressed bytes, in TBI or CSI layout.
    ///
    /// # Errors
    /// Returns [`SeqError::Unsupported`] for a TBI index without text configuration, or
    /// an I/O error.
    pub fn write_to<W: Write>(&self, w: W) -> Result<()> {
        let bytes = self.to_bytes()?;
        let mut out = BgzfWriter::new(w);
        out.write_all(&bytes)?;
        out.finish()?;
        Ok(())
    }

    /// Reads a TBI or CSI index. BGZF-compressed and uncompressed input are both
    /// accepted.
    ///
    /// # Errors
    /// Returns [`SeqError::Malformed`] for a bad magic or truncated data, or an I/O
    /// error.
    pub fn read_from<R: Read>(mut r: R) -> Result<Self> {
        let mut raw = Vec::new();
        r.read_to_end(&mut raw)?;
        let data = if raw.starts_with(&[0x1f, 0x8b]) {
            BgzfReader::new(&raw[..]).read_to_vec()?
        } else {
            raw
        };
        parse_index(&data)
    }

    fn to_bytes(&self) -> Result<Vec<u8>> {
        self.kind.validate()?;
        let mut out = Vec::new();
        out.extend_from_slice(self.kind.magic());
        match self.kind {
            IndexKind::Tbi => {
                let cfg = self.config.ok_or_else(|| {
                    SeqError::Unsupported("a TBI index needs a text configuration".into())
                })?;
                put_i32(&mut out, count_i32(self.refs.len())?);
                put_tabix_block(&mut out, cfg, &self.names)?;
            }
            IndexKind::Csi { min_shift, depth } => {
                let aux = match self.config {
                    Some(cfg) => {
                        let mut block = Vec::new();
                        put_tabix_block(&mut block, cfg, &self.names)?;
                        block
                    }
                    None => Vec::new(),
                };
                put_i32(&mut out, count_i32(min_shift as usize)?);
                put_i32(&mut out, count_i32(depth as usize)?);
                put_i32(&mut out, count_i32(aux.len())?);
                out.extend_from_slice(&aux);
                put_i32(&mut out, count_i32(self.refs.len())?);
            }
        }
        let meta_bin = self.kind.meta_bin();
        let csi = matches!(self.kind, IndexKind::Csi { .. });
        for r in &self.refs {
            put_i32(&mut out, count_i32(r.bins.len() + usize::from(r.meta.is_some()))?);
            for (id, bin) in &r.bins {
                put_u32(&mut out, *id);
                if csi {
                    put_u64(&mut out, bin.loffset.raw());
                }
                put_i32(&mut out, count_i32(bin.chunks.len())?);
                for c in &bin.chunks {
                    put_u64(&mut out, c.beg.raw());
                    put_u64(&mut out, c.end.raw());
                }
            }
            if let Some(m) = r.meta {
                put_u32(&mut out, meta_bin);
                if csi {
                    put_u64(&mut out, 0);
                }
                put_i32(&mut out, 2);
                put_u64(&mut out, m.beg.raw());
                put_u64(&mut out, m.end.raw());
                put_u64(&mut out, m.n_mapped);
                put_u64(&mut out, m.n_unmapped);
            }
            if !csi {
                put_i32(&mut out, count_i32(r.linear.len())?);
                for v in &r.linear {
                    put_u64(&mut out, v.raw());
                }
            }
        }
        put_u64(&mut out, self.n_no_coor);
        Ok(out)
    }
}

/// Incrementally builds an index from records added in file order.
#[derive(Debug)]
pub(crate) struct IndexBuilder {
    kind: IndexKind,
    config: Option<TabixConfig>,
    names: Vec<String>,
    name_ids: HashMap<String, usize>,
    refs: Vec<RefIndex>,
    linear: Vec<Vec<Option<VirtualOffset>>>,
    done: Vec<bool>,
    current: Option<usize>,
    last_beg: u64,
}

impl IndexBuilder {
    pub(crate) fn new(kind: IndexKind, config: Option<TabixConfig>) -> Self {
        Self {
            kind,
            config,
            names: Vec::new(),
            name_ids: HashMap::new(),
            refs: Vec::new(),
            linear: Vec::new(),
            done: Vec::new(),
            current: None,
            last_beg: 0,
        }
    }

    /// Ensures at least `n` references exist (used to cover every BCF contig).
    pub(crate) fn ensure_refs(&mut self, n: usize) {
        while self.refs.len() < n {
            self.refs.push(RefIndex::default());
            self.linear.push(Vec::new());
            self.done.push(false);
        }
    }

    /// Maps a sequence name to its reference index, adding it if new.
    pub(crate) fn resolve_name(&mut self, name: &str) -> usize {
        if let Some(&id) = self.name_ids.get(name) {
            return id;
        }
        let id = self.names.len();
        self.names.push(name.to_string());
        self.name_ids.insert(name.to_string(), id);
        self.ensure_refs(id + 1);
        id
    }

    /// Adds one record, given its extent and its virtual offsets before and after it.
    pub(crate) fn add(
        &mut self,
        ref_id: usize,
        beg0: u64,
        end0: u64,
        begin: VirtualOffset,
        end: VirtualOffset,
    ) -> Result<()> {
        let max = self.kind.max_position();
        if beg0 >= max || end0 > max {
            return Err(SeqError::malformed(
                FORMAT,
                format!("record {beg0}..{end0} exceeds the index range {max}"),
            ));
        }
        let end0 = end0.max(beg0 + 1);
        self.ensure_refs(ref_id + 1);
        if self.current != Some(ref_id) {
            if self.done[ref_id] {
                return Err(SeqError::malformed(FORMAT, "input is not sorted by sequence"));
            }
            if let Some(prev) = self.current {
                self.done[prev] = true;
            }
            self.current = Some(ref_id);
            self.last_beg = 0;
        }
        if beg0 < self.last_beg {
            return Err(SeqError::malformed(FORMAT, "input is not sorted by position"));
        }
        self.last_beg = beg0;

        let (min_shift, depth) = self.kind.geometry();
        let bin_id = reg2bin(beg0, end0, min_shift, depth);
        let loffset = match self.kind {
            IndexKind::Tbi => VirtualOffset::from_raw(0),
            IndexKind::Csi { .. } => begin,
        };
        let r = &mut self.refs[ref_id];
        let bin = r.bins.entry(bin_id).or_insert_with(|| Bin { loffset, chunks: Vec::new() });
        match bin.chunks.last_mut() {
            Some(last) if begin <= last.end => {
                if end > last.end {
                    last.end = end;
                }
            }
            _ => bin.chunks.push(Chunk { beg: begin, end }),
        }
        if let Some(m) = r.meta.as_mut() {
            m.end = end;
            m.n_mapped += 1;
        } else {
            r.meta = Some(RefMeta { beg: begin, end, n_mapped: 1, n_unmapped: 0 });
        }
        if self.kind == IndexKind::Tbi {
            let lin = &mut self.linear[ref_id];
            let first = (beg0 >> LINEAR_SHIFT) as usize;
            let last = ((end0 - 1) >> LINEAR_SHIFT) as usize;
            if lin.len() <= last {
                lin.resize(last + 1, None);
            }
            for slot in &mut lin[first..=last] {
                if slot.is_none() {
                    *slot = Some(begin);
                }
            }
        }
        Ok(())
    }

    /// Finishes the index, filling unset linear windows from their predecessors.
    pub(crate) fn finish(self) -> TabixIndex {
        let kind = self.kind;
        let refs = self
            .refs
            .into_iter()
            .zip(self.linear)
            .map(|(mut r, lin)| {
                if kind == IndexKind::Tbi {
                    // As in htslib, an empty window takes the offset of the next window.
                    let mut next = VirtualOffset::from_raw(0);
                    let mut linear: Vec<VirtualOffset> = Vec::with_capacity(lin.len());
                    for v in lin.into_iter().rev() {
                        next = v.unwrap_or(next);
                        linear.push(next);
                    }
                    linear.reverse();
                    r.linear = linear;
                }
                r
            })
            .collect();
        TabixIndex { kind, config: self.config, names: self.names, refs, n_no_coor: 0 }
    }
}

fn merge_chunks(sorted: Vec<Chunk>) -> Vec<Chunk> {
    let mut merged: Vec<Chunk> = Vec::with_capacity(sorted.len());
    for c in sorted {
        match merged.last_mut() {
            Some(last) if c.beg <= last.end => {
                if c.end > last.end {
                    last.end = c.end;
                }
            }
            _ => merged.push(c),
        }
    }
    merged
}

fn check_region(start: u64, end: u64) -> Result<()> {
    if start < 1 || end < start {
        return Err(SeqError::malformed(FORMAT, "region start must be >= 1 and end >= start"));
    }
    Ok(())
}

/// Reads one line without its newline. Returns the bytes consumed, 0 at end of input.
pub(crate) fn read_line<R: Read>(r: &mut R, out: &mut Vec<u8>) -> Result<usize> {
    let mut byte = [0u8; 1];
    let mut n = 0;
    loop {
        match r.read(&mut byte) {
            Ok(0) => break,
            Ok(_) => {
                n += 1;
                if byte[0] == b'\n' {
                    break;
                }
                out.push(byte[0]);
            }
            Err(e) if e.kind() == std::io::ErrorKind::Interrupted => {}
            Err(e) => return Err(SeqError::Io(e)),
        }
    }
    Ok(n)
}

/// Returns the VCF or BED text lines whose extent overlaps a region.
///
/// `reader` must be a BGZF reader over the same file as `index`. `start` and `end`
/// are 1-based and inclusive. Lines come back in file order.
///
/// # Errors
/// Returns [`SeqError::NotFound`] for an unknown sequence, [`SeqError::Unsupported`] if
/// the index has no text configuration, [`SeqError::Malformed`] for a bad region or
/// corrupt line, or an I/O error.
pub fn fetch_text<R: Read + Seek>(
    reader: &mut BgzfReader<R>,
    index: &TabixIndex,
    chrom: &str,
    start: u64,
    end: u64,
) -> Result<Vec<String>> {
    let config = index
        .config
        .ok_or_else(|| SeqError::Unsupported("index has no text configuration".into()))?;
    check_region(start, end)?;
    let (q_beg, q_end) = (start - 1, end);
    let mut out = Vec::new();
    let mut line = Vec::new();
    for chunk in index.query_chunks(chrom, start, end)? {
        reader.seek_virtual(chunk.beg)?;
        while reader.virtual_position() < chunk.end {
            line.clear();
            if read_line(reader, &mut line)? == 0 {
                break;
            }
            if line.last() == Some(&b'\r') {
                line.pop();
            }
            if line.is_empty() || (config.meta != 0 && line[0] == config.meta) {
                continue;
            }
            let text = std::str::from_utf8(&line)
                .map_err(|_| SeqError::malformed(FORMAT, "indexed line is not UTF-8"))?;
            let (name, beg0, end0) = config.text_extent(text)?;
            if name == chrom && beg0 < q_end && end0 > q_beg {
                out.push(text.to_string());
            }
        }
    }
    Ok(out)
}

fn put_i32(out: &mut Vec<u8>, v: i32) {
    out.extend_from_slice(&v.to_le_bytes());
}

fn put_u32(out: &mut Vec<u8>, v: u32) {
    out.extend_from_slice(&v.to_le_bytes());
}

fn put_u64(out: &mut Vec<u8>, v: u64) {
    out.extend_from_slice(&v.to_le_bytes());
}

fn count_i32(n: usize) -> Result<i32> {
    i32::try_from(n).map_err(|_| SeqError::Unsupported("index field exceeds 32 bits".into()))
}

/// Writes the tabix block: format, columns, meta, skip, then NUL-terminated names.
fn put_tabix_block(out: &mut Vec<u8>, cfg: TabixConfig, names: &[String]) -> Result<()> {
    put_i32(out, cfg.format_bits() as i32);
    put_i32(out, count_i32(cfg.seq_col as usize)?);
    put_i32(out, count_i32(cfg.beg_col as usize)?);
    put_i32(out, count_i32(cfg.end_col as usize)?);
    put_i32(out, i32::from(cfg.meta));
    put_i32(out, count_i32(cfg.skip as usize)?);
    let mut block = Vec::new();
    for name in names {
        block.extend_from_slice(name.as_bytes());
        block.push(0);
    }
    put_i32(out, count_i32(block.len())?);
    out.extend_from_slice(&block);
    Ok(())
}

struct Parser<'a> {
    data: &'a [u8],
    pos: usize,
}

impl<'a> Parser<'a> {
    fn take(&mut self, n: usize) -> Result<&'a [u8]> {
        let end = self
            .pos
            .checked_add(n)
            .filter(|&e| e <= self.data.len())
            .ok_or_else(|| SeqError::malformed(FORMAT, "index is truncated"))?;
        let bytes = &self.data[self.pos..end];
        self.pos = end;
        Ok(bytes)
    }

    fn u32(&mut self) -> Result<u32> {
        let b = self.take(4)?;
        Ok(u32::from_le_bytes([b[0], b[1], b[2], b[3]]))
    }

    fn i32(&mut self) -> Result<i32> {
        Ok(self.u32()? as i32)
    }

    fn u64(&mut self) -> Result<u64> {
        let b = self.take(8)?;
        let mut a = [0u8; 8];
        a.copy_from_slice(b);
        Ok(u64::from_le_bytes(a))
    }

    fn count(&mut self) -> Result<usize> {
        usize::try_from(self.i32()?).map_err(|_| SeqError::malformed(FORMAT, "negative count"))
    }

    /// The tabix block: six integers, `l_nm`, and the names.
    fn tabix_block(&mut self) -> Result<(TabixConfig, Vec<String>)> {
        let format = self.u32()?;
        let seq = self.u32()?;
        let beg = self.u32()?;
        let end = self.u32()?;
        let meta = self.u32()?;
        let skip = self.u32()?;
        let l_nm = self.count()?;
        let raw = self.take(l_nm)?;
        let cfg = TabixConfig::from_fields(format, seq, beg, end, meta, skip)?;
        Ok((cfg, split_names(raw)?))
    }
}

fn split_names(raw: &[u8]) -> Result<Vec<String>> {
    let mut names = Vec::new();
    let mut start = 0;
    for (i, &b) in raw.iter().enumerate() {
        if b == 0 {
            let name = std::str::from_utf8(&raw[start..i])
                .map_err(|_| SeqError::malformed(FORMAT, "sequence name is not UTF-8"))?;
            names.push(name.to_string());
            start = i + 1;
        }
    }
    if start != raw.len() {
        return Err(SeqError::malformed(FORMAT, "names block is not NUL-terminated"));
    }
    Ok(names)
}

/// Reads the trailing count of records with no coordinate. Older writers omit it, so
/// its absence is not an error.
fn read_n_no_coor(p: &mut Parser<'_>) -> Result<u64> {
    if p.data.len() - p.pos < 8 {
        return Ok(0);
    }
    p.u64()
}

fn parse_index(data: &[u8]) -> Result<TabixIndex> {
    let mut p = Parser { data, pos: 0 };
    let magic = p.take(4)?;
    if magic == TBI_MAGIC {
        let kind = IndexKind::Tbi;
        let n_ref = p.count()?;
        let (cfg, names) = p.tabix_block()?;
        check_names(&names, n_ref)?;
        let refs = read_refs(&mut p, n_ref, kind)?;
        let n_no_coor = read_n_no_coor(&mut p)?;
        return Ok(TabixIndex { kind, config: Some(cfg), names, refs, n_no_coor });
    }
    if magic == CSI_MAGIC {
        let min_shift = u32::try_from(p.i32()?)
            .map_err(|_| SeqError::malformed(FORMAT, "negative min_shift"))?;
        let depth =
            u32::try_from(p.i32()?).map_err(|_| SeqError::malformed(FORMAT, "negative depth"))?;
        let kind = IndexKind::Csi { min_shift, depth };
        kind.validate()?;
        let l_aux = p.count()?;
        let aux = p.take(l_aux)?;
        let (config, names) = if aux.is_empty() {
            (None, Vec::new())
        } else {
            let (cfg, names) = Parser { data: aux, pos: 0 }.tabix_block()?;
            (Some(cfg), names)
        };
        let n_ref = p.count()?;
        if config.is_some() {
            check_names(&names, n_ref)?;
        }
        let refs = read_refs(&mut p, n_ref, kind)?;
        let n_no_coor = read_n_no_coor(&mut p)?;
        return Ok(TabixIndex { kind, config, names, refs, n_no_coor });
    }
    Err(SeqError::malformed(FORMAT, "missing TBI or CSI magic"))
}

fn check_names(names: &[String], n_ref: usize) -> Result<()> {
    if names.len() != n_ref {
        return Err(SeqError::malformed(FORMAT, "sequence count does not match n_ref"));
    }
    Ok(())
}

fn read_refs(p: &mut Parser<'_>, n_ref: usize, kind: IndexKind) -> Result<Vec<RefIndex>> {
    let meta_bin = kind.meta_bin();
    let csi = matches!(kind, IndexKind::Csi { .. });
    let mut refs = Vec::new();
    for _ in 0..n_ref {
        let n_bin = p.count()?;
        let mut r = RefIndex::default();
        for _ in 0..n_bin {
            let bin = p.u32()?;
            let loffset = if csi { p.u64()? } else { 0 };
            let n_chunk = p.count()?;
            let mut chunks = Vec::new();
            for _ in 0..n_chunk {
                let beg = VirtualOffset::from_raw(p.u64()?);
                let end = VirtualOffset::from_raw(p.u64()?);
                chunks.push(Chunk { beg, end });
            }
            if bin == meta_bin && n_chunk == 2 {
                r.meta = Some(RefMeta {
                    beg: chunks[0].beg,
                    end: chunks[0].end,
                    n_mapped: chunks[1].beg.raw(),
                    n_unmapped: chunks[1].end.raw(),
                });
            } else {
                r.bins.insert(bin, Bin { loffset: VirtualOffset::from_raw(loffset), chunks });
            }
        }
        if !csi {
            let n_intv = p.count()?;
            for _ in 0..n_intv {
                r.linear.push(VirtualOffset::from_raw(p.u64()?));
            }
        }
        refs.push(r);
    }
    Ok(refs)
}

/// Shared helpers for the crate's tests: deterministic synthetic VCF data and a brute
/// force reference for region queries.
#[cfg(test)]
pub(crate) mod test_support {
    use std::fmt::Write as _;

    /// A small deterministic generator (LCG) so tests need no external data.
    pub(crate) struct Lcg(u64);

    impl Lcg {
        pub(crate) fn new(seed: u64) -> Self {
            Self(seed.wrapping_mul(0x9E37_79B9_7F4A_7C15) | 1)
        }

        pub(crate) fn next(&mut self) -> u64 {
            self.0 = self
                .0
                .wrapping_mul(6_364_136_223_846_793_005)
                .wrapping_add(1_442_695_040_888_963_407);
            self.0 >> 33
        }
    }

    /// Two contigs with `n` sorted records each. Records have canonical text (integer
    /// QUAL, `PASS`, `DP` and optional `END`), so BCF decoding reproduces them exactly.
    pub(crate) fn synthetic_vcf(n: usize, seed: u64) -> String {
        let mut rng = Lcg::new(seed);
        let mut out = String::from(
            "##fileformat=VCFv4.3\n\
##contig=<ID=chr1,length=250000000>\n\
##contig=<ID=chr2,length=250000000>\n\
##FILTER=<ID=PASS,Description=\"All filters passed\">\n\
##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Depth\">\n\
##INFO=<ID=END,Number=1,Type=Integer,Description=\"End\">\n\
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n\
##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Read depth\">\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1\ts2\n",
        );
        let genotypes = ["0/0", "0/1", "1|1", "./.", "1/0"];
        for chrom in ["chr1", "chr2"] {
            let mut pos = 1u64;
            for _ in 0..n {
                pos += 1 + rng.next() % 900;
                let ref_len = 1 + rng.next() % 12;
                let reference: String = "ACGT".chars().cycle().take(ref_len as usize).collect();
                let info = if rng.next().is_multiple_of(40) {
                    format!("DP={};END={}", rng.next() % 500, pos + rng.next() % 60_000)
                } else {
                    format!("DP={}", rng.next() % 500)
                };
                let g1 = genotypes[(rng.next() % 5) as usize];
                let g2 = genotypes[(rng.next() % 5) as usize];
                let _ = writeln!(
                    out,
                    "{chrom}\t{pos}\t.\t{reference}\tT\t{}\tPASS\t{info}\tGT:DP\t{g1}:{}\t{g2}:{}",
                    rng.next() % 99,
                    rng.next() % 60,
                    rng.next() % 60
                );
            }
        }
        out
    }

    /// The data lines of a VCF text (no header).
    pub(crate) fn data_lines(text: &str) -> Vec<String> {
        text.lines().filter(|l| !l.starts_with('#')).map(str::to_string).collect()
    }

    /// Lines whose extent overlaps the 1-based inclusive region, in file order.
    pub(crate) fn brute_force(lines: &[String], chrom: &str, start: u64, end: u64) -> Vec<String> {
        lines
            .iter()
            .filter(|l| {
                let cols: Vec<&str> = l.split('\t').collect();
                if cols[0] != chrom {
                    return false;
                }
                let beg0: u64 = cols[1].parse::<u64>().unwrap() - 1;
                let end0 = cols[7]
                    .split(';')
                    .find_map(|kv| kv.strip_prefix("END=").and_then(|v| v.parse::<u64>().ok()))
                    .unwrap_or(beg0 + cols[3].len() as u64);
                beg0 < end && end0.max(beg0 + 1) > start - 1
            })
            .cloned()
            .collect()
    }
}

#[cfg(test)]
mod tests {
    use super::test_support::{Lcg, brute_force, data_lines, synthetic_vcf};
    use super::*;
    use std::io::Cursor;

    fn bgzf_bytes(text: &str) -> Vec<u8> {
        let mut out = Vec::new();
        let mut w = BgzfWriter::new(&mut out);
        w.write_all(text.as_bytes()).unwrap();
        w.finish().unwrap();
        out
    }

    fn check_queries(text: &str, index: &TabixIndex, n_queries: usize, seed: u64) {
        let gz = bgzf_bytes(text);
        let lines = data_lines(text);
        let mut rng = Lcg::new(seed);
        for _ in 0..n_queries {
            let chrom = if rng.next().is_multiple_of(2) { "chr1" } else { "chr2" };
            let start = 1 + rng.next() % 1_200_000;
            let span = match rng.next() % 3 {
                0 => rng.next() % 50,
                1 => rng.next() % 5_000,
                _ => rng.next() % 150_000,
            };
            let end = start + span;
            let expected = brute_force(&lines, chrom, start, end);
            let mut src = BgzfReader::new(Cursor::new(gz.clone()));
            let got = fetch_text(&mut src, index, chrom, start, end).unwrap();
            assert_eq!(got, expected, "region {chrom}:{start}-{end}");
        }
    }

    #[test]
    fn reg2bins_matches_bai_for_tbi_geometry() {
        let mut rng = Lcg::new(3);
        for _ in 0..500 {
            let beg = rng.next() * 997 % (1 << 28);
            let end = beg + 1 + rng.next() % 300_000;
            let mut ours = reg2bins(beg, end, TBI_MIN_SHIFT, TBI_DEPTH);
            let mut theirs = crate::bai::reg2bins(beg, end);
            ours.sort_unstable();
            theirs.sort_unstable();
            assert_eq!(ours, theirs, "{beg}..{end}");
        }
        assert_eq!(reg2bin(0, 1, 14, 5), 4681);
        assert_eq!(IndexKind::Tbi.meta_bin(), 37450);
    }

    #[test]
    fn tbi_build_serialize_and_query() {
        let text = synthetic_vcf(3000, 11);
        let gz = bgzf_bytes(&text);
        let index = TabixIndex::build_text(&gz[..], TabixConfig::vcf(), IndexKind::Tbi).unwrap();
        assert_eq!(index.names(), ["chr1", "chr2"]);

        let mut file = Vec::new();
        index.write_to(&mut file).unwrap();
        let raw = BgzfReader::new(&file[..]).read_to_vec().unwrap();
        assert_eq!(&raw[..4], b"TBI\x01");
        let back = TabixIndex::read_from(&file[..]).unwrap();
        assert_eq!(back, index);

        check_queries(&text, &back, 150, 21);
    }

    #[test]
    fn csi_build_serialize_and_query_with_other_geometry() {
        let text = synthetic_vcf(2000, 5);
        let gz = bgzf_bytes(&text);
        let kind = IndexKind::Csi { min_shift: 12, depth: 4 };
        let index = TabixIndex::build_text(&gz[..], TabixConfig::vcf(), kind).unwrap();
        let mut file = Vec::new();
        index.write_to(&mut file).unwrap();
        let raw = BgzfReader::new(&file[..]).read_to_vec().unwrap();
        assert_eq!(&raw[..4], b"CSI\x01");
        let back = TabixIndex::read_from(&file[..]).unwrap();
        assert_eq!(back, index);
        check_queries(&text, &back, 150, 77);
    }

    #[test]
    fn bed_input_uses_zero_based_starts() {
        let bed = "chr1\t0\t10\nchr1\t10\t20\nchr1\t20\t21\nchr2\t5\t6\n";
        let gz = bgzf_bytes(bed);
        let index = TabixIndex::build_text(&gz[..], TabixConfig::bed(), IndexKind::Tbi).unwrap();
        let mut src = BgzfReader::new(Cursor::new(gz));
        // 1-based 10 is 0-based 9, inside [0,10) only. 1-based 11 is 0-based 10, inside [10,20).
        let got = fetch_text(&mut src, &index, "chr1", 10, 10).unwrap();
        assert_eq!(got, ["chr1\t0\t10"]);
        let got = fetch_text(&mut src, &index, "chr1", 11, 11).unwrap();
        assert_eq!(got, ["chr1\t10\t20"]);
        let got = fetch_text(&mut src, &index, "chr1", 21, 21).unwrap();
        assert_eq!(got, ["chr1\t20\t21"]);
        assert!(index.config().unwrap().zero_based);
    }

    #[test]
    fn rejects_unsorted_input() {
        let text = "chr1\t50\t.\tA\tG\t.\t.\t.\nchr1\t10\t.\tA\tG\t.\t.\t.\n";
        let gz = bgzf_bytes(text);
        let err = TabixIndex::build_text(&gz[..], TabixConfig::vcf(), IndexKind::Tbi).unwrap_err();
        assert!(err.to_string().contains("not sorted"), "{err}");

        let text =
            "chr1\t5\t.\tA\tG\t.\t.\t.\nchr2\t5\t.\tA\tG\t.\t.\t.\nchr1\t9\t.\tA\tG\t.\t.\t.\n";
        let gz = bgzf_bytes(text);
        assert!(TabixIndex::build_text(&gz[..], TabixConfig::vcf(), IndexKind::Tbi).is_err());
    }

    #[test]
    fn unknown_sequence_and_bad_region_are_errors() {
        let gz = bgzf_bytes(&synthetic_vcf(10, 1));
        let index = TabixIndex::build_text(&gz[..], TabixConfig::vcf(), IndexKind::Tbi).unwrap();
        assert!(matches!(index.query_chunks("chrZ", 1, 10), Err(SeqError::NotFound(_))));
        assert!(index.query_chunks("chr1", 0, 10).is_err());
        assert!(index.query_chunks("chr1", 10, 5).is_err());
    }

    #[test]
    fn read_rejects_bad_magic_and_truncation() {
        assert!(TabixIndex::read_from(&b"NOPE\x00\x00"[..]).is_err());
        let gz = bgzf_bytes(&synthetic_vcf(50, 2));
        let index = TabixIndex::build_text(&gz[..], TabixConfig::vcf(), IndexKind::Tbi).unwrap();
        let mut file = Vec::new();
        index.write_to(&mut file).unwrap();
        let raw = BgzfReader::new(&file[..]).read_to_vec().unwrap();
        assert!(TabixIndex::read_from(&raw[..raw.len() - 9]).is_err());
    }

    #[test]
    fn merge_adjacent_query_chunks() {
        let c = |a: u64, b: u64| Chunk {
            beg: VirtualOffset::from_raw(a),
            end: VirtualOffset::from_raw(b),
        };
        assert_eq!(merge_chunks(vec![c(1, 5), c(5, 9), c(12, 14)]), vec![c(1, 9), c(12, 14)]);
    }
}
