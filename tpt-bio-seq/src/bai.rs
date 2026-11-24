//! BAI index: build from a coordinate-sorted BAM, read and write, and query by region.
//!
//! The index maps each reference to UCSC-binned chunks of virtual offsets, plus
//! a linear index with one virtual offset per 16 KiB window. Queries collect the
//! chunks for every bin that overlaps the region, starting no earlier than the
//! linear-index entry, and then filter records by exact overlap.

use std::collections::BTreeMap;
use std::io::{Read, Seek, Write};

use crate::bam::{BamReader, RecordView, reg2bin};
use crate::bgzf::VirtualOffset;
use crate::error::{Result, SeqError};
use crate::record::Record;

const MAGIC: &[u8; 4] = b"BAI\x01";
const LINEAR_SHIFT: u32 = 14;

/// A contiguous run of records in the compressed file.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Chunk {
    /// Start of the first record.
    pub beg: VirtualOffset,
    /// End of the last record.
    pub end: VirtualOffset,
}

/// Index data for one reference.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct RefIndex {
    bins: BTreeMap<u32, Vec<Chunk>>,
    linear: Vec<VirtualOffset>,
    meta: Option<RefMeta>,
    /// The open run of consecutive records in one bin. It becomes a chunk when the
    /// bin changes, as in htslib.
    run: Option<Run>,
}

/// A run of consecutive records that fall in the same bin.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct Run {
    bin: u32,
    beg: VirtualOffset,
    end: VirtualOffset,
}

/// Deepest binning level for BAI, which uses a `min_shift` of 14.
const N_LVLS: u32 = 5;
/// Bins below this are real bins. The pseudo-bin 37450 holds metadata and is excluded.
const N_BINS: u32 = 37_450;
/// A bin whose chunks span less than this many BGZF blocks is merged into its parent
/// (htslib's `HTS_MIN_MARKER_DIST`, in units of 64 KiB blocks).
const MIN_MARKER_DIST: u64 = 0x1_0000;

impl RefIndex {
    fn add_to_run(&mut self, bin: u32, beg: VirtualOffset, end: VirtualOffset) {
        match &mut self.run {
            Some(run) if run.bin == bin => run.end = end,
            _ => {
                self.flush_run();
                self.run = Some(Run { bin, beg, end });
            }
        }
    }

    fn flush_run(&mut self) {
        if let Some(run) = self.run.take() {
            self.bins.entry(run.bin).or_default().push(Chunk { beg: run.beg, end: run.end });
        }
    }

    /// Applies htslib's compression of the binning index: small bins move to their
    /// parent, chunks are sorted, and chunks that touch the same BGZF block are merged.
    fn compress(&mut self) {
        for level in (1..=N_LVLS).rev() {
            let start = (8u32.pow(level) - 1) / 7;
            let keys: Vec<u32> = self.bins.range(start..N_BINS).map(|(&k, _)| k).collect();
            for key in keys {
                let Some(list) = self.bins.get_mut(&key) else { continue };
                if level < N_LVLS {
                    list.sort_by_key(|c| (c.beg.raw(), c.end.raw()));
                }
                let (Some(head), Some(tail)) = (list.first(), list.last()) else { continue };
                let span = (tail.end.raw() >> 16).wrapping_sub(head.beg.raw() >> 16);
                let parent = (key - 1) >> 3;
                if span < MIN_MARKER_DIST
                    && self.bins.contains_key(&parent)
                    && let Some(moved) = self.bins.remove(&key)
                {
                    self.bins.entry(parent).or_default().extend(moved);
                }
            }
        }
        for list in self.bins.values_mut() {
            list.sort_by_key(|c| (c.beg.raw(), c.end.raw()));
            let mut m = 0;
            for l in 1..list.len() {
                let next = list[l];
                if list[m].end.raw() >> 16 >= next.beg.raw() >> 16 {
                    if list[m].end.raw() < next.end.raw() {
                        list[m].end = next.end;
                    }
                } else {
                    m += 1;
                    list[m] = next;
                }
            }
            list.truncate(if list.is_empty() { 0 } else { m + 1 });
        }
    }
}

/// Per-reference summary stored in the pseudo-bin 37450 that samtools writes: the
/// virtual span of the reference's records and their mapped and unmapped counts.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct RefMeta {
    beg: VirtualOffset,
    end: VirtualOffset,
    n_mapped: u64,
    n_unmapped: u64,
}

/// The pseudo-bin that holds [`RefMeta`]. It is the bin after every real bin.
const META_BIN: u32 = 37_450;

/// A BAI index for a whole BAM file.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct Index {
    refs: Vec<RefIndex>,
    n_no_coor: Option<u64>,
}

impl Index {
    /// Number of references covered.
    #[must_use]
    pub fn n_refs(&self) -> usize {
        self.refs.len()
    }

    /// Records with no coordinate (`ref_id == -1`), when the index stores the count.
    #[must_use]
    pub fn unplaced_count(&self) -> Option<u64> {
        self.n_no_coor
    }

    /// Adds one record's extent to the index.
    ///
    /// `begin` and `end` are the virtual offsets before and after the record.
    /// Records must be added in file order.
    pub fn add_record(&mut self, view: &RecordView<'_>, begin: VirtualOffset, end: VirtualOffset) {
        let ref_id = view.ref_id();
        if ref_id < 0 || view.pos() < 0 {
            *self.n_no_coor.get_or_insert(0) += 1;
            return;
        }
        let rid = ref_id as usize;
        if self.refs.len() <= rid {
            self.refs.resize_with(rid + 1, RefIndex::default);
        }
        let beg = view.pos();
        let end_pos = view.cigar().map_or(beg + 1, |c| {
            let span = crate::record::reference_length(&c) as i64;
            if span == 0 { beg + 1 } else { beg + span }
        });
        let bin = reg2bin(beg, end_pos);
        let r = &mut self.refs[rid];
        r.add_to_run(u32::from(bin), begin, end);
        let meta = r.meta.get_or_insert(RefMeta { beg: begin, end, n_mapped: 0, n_unmapped: 0 });
        meta.end = end;
        if view.flags().0 & crate::record::Flags::UNMAPPED != 0 {
            meta.n_unmapped += 1;
        } else {
            meta.n_mapped += 1;
        }

        let first = (beg >> LINEAR_SHIFT) as usize;
        let last = ((end_pos - 1) >> LINEAR_SHIFT) as usize;
        if r.linear.len() <= last {
            r.linear.resize(last + 1, VirtualOffset::from_raw(0));
        }
        for w in first..=last {
            if r.linear[w].raw() == 0 {
                r.linear[w] = begin;
            }
        }
    }

    /// Fills unset linear-index windows from their predecessors. Call once after
    /// all records are added.
    pub fn finish(&mut self) {
        for r in &mut self.refs {
            r.flush_run();
            r.compress();
            // htslib fills an empty window from the next window that has an offset, so
            // walk backwards. The last window is always set.
            for i in (0..r.linear.len().saturating_sub(1)).rev() {
                if r.linear[i].raw() == 0 {
                    r.linear[i] = r.linear[i + 1];
                }
            }
        }
    }

    /// Chunks that may contain records overlapping `[beg, end)` on `ref_id`.
    /// `beg` and `end` are 0-based half-open.
    #[must_use]
    pub fn chunks_for(&self, ref_id: usize, beg: u64, end: u64) -> Vec<Chunk> {
        let Some(r) = self.refs.get(ref_id) else { return Vec::new() };
        if end <= beg {
            return Vec::new();
        }
        let window = (beg >> LINEAR_SHIFT) as usize;
        let min_offset = match r.linear.get(window) {
            Some(v) => *v,
            None => return Vec::new(),
        };
        let mut hits = Vec::new();
        for bin in reg2bins(beg, end) {
            if let Some(chunks) = r.bins.get(&bin) {
                hits.extend(chunks.iter().copied().filter(|c| c.end > min_offset));
            }
        }
        hits.sort_by_key(|c| c.beg);
        let mut merged: Vec<Chunk> = Vec::with_capacity(hits.len());
        for c in hits {
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

    /// Serializes in the BAI binary format.
    ///
    /// # Errors
    /// Returns an error on I/O failure.
    pub fn write_to<W: Write>(&self, mut w: W) -> Result<()> {
        w.write_all(MAGIC)?;
        w.write_all(&(self.refs.len() as i32).to_le_bytes())?;
        for r in &self.refs {
            let n_bins = r.bins.len() + usize::from(r.meta.is_some());
            w.write_all(&(n_bins as i32).to_le_bytes())?;
            for (bin, chunks) in &r.bins {
                w.write_all(&bin.to_le_bytes())?;
                w.write_all(&(chunks.len() as i32).to_le_bytes())?;
                for c in chunks {
                    w.write_all(&c.beg.raw().to_le_bytes())?;
                    w.write_all(&c.end.raw().to_le_bytes())?;
                }
            }
            if let Some(m) = r.meta {
                w.write_all(&META_BIN.to_le_bytes())?;
                w.write_all(&2i32.to_le_bytes())?;
                w.write_all(&m.beg.raw().to_le_bytes())?;
                w.write_all(&m.end.raw().to_le_bytes())?;
                w.write_all(&m.n_mapped.to_le_bytes())?;
                w.write_all(&m.n_unmapped.to_le_bytes())?;
            }
            w.write_all(&(r.linear.len() as i32).to_le_bytes())?;
            for v in &r.linear {
                w.write_all(&v.raw().to_le_bytes())?;
            }
        }
        // htslib always writes the unplaced-read count, even when it is zero.
        w.write_all(&self.n_no_coor.unwrap_or(0).to_le_bytes())?;
        Ok(())
    }

    /// Parses the BAI binary format.
    ///
    /// # Errors
    /// Returns [`SeqError::Malformed`] for a bad magic number or truncated data.
    pub fn read_from<R: Read>(mut r: R) -> Result<Self> {
        let mut magic = [0u8; 4];
        r.read_exact(&mut magic).map_err(|_| SeqError::malformed("BAI", "truncated magic"))?;
        if &magic != MAGIC {
            return Err(SeqError::malformed("BAI", "missing BAI magic"));
        }
        let n_ref = read_i32(&mut r)?;
        let n_ref =
            usize::try_from(n_ref).map_err(|_| SeqError::malformed("BAI", "negative n_ref"))?;
        let mut refs = Vec::with_capacity(n_ref.min(1 << 16));
        for _ in 0..n_ref {
            let n_bin = read_i32(&mut r)?;
            let mut bins = BTreeMap::new();
            for _ in 0..n_bin.max(0) {
                let bin = read_u32(&mut r)?;
                let n_chunk = read_i32(&mut r)?;
                let mut chunks = Vec::new();
                for _ in 0..n_chunk.max(0) {
                    let beg = VirtualOffset::from_raw(read_u64(&mut r)?);
                    let end = VirtualOffset::from_raw(read_u64(&mut r)?);
                    chunks.push(Chunk { beg, end });
                }
                bins.insert(bin, chunks);
            }
            let n_intv = read_i32(&mut r)?;
            let mut linear = Vec::new();
            for _ in 0..n_intv.max(0) {
                linear.push(VirtualOffset::from_raw(read_u64(&mut r)?));
            }
            refs.push(RefIndex { bins, linear, meta: None, run: None });
        }
        // The trailing unplaced-read count is optional in the format.
        let mut tail = [0u8; 8];
        let n_no_coor = match r.read_exact(&mut tail) {
            Ok(()) => Some(u64::from_le_bytes(tail)),
            Err(_) => None,
        };
        Ok(Self { refs, n_no_coor })
    }
}

/// All UCSC bins that may hold records overlapping `[beg, end)`.
#[must_use]
pub fn reg2bins(beg: u64, end: u64) -> Vec<u32> {
    let beg = beg as i64;
    let end = end as i64 - 1;
    let mut bins = vec![0u32];
    for (offset, shift) in [(1i64, 26u32), (9, 23), (73, 20), (585, 17), (4681, 14)] {
        let lo = offset + (beg >> shift);
        let hi = offset + (end >> shift);
        bins.extend((lo..=hi).map(|b| b as u32));
    }
    bins
}

fn read_i32<R: Read>(r: &mut R) -> Result<i32> {
    let mut b = [0u8; 4];
    r.read_exact(&mut b).map_err(|_| SeqError::malformed("BAI", "truncated index"))?;
    Ok(i32::from_le_bytes(b))
}

fn read_u32<R: Read>(r: &mut R) -> Result<u32> {
    let mut b = [0u8; 4];
    r.read_exact(&mut b).map_err(|_| SeqError::malformed("BAI", "truncated index"))?;
    Ok(u32::from_le_bytes(b))
}

fn read_u64<R: Read>(r: &mut R) -> Result<u64> {
    let mut b = [0u8; 8];
    r.read_exact(&mut b).map_err(|_| SeqError::malformed("BAI", "truncated index"))?;
    Ok(u64::from_le_bytes(b))
}

/// Builds an index by scanning a coordinate-sorted BAM from its first record.
///
/// # Errors
/// Returns [`SeqError::Malformed`] if the BAM is corrupt or not coordinate-sorted.
pub fn build_index<R: Read>(reader: &mut BamReader<R>) -> Result<Index> {
    let mut index = Index::default();
    let mut last: Option<(i32, i64)> = None;
    while let Some((begin, end, view)) = reader.next_view_with_offsets()? {
        if view.ref_id() >= 0 {
            let key = (view.ref_id(), view.pos());
            if let Some(prev) = last
                && (key.0 < prev.0 || (key.0 == prev.0 && key.1 < prev.1))
            {
                return Err(SeqError::malformed("BAM", "input is not coordinate-sorted"));
            }
            last = Some(key);
        }
        index.add_record(&view, begin, end);
    }
    index.finish();
    Ok(index)
}

/// Fetches every record overlapping `[beg, end)` (0-based half-open) on `ref_id`.
///
/// # Errors
/// Returns an error on I/O or decoding failure.
pub fn fetch<R: Read + Seek>(
    reader: &mut BamReader<R>,
    index: &Index,
    ref_id: usize,
    beg: u64,
    end: u64,
) -> Result<Vec<Record>> {
    let mut out = Vec::new();
    for chunk in index.chunks_for(ref_id, beg, end) {
        reader.seek_virtual(chunk.beg)?;
        loop {
            if reader.virtual_position() >= chunk.end {
                break;
            }
            let Some(view) = reader.next_view()? else { break };
            if view.ref_id() != ref_id as i32 {
                if view.ref_id() > ref_id as i32 {
                    break;
                }
                continue;
            }
            let rstart = view.pos().max(0) as u64;
            // A record with no CIGAR, or a malformed one, is treated as one base long.
            let rend = rstart + view.reference_span().unwrap_or(0).max(1);
            if rstart < end && rend > beg {
                out.push(view.to_record()?);
            }
        }
    }
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn reg2bins_contains_leaf_bin_of_point() {
        let bins = reg2bins(100_000, 100_001);
        assert!(bins.contains(&u32::from(reg2bin(100_000, 100_001))));
        assert!(bins.contains(&0));
    }

    #[test]
    fn index_serialization_roundtrip() {
        // Files always carry the unplaced-read count, so start from a stored zero.
        let idx = Index { n_no_coor: Some(0), ..Index::default() };
        let mut bytes = Vec::new();
        idx.write_to(&mut bytes).unwrap();
        let back = Index::read_from(&bytes[..]).unwrap();
        assert_eq!(back, idx);
    }

    #[test]
    fn index_rejects_bad_magic_and_truncation() {
        assert!(Index::read_from(&b"NOPE"[..]).is_err());
        assert!(Index::read_from(&b"BAI\x01\x05\x00"[..]).is_err());
    }
}
