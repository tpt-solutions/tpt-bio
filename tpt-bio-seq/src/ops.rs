//! Record-level operations: sorting, merging sorted streams, and flagstat counts.

use std::cmp::Ordering;
use std::collections::BinaryHeap;

use crate::error::Result;
use crate::record::{Flags, Record};

/// Sort key for coordinate order. Unplaced reads sort last, as in samtools.
fn coordinate_key(rec: &Record) -> (u32, i64) {
    let ref_key = if rec.ref_id < 0 { u32::MAX } else { rec.ref_id as u32 };
    let pos = if rec.pos < 0 { i64::MAX } else { rec.pos };
    (ref_key, pos)
}

/// Sorts records by reference and position. The sort is stable.
pub fn sort_coordinate(records: &mut [Record]) {
    records.sort_by_key(coordinate_key);
}

/// Sorts records by read name, then by segment (READ1 before READ2).
///
/// Name order is byte-wise. It is not the natural order used by some tools,
/// so fixmate-style pairing should compare names with this same ordering.
pub fn sort_name(records: &mut [Record]) {
    records.sort_by(|a, b| a.name.cmp(&b.name).then_with(|| segment_rank(a).cmp(&segment_rank(b))));
}

fn segment_rank(rec: &Record) -> u8 {
    if rec.flags.contains(Flags::READ1) {
        0
    } else if rec.flags.contains(Flags::READ2) {
        1
    } else {
        2
    }
}

/// A record wrapped so `BinaryHeap` pops the smallest coordinate first.
struct HeapItem {
    key: (u32, i64),
    source: usize,
    record: Record,
}

impl PartialEq for HeapItem {
    fn eq(&self, other: &Self) -> bool {
        self.key == other.key && self.source == other.source
    }
}
impl Eq for HeapItem {}
impl PartialOrd for HeapItem {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}
impl Ord for HeapItem {
    fn cmp(&self, other: &Self) -> Ordering {
        // Reversed, so the max-heap yields the minimum key. Ties go to the lower
        // source index, which keeps the merge stable.
        other.key.cmp(&self.key).then_with(|| other.source.cmp(&self.source))
    }
}

/// Merges coordinate-sorted record streams into one coordinate-sorted stream.
///
/// Each input must already be sorted by [`sort_coordinate`]'s order. Ties keep
/// the input order, so the output is stable across sources.
pub struct MergeSorted<I: Iterator<Item = Result<Record>>> {
    sources: Vec<I>,
    heap: BinaryHeap<HeapItem>,
    primed: bool,
    failed: Option<crate::error::SeqError>,
}

impl<I: Iterator<Item = Result<Record>>> std::fmt::Debug for MergeSorted<I> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("MergeSorted").field("sources", &self.sources.len()).finish_non_exhaustive()
    }
}

impl<I: Iterator<Item = Result<Record>>> MergeSorted<I> {
    /// Creates a merge over the given sources.
    pub fn new(sources: Vec<I>) -> Self {
        Self { sources, heap: BinaryHeap::new(), primed: false, failed: None }
    }

    fn refill(&mut self, source: usize) -> Result<()> {
        if let Some(next) = self.sources[source].next() {
            let record = next?;
            self.heap.push(HeapItem { key: coordinate_key(&record), source, record });
        }
        Ok(())
    }
}

impl<I: Iterator<Item = Result<Record>>> Iterator for MergeSorted<I> {
    type Item = Result<Record>;

    fn next(&mut self) -> Option<Self::Item> {
        if let Some(e) = self.failed.take() {
            return Some(Err(e));
        }
        if !self.primed {
            self.primed = true;
            for s in 0..self.sources.len() {
                if let Err(e) = self.refill(s) {
                    return Some(Err(e));
                }
            }
        }
        let item = self.heap.pop()?;
        if let Err(e) = self.refill(item.source) {
            self.failed = Some(e);
        }
        Some(Ok(item.record))
    }
}

/// Counts in the style of `samtools flagstat`.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct Flagstat {
    /// All records.
    pub total: u64,
    /// Records failing QC (FLAG 0x200).
    pub qc_failed: u64,
    /// Secondary alignments (0x100).
    pub secondary: u64,
    /// Supplementary alignments (0x800).
    pub supplementary: u64,
    /// Duplicates (0x400).
    pub duplicates: u64,
    /// Mapped records (not 0x4).
    pub mapped: u64,
    /// Records flagged as paired (0x1).
    pub paired: u64,
    /// First segments (0x40).
    pub read1: u64,
    /// Last segments (0x80).
    pub read2: u64,
    /// Properly paired (0x2).
    pub properly_paired: u64,
    /// Paired, mapped, and mate mapped.
    pub both_mapped: u64,
    /// Paired, mapped, but mate unmapped.
    pub singletons: u64,
}

impl Flagstat {
    /// Adds one record.
    pub fn add(&mut self, rec: &Record) {
        let f = rec.flags;
        self.total += 1;
        if f.contains(Flags::QC_FAIL) {
            self.qc_failed += 1;
        }
        if f.contains(Flags::SECONDARY) {
            self.secondary += 1;
        }
        if f.contains(Flags::SUPPLEMENTARY) {
            self.supplementary += 1;
        }
        if f.contains(Flags::DUPLICATE) {
            self.duplicates += 1;
        }
        let mapped = !f.contains(Flags::UNMAPPED);
        if mapped {
            self.mapped += 1;
        }
        if f.contains(Flags::PAIRED) {
            self.paired += 1;
            if f.contains(Flags::READ1) {
                self.read1 += 1;
            }
            if f.contains(Flags::READ2) {
                self.read2 += 1;
            }
            if f.contains(Flags::PROPER_PAIR) {
                self.properly_paired += 1;
            }
            if mapped {
                if f.contains(Flags::MATE_UNMAPPED) {
                    self.singletons += 1;
                } else {
                    self.both_mapped += 1;
                }
            }
        }
    }

    /// Counts every record in an iterator.
    pub fn from_records<'a>(records: impl IntoIterator<Item = &'a Record>) -> Self {
        let mut s = Self::default();
        for r in records {
            s.add(r);
        }
        s
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::record::{CigarKind, CigarOp};

    fn rec(ref_id: i32, pos: i64, name: &str, flags: u16) -> Record {
        Record {
            name: name.as_bytes().to_vec(),
            ref_id,
            pos,
            flags: Flags(flags),
            cigar: vec![CigarOp { len: 10, kind: CigarKind::Match }],
            ..Record::default()
        }
    }

    #[test]
    fn coordinate_sort_puts_unplaced_last() {
        let mut v =
            vec![rec(-1, -1, "u", 4), rec(1, 50, "b", 0), rec(0, 900, "c", 0), rec(1, 10, "a", 0)];
        sort_coordinate(&mut v);
        let names: Vec<_> =
            v.iter().map(|r| String::from_utf8_lossy(&r.name).into_owned()).collect();
        assert_eq!(names, ["c", "a", "b", "u"]);
    }

    #[test]
    fn name_sort_orders_mates_together() {
        let mut v = vec![rec(0, 1, "r2", 0x80), rec(0, 2, "r1", 0x80), rec(0, 3, "r1", 0x40)];
        sort_name(&mut v);
        assert_eq!(v[0].name, b"r1");
        assert!(v[0].flags.contains(Flags::READ1));
        assert!(v[1].flags.contains(Flags::READ2));
        assert_eq!(v[2].name, b"r2");
    }

    #[test]
    fn merge_produces_sorted_output() {
        let a = vec![rec(0, 1, "a1", 0), rec(0, 5, "a2", 0), rec(1, 2, "a3", 0)];
        let b = vec![rec(0, 3, "b1", 0), rec(1, 1, "b2", 0)];
        let sources: Vec<_> = vec![a, b].into_iter().map(|v| v.into_iter().map(Ok)).collect();
        let merged: Vec<Record> = MergeSorted::new(sources).map(|r| r.unwrap()).collect();
        let names: Vec<_> =
            merged.iter().map(|r| String::from_utf8_lossy(&r.name).into_owned()).collect();
        assert_eq!(names, ["a1", "b1", "a2", "b2", "a3"]);
    }

    #[test]
    fn merge_propagates_errors() {
        let good: Vec<Result<Record>> = vec![Ok(rec(0, 1, "a", 0))];
        let bad: Vec<Result<Record>> = vec![Err(crate::error::SeqError::malformed("SAM", "boom"))];
        let mut m = MergeSorted::new(vec![good.into_iter(), bad.into_iter()]);
        assert!(m.next().unwrap().is_err());
    }

    #[test]
    fn flagstat_counts() {
        let records = vec![
            rec(0, 1, "p", Flags::PAIRED | Flags::READ1 | Flags::PROPER_PAIR),
            rec(0, 2, "p", Flags::PAIRED | Flags::READ2 | Flags::MATE_UNMAPPED),
            rec(-1, -1, "u", Flags::UNMAPPED | Flags::DUPLICATE),
            rec(0, 3, "s", Flags::SECONDARY | Flags::QC_FAIL),
        ];
        let s = Flagstat::from_records(&records);
        assert_eq!(s.total, 4);
        assert_eq!(s.mapped, 3);
        assert_eq!(s.paired, 2);
        assert_eq!(s.read1, 1);
        assert_eq!(s.read2, 1);
        assert_eq!(s.properly_paired, 1);
        assert_eq!(s.singletons, 1);
        assert_eq!(s.duplicates, 1);
        assert_eq!(s.secondary, 1);
        assert_eq!(s.qc_failed, 1);
    }
}
