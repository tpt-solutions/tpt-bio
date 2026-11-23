//! Alignment record model shared by SAM and BAM.
//!
//! Coordinates are stored as in BAM: 0-based positions, reference IDs as indexes
//! into the header's reference list, and -1 for "unset". SAM text is converted
//! to and from this form by [`crate::sam`], and BAM binary by [`crate::bam`].

use std::fmt;

/// SAM FLAG bits.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Default)]
pub struct Flags(pub u16);

impl Flags {
    /// Template has multiple segments.
    pub const PAIRED: u16 = 0x1;
    /// Each segment is properly aligned.
    pub const PROPER_PAIR: u16 = 0x2;
    /// Segment is unmapped.
    pub const UNMAPPED: u16 = 0x4;
    /// Next segment is unmapped.
    pub const MATE_UNMAPPED: u16 = 0x8;
    /// Sequence is reverse complemented.
    pub const REVERSE: u16 = 0x10;
    /// Next segment is reverse complemented.
    pub const MATE_REVERSE: u16 = 0x20;
    /// First segment in the template.
    pub const READ1: u16 = 0x40;
    /// Last segment in the template.
    pub const READ2: u16 = 0x80;
    /// Secondary alignment.
    pub const SECONDARY: u16 = 0x100;
    /// Failed quality checks.
    pub const QC_FAIL: u16 = 0x200;
    /// PCR or optical duplicate.
    pub const DUPLICATE: u16 = 0x400;
    /// Supplementary alignment.
    pub const SUPPLEMENTARY: u16 = 0x800;

    /// Whether a bit is set.
    #[must_use]
    pub const fn contains(self, bit: u16) -> bool {
        self.0 & bit != 0
    }
}

/// CIGAR operation kind.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum CigarKind {
    /// Alignment match (M).
    Match,
    /// Insertion to the reference (I).
    Insertion,
    /// Deletion from the reference (D).
    Deletion,
    /// Skipped region from the reference (N).
    Skip,
    /// Soft clip (S).
    SoftClip,
    /// Hard clip (H).
    HardClip,
    /// Padding (P).
    Pad,
    /// Sequence match (=).
    SeqMatch,
    /// Sequence mismatch (X).
    SeqMismatch,
}

impl CigarKind {
    /// The 4-bit BAM code (0..=8).
    #[must_use]
    pub const fn to_code(self) -> u8 {
        match self {
            Self::Match => 0,
            Self::Insertion => 1,
            Self::Deletion => 2,
            Self::Skip => 3,
            Self::SoftClip => 4,
            Self::HardClip => 5,
            Self::Pad => 6,
            Self::SeqMatch => 7,
            Self::SeqMismatch => 8,
        }
    }

    /// Decodes a 4-bit BAM code.
    #[must_use]
    pub const fn from_code(c: u8) -> Option<Self> {
        Some(match c {
            0 => Self::Match,
            1 => Self::Insertion,
            2 => Self::Deletion,
            3 => Self::Skip,
            4 => Self::SoftClip,
            5 => Self::HardClip,
            6 => Self::Pad,
            7 => Self::SeqMatch,
            8 => Self::SeqMismatch,
            _ => return None,
        })
    }

    /// The SAM character for this operation.
    #[must_use]
    pub const fn to_char(self) -> char {
        b"MIDNSHP=X"[self.to_code() as usize] as char
    }

    /// Parses a SAM CIGAR operation character.
    #[must_use]
    pub const fn from_char(c: u8) -> Option<Self> {
        Some(match c {
            b'M' => Self::Match,
            b'I' => Self::Insertion,
            b'D' => Self::Deletion,
            b'N' => Self::Skip,
            b'S' => Self::SoftClip,
            b'H' => Self::HardClip,
            b'P' => Self::Pad,
            b'=' => Self::SeqMatch,
            b'X' => Self::SeqMismatch,
            _ => return None,
        })
    }

    /// Whether the operation advances along the reference.
    #[must_use]
    pub const fn consumes_reference(self) -> bool {
        matches!(
            self,
            Self::Match | Self::Deletion | Self::Skip | Self::SeqMatch | Self::SeqMismatch
        )
    }

    /// Whether the operation advances along the query sequence.
    #[must_use]
    pub const fn consumes_query(self) -> bool {
        matches!(
            self,
            Self::Match | Self::Insertion | Self::SoftClip | Self::SeqMatch | Self::SeqMismatch
        )
    }
}

/// One CIGAR operation.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct CigarOp {
    /// Operation length.
    pub len: u32,
    /// Operation kind.
    pub kind: CigarKind,
}

/// Reference bases spanned by a CIGAR.
#[must_use]
pub fn reference_length(cigar: &[CigarOp]) -> u64 {
    cigar.iter().filter(|op| op.kind.consumes_reference()).map(|op| u64::from(op.len)).sum()
}

/// Query bases described by a CIGAR (including soft clips).
#[must_use]
pub fn query_length(cigar: &[CigarOp]) -> u64 {
    cigar.iter().filter(|op| op.kind.consumes_query()).map(|op| u64::from(op.len)).sum()
}

/// Formats a CIGAR as text; an empty CIGAR is `*`.
#[must_use]
pub fn format_cigar(cigar: &[CigarOp]) -> String {
    if cigar.is_empty() {
        return "*".to_string();
    }
    let mut s = String::new();
    for op in cigar {
        s.push_str(&op.len.to_string());
        s.push(op.kind.to_char());
    }
    s
}

/// A typed auxiliary tag value.
#[derive(Debug, Clone, PartialEq)]
pub enum TagValue {
    /// Single printable character (`A`).
    Char(u8),
    /// Integer (`i`). BAM stores it in the smallest type that fits.
    Int(i64),
    /// Single-precision float (`f`).
    Float(f32),
    /// Text (`Z`).
    String(Vec<u8>),
    /// Hex-encoded bytes (`H`), stored as uppercase hex text.
    Hex(Vec<u8>),
    /// Integer array (`B` with an integer subtype).
    IntArray(Vec<i64>),
    /// Float array (`B:f`).
    FloatArray(Vec<f32>),
}

/// An auxiliary tag: a two-character key and a typed value.
#[derive(Debug, Clone, PartialEq)]
pub struct Tag {
    /// Two-byte tag name, such as `NM` or `RG`.
    pub key: [u8; 2],
    /// Tag value.
    pub value: TagValue,
}

/// An alignment record in BAM-native form.
#[derive(Debug, Clone, PartialEq)]
pub struct Record {
    /// Read name. Empty means `*`.
    pub name: Vec<u8>,
    /// SAM FLAG.
    pub flags: Flags,
    /// Reference index, or -1 for unmapped.
    pub ref_id: i32,
    /// 0-based leftmost position, or -1 when unset.
    pub pos: i64,
    /// Mapping quality (255 means unavailable).
    pub mapq: u8,
    /// CIGAR operations. Empty means `*`.
    pub cigar: Vec<CigarOp>,
    /// Mate reference index, or -1.
    pub next_ref_id: i32,
    /// Mate 0-based position, or -1.
    pub next_pos: i64,
    /// Observed template length.
    pub tlen: i64,
    /// Bases as ASCII. Empty means `*`.
    pub seq: Vec<u8>,
    /// Raw Phred scores (not +33). Empty means `*`.
    pub qual: Vec<u8>,
    /// Auxiliary tags, in order.
    pub tags: Vec<Tag>,
}

impl Default for Record {
    fn default() -> Self {
        Self {
            name: Vec::new(),
            flags: Flags::default(),
            ref_id: -1,
            pos: -1,
            mapq: 255,
            cigar: Vec::new(),
            next_ref_id: -1,
            next_pos: -1,
            tlen: 0,
            seq: Vec::new(),
            qual: Vec::new(),
            tags: Vec::new(),
        }
    }
}

impl Record {
    /// Whether the read is mapped (FLAG 0x4 not set).
    #[must_use]
    pub fn is_mapped(&self) -> bool {
        !self.flags.contains(Flags::UNMAPPED)
    }

    /// The 0-based exclusive end of the alignment on the reference, or `None`
    /// when unmapped or without a CIGAR. Unmapped reads with a position use the position.
    #[must_use]
    pub fn alignment_end(&self) -> Option<i64> {
        if self.pos < 0 {
            return None;
        }
        if self.cigar.is_empty() {
            return Some(self.pos + 1);
        }
        Some(self.pos + reference_length(&self.cigar) as i64)
    }

    /// Looks up a tag by key.
    #[must_use]
    pub fn tag(&self, key: &[u8; 2]) -> Option<&TagValue> {
        self.tags.iter().find(|t| &t.key == key).map(|t| &t.value)
    }
}

impl fmt::Display for TagValue {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Char(c) => write!(f, "A:{}", *c as char),
            Self::Int(v) => write!(f, "i:{v}"),
            Self::Float(v) => write!(f, "f:{v}"),
            Self::String(s) => write!(f, "Z:{}", String::from_utf8_lossy(s)),
            Self::Hex(s) => write!(f, "H:{}", String::from_utf8_lossy(s)),
            Self::IntArray(values) => {
                write!(f, "B:i")?;
                for v in values {
                    write!(f, ",{v}")?;
                }
                Ok(())
            }
            Self::FloatArray(values) => {
                write!(f, "B:f")?;
                for v in values {
                    write!(f, ",{v}")?;
                }
                Ok(())
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn cigar_char_roundtrip_and_lengths() {
        let text = b"MIDNSHP=X";
        let ops: Vec<CigarKind> = text.iter().map(|&c| CigarKind::from_char(c).unwrap()).collect();
        for (kind, &c) in ops.iter().zip(text) {
            assert_eq!(kind.to_char() as u8, c);
            assert_eq!(CigarKind::from_code(kind.to_code()), Some(*kind));
        }
        let cigar = [
            CigarOp { len: 5, kind: CigarKind::SoftClip },
            CigarOp { len: 10, kind: CigarKind::Match },
            CigarOp { len: 2, kind: CigarKind::Deletion },
            CigarOp { len: 3, kind: CigarKind::Insertion },
        ];
        assert_eq!(reference_length(&cigar), 12);
        assert_eq!(query_length(&cigar), 18);
        assert_eq!(format_cigar(&cigar), "5S10M2D3I");
        assert_eq!(format_cigar(&[]), "*");
    }

    #[test]
    fn alignment_end_uses_reference_span() {
        let rec = Record {
            pos: 100,
            cigar: vec![
                CigarOp { len: 4, kind: CigarKind::Match },
                CigarOp { len: 1, kind: CigarKind::Insertion },
            ],
            ..Record::default()
        };
        assert_eq!(rec.alignment_end(), Some(104));
        assert_eq!(Record::default().alignment_end(), None);
    }
}
