//! Pure Rust sequencing I/O for tpt-bio: BGZF, FASTQ, SAM, BAM and BAI.
//!
//! No C code is linked. Compression uses `flate2` with the `miniz_oxide` backend.
//! The only `unsafe` code is the memory-mapping wrapper in [`mmap`], enabled by
//! the default `mmap` feature.

pub mod bai;
pub mod bam;
pub mod bcf;
pub mod bgzf;
pub mod cram;
pub mod error;
pub mod fasta;
pub mod fastq;
pub mod header;
#[cfg(feature = "mmap")]
pub mod mmap;
pub mod ops;
pub mod record;
pub mod sam;
pub mod tabix;
pub mod table;
pub mod vcf;

pub use bai::{Index, build_index, fetch};
pub use bam::{BamReader, BamWriter, RecordView};
pub use bcf::{BcfHeader, BcfReader, BcfWriter};
pub use bgzf::{BgzfReader, BgzfWriter, VirtualOffset};
pub use cram::reader::CramReader;
pub use cram::reference::CramReference;
pub use cram::writer::CramWriter;
pub use error::{Result, SeqError};
pub use fasta::{Fai, FastaRecord};
pub use fastq::{FastqBuf, FastqParser, FastqReader, FastqRecord, FastqWriter};
pub use header::{Header, Reference};
pub use ops::{Flagstat, MergeSorted, sort_coordinate, sort_name};
pub use record::{CigarKind, CigarOp, Flags, Record, Tag, TagValue};
pub use sam::{SamReader, SamWriter};
pub use tabix::{IndexKind, TabixConfig, TabixIndex};
pub use table::{records_to_csv, records_to_jsonl, records_to_table};
pub use vcf::{MetaLine, VcfHeader, VcfReader, VcfRecord, VcfWriter};
