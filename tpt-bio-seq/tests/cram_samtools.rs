//! Checks CRAM against files written by samtools 1.22.
//!
//! The fixtures under `tests/data/samtools/` are CRAM 3.0 output of `samtools view -C`
//! for 250 generated alignments (mismatches, insertions, deletions, soft clips, reverse
//! strands, unmapped reads), with the reference FASTA beside them. `expected_records.txt`
//! is samtools' own decode of that file. The tests need no samtools installation.

use std::fs::File;
use std::io::BufReader;
use std::path::PathBuf;

use tpt_bio_seq::{CramReader, CramReference, CramWriter, Header, SamWriter};

fn fixture(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/data/samtools").join(name)
}

/// Keeps the eleven mandatory columns, with the tags sorted so that order does not matter.
fn project(line: &str) -> String {
    let fields: Vec<&str> = line.split('\t').collect();
    let mut tags: Vec<&str> = fields[11..].to_vec();
    tags.sort_unstable();
    let mut out: Vec<&str> = fields[..11].to_vec();
    out.extend(tags);
    out.join("\t")
}

#[test]
fn decodes_samtools_cram_30_to_samtools_records() {
    let reference = CramReference::from_path(fixture("ref.fa")).unwrap();
    let file = BufReader::new(File::open(fixture("samtools_v30.cram")).unwrap());
    let mut reader = CramReader::with_reference(file, reference).unwrap();
    let mut sam = SamWriter::new(Vec::new(), reader.header().clone()).unwrap();
    while let Some(record) = reader.next_record().unwrap() {
        sam.write_record(&record).unwrap();
    }
    let text = String::from_utf8(sam.finish().unwrap()).unwrap();

    let mut got: Vec<String> = text.lines().filter(|l| !l.starts_with('@')).map(project).collect();
    got.sort();
    let expected = std::fs::read_to_string(fixture("expected_records.txt")).unwrap();
    let want: Vec<&str> = expected.lines().collect();
    assert_eq!(got.len(), want.len(), "record count");
    for (g, w) in got.iter().zip(want) {
        assert_eq!(g, w);
    }
}

#[test]
fn decodes_samtools_cram_31_to_samtools_records() {
    // The same 250 alignments as the 3.0 file, written by samtools with its default
    // codecs: rANS Nx16 with PACK and CAT, and plain order-0 blocks.
    let reference = CramReference::from_path(fixture("ref.fa")).unwrap();
    let file = BufReader::new(File::open(fixture("samtools_v31.cram")).unwrap());
    let mut reader = CramReader::with_reference(file, reference).unwrap();
    let mut sam = SamWriter::new(Vec::new(), reader.header().clone()).unwrap();
    while let Some(record) = reader.next_record().unwrap() {
        sam.write_record(&record).unwrap();
    }
    let text = String::from_utf8(sam.finish().unwrap()).unwrap();
    let mut got: Vec<String> = text.lines().filter(|l| !l.starts_with('@')).map(project).collect();
    got.sort();
    let expected = std::fs::read_to_string(fixture("expected_records.txt")).unwrap();
    let want: Vec<&str> = expected.lines().collect();
    assert_eq!(got.len(), want.len(), "record count");
    for (g, w) in got.iter().zip(want) {
        assert_eq!(g, w);
    }
}

/// The end-of-file container samtools writes for CRAM 3.0, byte for byte.
const SAMTOOLS_EOF_V3: [u8; 38] = [
    0x0f, 0x00, 0x00, 0x00, 0xff, 0xff, 0xff, 0xff, 0x0f, 0xe0, 0x45, 0x4f, 0x46, 0x00, 0x00, 0x00,
    0x00, 0x01, 0x00, 0x05, 0xbd, 0xd9, 0x4f, 0x00, 0x01, 0x00, 0x06, 0x06, 0x01, 0x00, 0x01, 0x00,
    0x01, 0x00, 0xee, 0x63, 0x01, 0x4b,
];

#[test]
fn writer_ends_with_samtools_eof_container() {
    let header = Header::parse("@HD\tVN:1.6\tSO:unsorted\n@SQ\tSN:chr1\tLN:100\n").unwrap();
    let bytes = CramWriter::new(Vec::new(), header, None).unwrap().finish().unwrap();
    assert_eq!(&bytes[bytes.len() - 38..], &SAMTOOLS_EOF_V3[..]);
}
