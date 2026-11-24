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

/// FNV-1a 64 over the projected records, each followed by a newline, sorted.
fn digest_records(lines: &mut [String]) -> (usize, u64) {
    lines.sort();
    let mut hash: u64 = 0xcbf2_9ce4_8422_2325;
    for line in lines.iter() {
        for &b in line.as_bytes().iter().chain(std::iter::once(&b'\n')) {
            hash ^= u64::from(b);
            hash = hash.wrapping_mul(0x0100_0000_01b3);
        }
    }
    (lines.len(), hash)
}

/// Decodes a CRAM file, with the reference when one is given, and digests its records.
fn digest_cram(name: &str, with_reference: bool) -> (usize, u64) {
    let file = BufReader::new(File::open(fixture(name)).unwrap());
    let mut reader = if with_reference {
        let reference = CramReference::from_path(fixture("ref.fa")).unwrap();
        CramReader::with_reference(file, reference).unwrap()
    } else {
        CramReader::new(file).unwrap()
    };
    let mut sam = SamWriter::new(Vec::new(), reader.header().clone()).unwrap();
    while let Some(record) = reader.next_record().unwrap() {
        sam.write_record(&record).unwrap();
    }
    let text = String::from_utf8(sam.finish().unwrap()).unwrap();
    let mut lines: Vec<String> =
        text.lines().filter(|l| !l.starts_with('@')).map(project).collect();
    digest_records(&mut lines)
}

// The two files below hold the same 3000 alignments. The digests were computed from
// `samtools view -h -T ref.fa FILE`, projected as above, sorted, and hashed with FNV-1a 64.
// samtools 1.22.1 wrote them with `use_tok=0` (names are not tokenized), and
// `level=1` (N32 blocks, order 0 and order 1) and `level=5` (PACK, RLE, N32 with order 1).

#[test]
fn decodes_samtools_cram_31_n32_blocks() {
    assert_eq!(digest_cram("samtools_v31_n32.cram", true), (3000, 0xde3d_334c_10df_af0d));
}

#[test]
fn decodes_samtools_cram_31_pack_n32_order1_blocks() {
    assert_eq!(digest_cram("samtools_v31_pack_n32.cram", true), (3000, 0xde3d_334c_10df_af0d));
}

#[test]
fn decodes_reference_free_samtools_cram_31() {
    // Written by `samtools view -C --output-fmt-option version=3.1 --output-fmt-option
    // no_ref=1` from the 250 alignments, with no reference. Read without one.
    assert_eq!(digest_cram("samtools_v31_noref.cram", false), (250, 0xb713_34f4_f981_0ce9));
}

#[test]
fn decodes_samtools_cram_31_tok3_names_to_samtools_records() {
    // Written by `samtools view -C -T ref.fa --output-fmt-option version=3.1` from the
    // same 250 alignments, with Illumina-style read names. samtools 1.22 stores these
    // names with the TOK3 method (block method 8), which the reader decodes.
    // `expected_tok3_records.txt` is samtools' own decode of that file, projected and sorted.
    let reference = CramReference::from_path(fixture("ref.fa")).unwrap();
    let file = BufReader::new(File::open(fixture("samtools_v31_tok3.cram")).unwrap());
    let mut reader = CramReader::with_reference(file, reference).unwrap();
    let mut sam = SamWriter::new(Vec::new(), reader.header().clone()).unwrap();
    while let Some(record) = reader.next_record().unwrap() {
        sam.write_record(&record).unwrap();
    }
    let text = String::from_utf8(sam.finish().unwrap()).unwrap();
    let mut got: Vec<String> = text.lines().filter(|l| !l.starts_with('@')).map(project).collect();
    got.sort();
    let expected = std::fs::read_to_string(fixture("expected_tok3_records.txt")).unwrap();
    let want: Vec<&str> = expected.lines().collect();
    assert_eq!(got.len(), want.len(), "record count");
    for (g, w) in got.iter().zip(want) {
        assert_eq!(g, w);
    }
}

/// Decodes a TOK3 fixture and checks its records against samtools' own projected decode.
fn assert_tok3_fixture_matches(cram: &str, expected: &str) {
    let reference = CramReference::from_path(fixture("ref.fa")).unwrap();
    let file = BufReader::new(File::open(fixture(cram)).unwrap());
    let mut reader = CramReader::with_reference(file, reference).unwrap();
    let mut sam = SamWriter::new(Vec::new(), reader.header().clone()).unwrap();
    while let Some(record) = reader.next_record().unwrap() {
        sam.write_record(&record).unwrap();
    }
    let text = String::from_utf8(sam.finish().unwrap()).unwrap();
    let mut got: Vec<String> = text.lines().filter(|l| !l.starts_with('@')).map(project).collect();
    got.sort();
    let expected = std::fs::read_to_string(fixture(expected)).unwrap();
    let want: Vec<&str> = expected.lines().collect();
    assert_eq!(got.len(), want.len(), "record count");
    for (g, w) in got.iter().zip(want) {
        assert_eq!(g, w);
    }
}

#[test]
fn decodes_samtools_cram_31_tok3_illumina_old_names() {
    // `HWI-ST123:4:1101:1234:5678`-style names; same 250 alignments as the base TOK3 fixture.
    assert_tok3_fixture_matches(
        "samtools_v31_tok3_illum_old.cram",
        "expected_tok3_illum_old_records.txt",
    );
}

#[test]
fn decodes_samtools_cram_31_tok3_illumina_new_names() {
    // `A00123:45:HXXXXDSXX:1:1101:12345:1000`-style names.
    assert_tok3_fixture_matches(
        "samtools_v31_tok3_illum_new.cram",
        "expected_tok3_illum_new_records.txt",
    );
}

#[test]
fn decodes_samtools_cram_31_tok3_sra_names() {
    // `SRR1234567.1234` names with increasing numbers.
    assert_tok3_fixture_matches("samtools_v31_tok3_sra.cram", "expected_tok3_sra_records.txt");
}

#[test]
fn decodes_samtools_cram_31_tok3_sra_mate_names() {
    // `SRR1234567.1234 1234/1` names: mates share a base number and differ in the suffix.
    assert_tok3_fixture_matches(
        "samtools_v31_tok3_sra_mate.cram",
        "expected_tok3_sra_mate_records.txt",
    );
}

#[test]
fn decodes_samtools_cram_31_tok3_numeric_names() {
    // Purely numeric names `1`, `2`, ...
    assert_tok3_fixture_matches(
        "samtools_v31_tok3_numeric.cram",
        "expected_tok3_numeric_records.txt",
    );
}

#[test]
fn decodes_samtools_cram_31_tok3_duplicate_and_zero_delta_names() {
    // `pair:0001`-style names, each one repeated at distance 2 (`pair:0001` twice, then
    // `pair:0002` twice, ...). The zero-padded counter steps by one per pair. samtools
    // stores the names with TOK3 and they use the DUP (whole-name repeat), DDELTA0
    // (zero-padded numeric delta) and MATCH token types.
    assert_tok3_fixture_matches("samtools_v31_tok3_dup.cram", "expected_tok3_dup_records.txt");
}

#[test]
fn decodes_samtools_cram_31_arithmetic_tok3_names() {
    // Written by `samtools view -C -T ref.fa --output-fmt-option version=3.1 --output-fmt-option
    // use_arith=1 names_illum.sam` (same names as the base TOK3 fixture). The name blocks are
    // TOK3 with the arithmetic flag set, and the data series are arithmetic blocks (method 6).
    assert_tok3_fixture_matches("samtools_v31_tok3_arith.cram", "expected_tok3_arith_records.txt");
}
