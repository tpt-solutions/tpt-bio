//! Malformed-input corpus: every parser must return `Ok` or `Err`, never panic, and must
//! reject oversized length fields instead of allocating for them.
//!
//! Inputs are valid files built in memory. Each one is mutated deterministically:
//! truncation at every prefix length (up to a cap), bit flips in the first bytes, and
//! `u32` length fields set to huge values. For BGZF-wrapped formats the decompressed
//! payload is also mutated and re-wrapped, so the CRC check does not hide parser bugs.

mod common;

use std::collections::BTreeMap;
use std::fmt::Write as _;
use std::panic::{AssertUnwindSafe, catch_unwind};

use common::*;
use tpt_bio_seq::bgzf::BgzfReader;
use tpt_bio_seq::fasta::Fai;
use tpt_bio_seq::{BamReader, BcfReader, FastqParser, SamReader, TabixIndex, VcfReader};

/// Longest prefix tried for each input (truncation cases).
const PREFIX_CAP: usize = 4096;
/// Number of leading bytes that are each XOR-flipped.
const FLIP_BYTES: usize = 256;
/// Leading offsets at which a `u32` huge length field is written.
const FIELD_OFFSETS: usize = 64;

/// Deterministic mutations of `valid`, as (description, bytes).
fn mutations(valid: &[u8]) -> Vec<(String, Vec<u8>)> {
    let n = valid.len();
    let mut out = Vec::new();
    for len in 0..n.min(PREFIX_CAP) {
        out.push((format!("truncated to {len} of {n} bytes"), valid[..len].to_vec()));
    }
    for i in 0..n.min(FLIP_BYTES) {
        for mask in [0xFFu8, 0x80] {
            let mut v = valid.to_vec();
            v[i] ^= mask;
            out.push((format!("byte {i} xor {mask:#04x}"), v));
        }
    }
    for off in 0..n.saturating_sub(3).min(FIELD_OFFSETS) {
        for value in [u32::MAX, 0x7FFF_FFFF] {
            let mut v = valid.to_vec();
            v[off..off + 4].copy_from_slice(&value.to_le_bytes());
            out.push((format!("u32 at offset {off} set to {value:#010x}"), v));
        }
    }
    out
}

/// Panics seen so far, grouped by parser name.
#[derive(Default)]
struct Failures {
    by_parser: BTreeMap<&'static str, (usize, String)>,
}

impl Failures {
    fn record(&mut self, parser: &'static str, case: &str, message: &str) {
        let entry = self.by_parser.entry(parser).or_insert((0, String::new()));
        entry.0 += 1;
        if entry.1.is_empty() {
            entry.1 = format!("{case}: {message}");
        }
    }

    fn assert_clean(&self) {
        if self.by_parser.is_empty() {
            return;
        }
        let mut msg = String::from("parsers panicked on malformed input:\n");
        for (parser, (count, first)) in &self.by_parser {
            let _ = writeln!(msg, "  {parser}: {count} panicking case(s); first: {first}");
        }
        panic!("{msg}");
    }
}

fn panic_message(payload: &(dyn std::any::Any + Send)) -> String {
    if let Some(s) = payload.downcast_ref::<&str>() {
        (*s).to_string()
    } else if let Some(s) = payload.downcast_ref::<String>() {
        s.clone()
    } else {
        "non-string panic payload".to_string()
    }
}

/// Runs each named parser on `bytes`, recording any panic against `label` and `case`.
fn run_parsers(label: &str, case: &str, bytes: &[u8], parsers: &[&'static str], f: &mut Failures) {
    for &name in parsers {
        let parse = parser(name);
        if let Err(payload) = catch_unwind(AssertUnwindSafe(|| parse(bytes))) {
            f.record(name, &format!("{label}, {case}"), &panic_message(payload.as_ref()));
        }
    }
}

/// Checks raw mutations of `valid`, and payload mutations when it is BGZF.
fn check_input(
    label: &str,
    valid: &[u8],
    raw_parsers: &[&'static str],
    payload_parsers: &[&'static str],
    f: &mut Failures,
) {
    for (case, bytes) in mutations(valid) {
        run_parsers(label, &case, &bytes, raw_parsers, f);
    }
    if let Some(payload) = bgzf_unwrap(valid) {
        for (case, bytes) in mutations(&payload) {
            let wrapped = bgzf_wrap(&bytes);
            let case = format!("decompressed payload, {case}");
            run_parsers(label, &case, &wrapped, payload_parsers, f);
        }
    }
}

#[test]
fn valid_inputs_parse_cleanly() {
    // Guards the fixtures: if a valid input stops parsing, the mutation results mean nothing.
    assert!(bgzf_unwrap(&bam_bytes()).is_some());
    let bam_data = bam_bytes();
    let mut bam = BamReader::new(&bam_data[..]).expect("BAM header");
    let mut n = 0;
    while bam.next_record().expect("BAM record").is_some() {
        n += 1;
    }
    assert_eq!(n, 3);
    let mut sam = SamReader::new(SAM_TEXT.as_bytes()).expect("SAM header");
    assert!(sam.next_record().expect("SAM record").is_some());
    let mut vcf = VcfReader::new(VCF_TEXT.as_bytes()).expect("VCF header");
    assert!(vcf.next_record().expect("VCF record").is_some());
    let bcf_data = bcf_bgzf();
    let mut bcf = BcfReader::new_bgzf(&bcf_data[..]).expect("BCF header");
    assert!(bcf.next_record().expect("BCF record").is_some());
    assert!(TabixIndex::read_from(&tbi_bytes()[..]).is_ok());
    assert!(TabixIndex::read_from(&csi_bytes()[..]).is_ok());
    assert!(tpt_bio_seq::bai::Index::read_from(&bai_bytes()[..]).is_ok());
    assert!(Fai::read_from(std::str::from_utf8(&fai_text()).unwrap()).is_ok());
    assert_eq!(FastqParser::new(FASTQ_TEXT).count(), 2);
}

/// Runs the full mutation set against every format.
///
/// Runtime note: the mutated BAM record size fields (`block_size` in `BamReader`) can claim
/// up to 4 GiB. The reader zero-fills that much before reading, so a few cases take about
/// 10 to 20 seconds each in a debug build and commit that much RAM. This is an
/// over-allocation finding; see the report, not a test bug.
#[test]
fn malformed_inputs_never_panic() {
    let mut f = Failures::default();

    check_input(
        "BAM",
        &bam_bytes(),
        &["bgzf::decompress_parallel", "bgzf::BgzfReader", "bam::BamReader"],
        &["bam::BamReader"],
        &mut f,
    );
    check_input("SAM", SAM_TEXT.as_bytes(), &["sam::SamReader"], &[], &mut f);
    check_input("FASTQ", FASTQ_TEXT, &["fastq::FastqParser"], &[], &mut f);
    check_input(
        "BGZF(FASTQ)",
        &bgzf_wrap(FASTQ_TEXT),
        &["bgzf::decompress_parallel", "bgzf::BgzfReader", "fastq::FastqParser"],
        &["fastq::FastqParser"],
        &mut f,
    );
    check_input("VCF", &vcf_plain(), &["vcf::VcfReader"], &[], &mut f);
    check_input(
        "VCF(BGZF)",
        &bgzf_wrap(VCF_TEXT.as_bytes()),
        &["vcf::VcfReader::new_bgzf", "bgzf::BgzfReader"],
        &["vcf::VcfReader"],
        &mut f,
    );
    check_input("BCF", &bcf_plain(), &["bcf::BcfReader"], &[], &mut f);
    check_input(
        "BCF(BGZF)",
        &bcf_bgzf(),
        &["bcf::BcfReader::new_bgzf", "bgzf::BgzfReader"],
        &["bcf::BcfReader"],
        &mut f,
    );
    check_input(
        "TBI",
        &tbi_bytes(),
        &["tabix::TabixIndex::read_from", "bgzf::BgzfReader"],
        &["tabix::TabixIndex::read_from"],
        &mut f,
    );
    check_input(
        "CSI",
        &csi_bytes(),
        &["tabix::TabixIndex::read_from"],
        &["tabix::TabixIndex::read_from"],
        &mut f,
    );
    check_input("BAI", &bai_bytes(), &["bai::Index::read_from"], &[], &mut f);
    check_input("FASTA", FASTA_TEXT, &["fasta::FastaReader"], &[], &mut f);
    check_input("FAI", &fai_text(), &["fasta::Fai::read_from"], &[], &mut f);

    f.assert_clean();
}

// ---------------------------------------------------------------------------
// Oversized length fields: each case must be rejected with an error.
//
// These check the outcome only. They do not measure allocation, because the crate forbids
// a custom global allocator in tests. The allocation behavior is recorded in the report.
// ---------------------------------------------------------------------------

fn le32(v: u32) -> [u8; 4] {
    v.to_le_bytes()
}

#[test]
fn bam_header_text_length_near_4gib_is_rejected() {
    let mut payload = b"BAM\x01".to_vec();
    payload.extend_from_slice(&le32(0xF000_0000));
    let wrapped = bgzf_wrap(&payload);
    assert!(BamReader::new(&wrapped[..]).is_err());
}

#[test]
fn bam_reference_name_length_near_4gib_is_rejected() {
    let mut payload = b"BAM\x01".to_vec();
    payload.extend_from_slice(&le32(0)); // l_text
    payload.extend_from_slice(&le32(1)); // n_ref
    payload.extend_from_slice(&le32(0xFFFF_FFF0)); // l_name
    let wrapped = bgzf_wrap(&payload);
    assert!(BamReader::new(&wrapped[..]).is_err());
}

#[test]
fn bcf_header_text_length_at_limit_is_rejected() {
    // The declared text is 256 MiB but the stream ends immediately.
    let mut data = b"BCF\x02\x02".to_vec();
    data.extend_from_slice(&le32(0x0FFF_FFFF));
    assert!(BcfReader::new(&data[..]).is_err());
}

// Ignored by default: the reader zero-fills a buffer of the claimed size (up to 2 GiB)
// before it reads, so this run commits about 2 GiB of RAM. Run with `--ignored`.
#[test]
#[ignore = "allocates up to 2 GiB before reading; see comment"]
fn bcf_record_length_near_2gib_is_rejected() {
    let plain = bcf_plain();
    let l_text = u32::from_le_bytes([plain[5], plain[6], plain[7], plain[8]]) as usize;
    let mut data = plain[..9 + l_text].to_vec();
    data.extend_from_slice(&le32(0x7FFF_FFFF)); // l_shared
    data.extend_from_slice(&le32(0)); // l_indiv
    let mut reader = BcfReader::new(&data[..]).expect("valid BCF header");
    assert!(reader.next_record().is_err());
}

#[test]
fn tabix_csi_aux_length_near_2gib_is_rejected() {
    let mut data = b"CSI\x01".to_vec();
    data.extend_from_slice(&le32(14)); // min_shift
    data.extend_from_slice(&le32(5)); // depth
    data.extend_from_slice(&le32(0x7FFF_FFFF)); // l_aux
    assert!(TabixIndex::read_from(&data[..]).is_err());
}

#[test]
fn bai_reference_count_near_2_billion_is_rejected() {
    let mut data = b"BAI\x01".to_vec();
    data.extend_from_slice(&le32(0x7FFF_FFFF));
    assert!(tpt_bio_seq::bai::Index::read_from(&data[..]).is_err());
}

#[test]
fn bgzf_block_claiming_maximum_size_is_rejected() {
    // BC subfield with BSIZE = 0xFFFF, above the 64 KiB block limit.
    let mut block = vec![0x1f, 0x8b, 8, 4, 0, 0, 0, 0, 0, 0xff, 6, 0, b'B', b'C', 2, 0];
    block.extend_from_slice(&[0xFF, 0xFF]);
    block.extend_from_slice(&[0u8; 64]);
    assert!(tpt_bio_seq::bgzf::decompress_parallel(&block).is_err());
    assert!(BgzfReader::new(&block[..]).read_to_vec().is_err());
}
