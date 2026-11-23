//! Property-style random testing: every parser must return `Ok` or `Err` on arbitrary bytes.
//!
//! Each parser gets 10,000 cases with inputs of up to 4 KiB, in two strategies:
//! uniformly random bytes, and random bytes after a format's magic number (so the parser
//! gets past its first check more often).

mod common;

use std::panic::{AssertUnwindSafe, catch_unwind};

use common::*;
use proptest::prelude::*;
use proptest::test_runner::{Config, TestCaseError, TestRunner};

const CASES: u32 = 10_000;
const MAX_LEN: usize = 4096;

/// Leading bytes of each format, to get past the first check.
const MAGICS: &[&[u8]] = &[
    b"\x1f\x8b\x08\x04",
    b"BAM\x01",
    b"BCF\x02\x02",
    b"TBI\x01",
    b"CSI\x01",
    b"BAI\x01",
    b"@r1\nACGT\n+\nIIII\n",
    b">chr1\nACGT\n",
    b"##fileformat=VCFv4.3\n#CHROM\tPOS\n",
];

fn check(name: &str, data: &[u8]) -> Result<(), TestCaseError> {
    let parse = parser(name);
    catch_unwind(AssertUnwindSafe(|| parse(data))).map_err(|payload| {
        let message = if let Some(s) = payload.downcast_ref::<&str>() {
            (*s).to_string()
        } else if let Some(s) = payload.downcast_ref::<String>() {
            s.clone()
        } else {
            "non-string panic payload".to_string()
        };
        TestCaseError::fail(format!("{name} panicked on {} bytes: {message}", data.len()))
    })
}

fn run_cases<S>(name: &'static str, strategy: S)
where
    S: Strategy<Value = Vec<u8>>,
{
    let config = Config { cases: CASES, failure_persistence: None, ..Config::default() };
    let mut runner = TestRunner::new(config);
    if let Err(e) = runner.run(&strategy, |data| check(name, &data)) {
        panic!("{name}: {e}");
    }
}

fn uniform_bytes() -> impl Strategy<Value = Vec<u8>> {
    proptest::collection::vec(any::<u8>(), 0..=MAX_LEN)
}

fn after_magic() -> impl Strategy<Value = Vec<u8>> {
    (0..MAGICS.len(), proptest::collection::vec(any::<u8>(), 0..=MAX_LEN)).prop_map(|(i, tail)| {
        let mut out = MAGICS[i].to_vec();
        out.extend_from_slice(&tail);
        out.truncate(MAX_LEN);
        out
    })
}

/// Runs both strategies against one parser.
fn fuzz_parser(name: &'static str) {
    run_cases(name, uniform_bytes());
    run_cases(name, after_magic());
}

#[test]
fn random_bytes_bgzf_parallel() {
    fuzz_parser("bgzf::decompress_parallel");
}

#[test]
fn random_bytes_bgzf_stream() {
    fuzz_parser("bgzf::BgzfReader");
}

#[test]
fn random_bytes_bam() {
    fuzz_parser("bam::BamReader");
}

#[test]
fn random_bytes_sam() {
    fuzz_parser("sam::SamReader");
}

#[test]
fn random_bytes_fastq() {
    fuzz_parser("fastq::FastqParser");
}

#[test]
fn random_bytes_vcf() {
    fuzz_parser("vcf::VcfReader");
}

#[test]
fn random_bytes_vcf_bgzf() {
    fuzz_parser("vcf::VcfReader::new_bgzf");
}

#[test]
fn random_bytes_bcf() {
    fuzz_parser("bcf::BcfReader");
}

#[test]
fn random_bytes_bcf_bgzf() {
    fuzz_parser("bcf::BcfReader::new_bgzf");
}

#[test]
fn random_bytes_tabix() {
    fuzz_parser("tabix::TabixIndex::read_from");
}

#[test]
fn random_bytes_bai() {
    fuzz_parser("bai::Index::read_from");
}

#[test]
fn random_bytes_fasta() {
    fuzz_parser("fasta::FastaReader");
}

#[test]
fn random_bytes_fai() {
    fuzz_parser("fasta::Fai::read_from");
}
