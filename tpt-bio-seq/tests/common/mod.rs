//! Shared fixtures and parser entry points for the malformed-input and random-input tests.
//!
//! Every entry point discards its result. A parser passes a test by returning, with `Ok`
//! or `Err`, and fails it by panicking.

#![allow(dead_code)]

use std::io::{Cursor, Write};

use tpt_bio_seq::bai;
use tpt_bio_seq::bcf::{BcfHeader, BcfWriter};
use tpt_bio_seq::bgzf::{BgzfReader, BgzfWriter};
use tpt_bio_seq::fasta::{Fai, FastaReader};
use tpt_bio_seq::sam::SamReader;
use tpt_bio_seq::tabix::{IndexKind, TabixConfig, TabixIndex};
use tpt_bio_seq::{BamReader, BamWriter, BcfReader, FastqParser, VcfReader, VcfRecord};

pub(crate) const VCF_TEXT: &str = "##fileformat=VCFv4.3\n\
##contig=<ID=chr1,length=1000>\n\
##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Depth\">\n\
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n\
##FILTER=<ID=q10,Description=\"Quality below 10\">\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1\ts2\n\
chr1\t100\trs1\tA\tG\t50\tPASS\tDP=10\tGT\t0/1\t1/1\n\
chr1\t200\t.\tAC\tA,T\t.\t.\tDP=7\tGT\t0|1\t.|.\n\
chr1\t300\t.\tG\tC\t20\tq10\tDP=3\tGT\t1/1\t0/0\n";

pub(crate) const SAM_TEXT: &str = "@HD\tVN:1.6\tSO:coordinate\n\
@SQ\tSN:chr1\tLN:1000\n\
@SQ\tSN:chr2\tLN:500\n\
r1\t0\tchr1\t10\t60\t4M\t*\t0\t0\tACGT\tIIII\tNM:i:0\n\
r2\t16\tchr1\t50\t30\t10M\t*\t0\t0\tACGTACGTAC\tIIIIIIIIII\tNM:i:1\n\
r3\t0\tchr2\t7\t60\t3M2I1M\t*\t0\t0\tACGTTA\tIIIIII\tXA:Z:foo\n";

pub(crate) const FASTQ_TEXT: &[u8] = b"@r1 extra\nACGT\n+\nIIII\n@r2\nNNA\n+\n!!#\n";

pub(crate) const FASTA_TEXT: &[u8] = b">chr1 test contig\nACGTACGT\nACGTAC\n>chr2\nGGGG\nCC\n";

/// Compresses `payload` as a BGZF stream (with EOF marker).
///
/// Uses stored (uncompressed) deflate blocks. The harness wraps thousands of mutated
/// payloads, and the debug-build compressor would dominate the run time.
pub(crate) fn bgzf_wrap(payload: &[u8]) -> Vec<u8> {
    let mut w = BgzfWriter::with_level(Vec::new(), flate2::Compression::none());
    w.write_all(payload).expect("in-memory write");
    w.finish().expect("in-memory finish")
}

/// Decompresses a BGZF stream. Returns `None` if the bytes are not valid BGZF.
pub(crate) fn bgzf_unwrap(data: &[u8]) -> Option<Vec<u8>> {
    BgzfReader::new(data).read_to_vec().ok()
}

pub(crate) fn bam_bytes() -> Vec<u8> {
    let mut reader = SamReader::new(SAM_TEXT.as_bytes()).expect("valid SAM");
    let header = reader.header().clone();
    let mut out = Vec::new();
    let mut w = BamWriter::new(&mut out, header).expect("BAM header");
    while let Some(rec) = reader.next_record().expect("valid SAM record") {
        w.write_record(&rec).expect("BAM record");
    }
    w.finish().expect("BAM finish");
    out
}

pub(crate) fn bai_bytes() -> Vec<u8> {
    let mut reader = BamReader::new(Cursor::new(bam_bytes())).expect("valid BAM");
    let index = bai::build_index(&mut reader).expect("BAI build");
    let mut out = Vec::new();
    index.write_to(&mut out).expect("BAI write");
    out
}

pub(crate) fn vcf_plain() -> Vec<u8> {
    VCF_TEXT.as_bytes().to_vec()
}

pub(crate) fn bcf_header() -> BcfHeader {
    let mut header_only = String::new();
    for line in VCF_TEXT.lines().filter(|l| l.starts_with('#')) {
        header_only.push_str(line);
        header_only.push('\n');
    }
    BcfHeader::parse(&header_only).expect("valid BCF header")
}

fn vcf_record_lines() -> Vec<&'static str> {
    VCF_TEXT.lines().filter(|l| !l.starts_with('#')).collect()
}

/// Uncompressed BCF.
pub(crate) fn bcf_plain() -> Vec<u8> {
    let header = bcf_header();
    let mut w = BcfWriter::new(Vec::new(), &header).expect("BCF header");
    for (i, line) in vcf_record_lines().iter().enumerate() {
        let rec = VcfRecord::parse(line, header.n_samples(), i as u64 + 1).expect("VCF record");
        w.write_record(&rec).expect("BCF record");
    }
    w.into_inner().expect("BCF finish")
}

/// BGZF-compressed BCF (the usual form).
pub(crate) fn bcf_bgzf() -> Vec<u8> {
    let header = bcf_header();
    let mut w = BcfWriter::new_bgzf(Vec::new(), &header).expect("BCF header");
    for (i, line) in vcf_record_lines().iter().enumerate() {
        let rec = VcfRecord::parse(line, header.n_samples(), i as u64 + 1).expect("VCF record");
        w.write_record(&rec).expect("BCF record");
    }
    w.finish().expect("BCF finish")
}

/// Tabix (TBI) index over the BGZF-compressed VCF.
pub(crate) fn tbi_bytes() -> Vec<u8> {
    let idx = TabixIndex::build_text(
        &bgzf_wrap(VCF_TEXT.as_bytes())[..],
        TabixConfig::vcf(),
        IndexKind::Tbi,
    )
    .expect("TBI build");
    let mut out = Vec::new();
    idx.write_to(&mut out).expect("TBI write");
    out
}

/// CSI index over the BGZF-compressed VCF.
pub(crate) fn csi_bytes() -> Vec<u8> {
    let kind = IndexKind::Csi { min_shift: 14, depth: 5 };
    let idx = TabixIndex::build_text(&bgzf_wrap(VCF_TEXT.as_bytes())[..], TabixConfig::vcf(), kind)
        .expect("CSI build");
    let mut out = Vec::new();
    idx.write_to(&mut out).expect("CSI write");
    out
}

pub(crate) fn fai_text() -> Vec<u8> {
    let fai = Fai::build(FASTA_TEXT).expect("FAI build");
    let mut out = Vec::new();
    fai.write_to(&mut out).expect("FAI write");
    out
}

// ---------------------------------------------------------------------------
// Parser entry points. Each one runs a parser to completion and discards results.
// ---------------------------------------------------------------------------

pub(crate) type ParseFn = fn(&[u8]);

pub(crate) fn parse_bgzf_parallel(d: &[u8]) {
    let _ = tpt_bio_seq::bgzf::decompress_parallel(d);
}

pub(crate) fn parse_bgzf_stream(d: &[u8]) {
    let _ = BgzfReader::new(d).read_to_vec();
}

pub(crate) fn parse_bam(d: &[u8]) {
    if let Ok(mut r) = BamReader::new(d) {
        while let Ok(Some(_)) = r.next_record() {}
    }
}

pub(crate) fn parse_sam(d: &[u8]) {
    if let Ok(mut r) = SamReader::new(d) {
        while let Ok(Some(_)) = r.next_record() {}
    }
}

pub(crate) fn parse_fastq(d: &[u8]) {
    for rec in FastqParser::new(d) {
        if rec.is_err() {
            break;
        }
    }
}

pub(crate) fn parse_vcf(d: &[u8]) {
    if let Ok(mut r) = VcfReader::new(d) {
        while let Ok(Some(rec)) = r.next_record() {
            let _ = rec.samples().count();
        }
    }
}

pub(crate) fn parse_vcf_bgzf(d: &[u8]) {
    if let Ok(mut r) = VcfReader::new_bgzf(d) {
        while let Ok(Some(_)) = r.next_record() {}
    }
}

pub(crate) fn parse_bcf(d: &[u8]) {
    if let Ok(mut r) = BcfReader::new(d) {
        while let Ok(Some(_)) = r.next_record() {}
    }
}

pub(crate) fn parse_bcf_bgzf(d: &[u8]) {
    if let Ok(mut r) = BcfReader::new_bgzf(d) {
        while let Ok(Some(_)) = r.next_record() {}
    }
}

pub(crate) fn parse_tabix(d: &[u8]) {
    let _ = TabixIndex::read_from(d);
}

pub(crate) fn parse_bai(d: &[u8]) {
    let _ = bai::Index::read_from(d);
}

pub(crate) fn parse_fasta(d: &[u8]) {
    let mut r = FastaReader::new(d);
    while let Ok(Some(_)) = r.next_record() {}
}

pub(crate) fn parse_fai(d: &[u8]) {
    if let Ok(text) = std::str::from_utf8(d) {
        let _ = Fai::read_from(text);
    }
}

/// Every parser, with a stable name for failure reports.
pub(crate) const PARSERS: &[(&str, ParseFn)] = &[
    ("bgzf::decompress_parallel", parse_bgzf_parallel),
    ("bgzf::BgzfReader", parse_bgzf_stream),
    ("bam::BamReader", parse_bam),
    ("sam::SamReader", parse_sam),
    ("fastq::FastqParser", parse_fastq),
    ("vcf::VcfReader", parse_vcf),
    ("vcf::VcfReader::new_bgzf", parse_vcf_bgzf),
    ("bcf::BcfReader", parse_bcf),
    ("bcf::BcfReader::new_bgzf", parse_bcf_bgzf),
    ("tabix::TabixIndex::read_from", parse_tabix),
    ("bai::Index::read_from", parse_bai),
    ("fasta::FastaReader", parse_fasta),
    ("fasta::Fai::read_from", parse_fai),
];

/// Looks up a parser by its name in [`PARSERS`].
pub(crate) fn parser(name: &str) -> ParseFn {
    PARSERS.iter().find(|(n, _)| *n == name).map(|(_, f)| *f).expect("known parser name")
}
