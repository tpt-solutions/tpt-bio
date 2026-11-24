//! Checks BCF, TBI and CSI against files written by bcftools 1.22 and tabix 1.22.
//!
//! The fixtures under `tests/data/bcftools/` come from `in.vcf`, a generated 240-record
//! VCF with three samples: `bcftools view -O b` made `oracle.bcf`, `bcftools index`
//! made its CSI, and `bgzip` and `tabix -p vcf` made `oracle.vcf.gz` and its TBI. The
//! `expected_*` files are the outputs of `bcftools view -H` and `tabix`. No tool needs
//! to be installed to run these tests.

use std::fs::File;
use std::io::{BufReader, Cursor, Read};
use std::path::PathBuf;

use tpt_bio_seq::bcf::{build_bcf_index, fetch_lines};
use tpt_bio_seq::tabix::fetch_text;
use tpt_bio_seq::{
    BcfHeader, BcfReader, BcfWriter, BgzfReader, IndexKind, TabixConfig, TabixIndex, VcfReader,
};

fn fixture(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/data/bcftools").join(name)
}

fn expected(name: &str) -> Vec<String> {
    std::fs::read_to_string(fixture(name)).unwrap().lines().map(str::to_owned).collect()
}

/// Decompresses a BGZF file completely.
fn gunzip(bytes: &[u8]) -> Vec<u8> {
    let mut out = Vec::new();
    BgzfReader::new(bytes).read_to_end(&mut out).unwrap();
    out
}

fn bcf_lines(bytes: Vec<u8>) -> Vec<String> {
    let mut reader = BcfReader::new_bgzf(Cursor::new(bytes)).unwrap();
    let mut lines = Vec::new();
    while let Some(rec) = reader.next_record().unwrap() {
        lines.push(rec.line().to_owned());
    }
    lines
}

#[test]
fn reads_bcftools_bcf_records() {
    let bytes = std::fs::read(fixture("oracle.bcf")).unwrap();
    assert_eq!(bcf_lines(bytes), expected("expected_bcf_records.txt"));
}

#[test]
fn our_bcf_writer_output_reads_back_as_bcftools_records() {
    let mut reader =
        VcfReader::new(BufReader::new(File::open(fixture("in.vcf")).unwrap())).unwrap();
    let header = BcfHeader::from_vcf(reader.header().clone()).unwrap();
    let mut writer = BcfWriter::new_bgzf(Vec::new(), &header).unwrap();
    while let Some(rec) = reader.next_record().unwrap() {
        writer.write_record(&rec).unwrap();
    }
    let bytes = writer.finish().unwrap();
    assert_eq!(bcf_lines(bytes), expected("expected_bcf_records.txt"));
}

#[test]
fn tbi_built_for_vcf_gz_matches_tabix_bytes() {
    let built = TabixIndex::build_text(
        BufReader::new(File::open(fixture("oracle.vcf.gz")).unwrap()),
        TabixConfig::vcf(),
        IndexKind::Tbi,
    )
    .unwrap();
    let mut ours = Vec::new();
    built.write_to(&mut ours).unwrap();
    let theirs = std::fs::read(fixture("oracle.vcf.gz.tbi")).unwrap();
    assert_eq!(gunzip(&ours), gunzip(&theirs));
}

#[test]
fn csi_built_for_bcf_matches_bcftools_bytes() {
    let kind = IndexKind::Csi { min_shift: 14, depth: 0 };
    let built =
        build_bcf_index(BufReader::new(File::open(fixture("oracle.bcf")).unwrap()), kind).unwrap();
    let mut ours = Vec::new();
    built.write_to(&mut ours).unwrap();
    let theirs = std::fs::read(fixture("oracle.bcf.csi")).unwrap();
    assert_eq!(gunzip(&ours), gunzip(&theirs));
}

#[test]
fn tbi_region_query_matches_tabix() {
    let index =
        TabixIndex::read_from(BufReader::new(File::open(fixture("oracle.vcf.gz.tbi")).unwrap()))
            .unwrap();
    let mut bgzf = BgzfReader::new(File::open(fixture("oracle.vcf.gz")).unwrap());
    let got = fetch_text(&mut bgzf, &index, "chr1", 1000, 3000).unwrap();
    assert_eq!(got, expected("expected_vcfgz_chr1_1000_3000.txt"));
}

#[test]
fn csi_region_query_matches_bcftools() {
    let index =
        TabixIndex::read_from(BufReader::new(File::open(fixture("oracle.bcf.csi")).unwrap()))
            .unwrap();
    let reader =
        BcfReader::new_bgzf(BufReader::new(File::open(fixture("oracle.bcf")).unwrap())).unwrap();
    let header = reader.header().clone();
    let mut bgzf = BgzfReader::new(File::open(fixture("oracle.bcf")).unwrap());
    let got = fetch_lines(&mut bgzf, &header, &index, "chr2", 100, 200).unwrap();
    assert_eq!(got, expected("expected_bcf_chr2_100_200.txt"));
}
