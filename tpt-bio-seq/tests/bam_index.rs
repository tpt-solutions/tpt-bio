//! End-to-end: write a coordinate-sorted BAM, build its BAI, and fetch regions.

use std::fmt::Write as _;
use std::io::Cursor;

use tpt_bio_seq::bai::{build_index, fetch};
use tpt_bio_seq::sam::SamReader;
use tpt_bio_seq::{BamReader, BamWriter};

/// Deterministic pseudo-random sequence so the test needs no fixtures.
fn seq_for(i: u64, len: usize) -> String {
    let bases = b"ACGT";
    let mut x = i.wrapping_mul(0x9E37_79B9_7F4A_7C15) | 1;
    (0..len)
        .map(|_| {
            x ^= x << 13;
            x ^= x >> 7;
            x ^= x << 17;
            bases[(x % 4) as usize] as char
        })
        .collect()
}

fn build_sam(n: u64) -> String {
    let mut s = String::from(
        "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:2000000\n@SQ\tSN:chr2\tLN:500000\n",
    );
    // chr1: a read every 500 bp, 100 bp each, so the 16 KiB linear index is exercised.
    for i in 0..n {
        let pos = 1 + i * 500;
        let seq = seq_for(i, 100);
        let qual = "I".repeat(100);
        let _ = writeln!(s, "read{i}\t0\tchr1\t{pos}\t60\t100M\t*\t0\t0\t{seq}\t{qual}\tNM:i:0");
    }
    for i in 0..10u64 {
        let pos = 1 + i * 1000;
        let seq = seq_for(i + 9999, 50);
        let _ =
            writeln!(s, "c2_{i}\t0\tchr2\t{pos}\t60\t50M\t*\t0\t0\t{seq}\t{}\n", "I".repeat(50));
    }
    s
}

fn write_bam(sam: &str) -> Vec<u8> {
    let mut reader = SamReader::new(sam.as_bytes()).unwrap();
    let header = reader.header().clone();
    let mut out = Vec::new();
    let mut w = BamWriter::new(&mut out, header).unwrap();
    while let Some(rec) = reader.next_record().unwrap() {
        w.write_record(&rec).unwrap();
    }
    w.finish().unwrap();
    out
}

fn brute_force(sam: &str, chrom: &str, beg: u64, end: u64) -> Vec<String> {
    let mut reader = SamReader::new(sam.as_bytes()).unwrap();
    let mut names = Vec::new();
    while let Some(rec) = reader.next_record().unwrap() {
        let rid = reader.header().ref_id(chrom).unwrap();
        if rec.ref_id != rid {
            continue;
        }
        let rstart = rec.pos as u64;
        let rend = rstart + tpt_bio_seq::record::reference_length(&rec.cigar);
        if rstart < end && rend > beg {
            names.push(String::from_utf8(rec.name.clone()).unwrap());
        }
    }
    names
}

#[test]
fn index_fetch_matches_brute_force() {
    let sam = build_sam(2000);
    let bam = write_bam(&sam);

    let mut reader = BamReader::new(Cursor::new(&bam)).unwrap();
    let index = build_index(&mut reader).unwrap();
    assert_eq!(index.n_refs(), 2);

    let mut indexed = BamReader::new(Cursor::new(&bam)).unwrap();
    let chr1 = indexed.header().ref_id("chr1").unwrap() as usize;
    let chr2 = indexed.header().ref_id("chr2").unwrap() as usize;

    // (ref id, name, beg, end, expect_reads): the last window lies past the data.
    let cases = [
        (chr1, "chr1", 0u64, 1_000u64, true),
        (chr1, "chr1", 250_000, 251_000, true),
        (chr1, "chr1", 999_400, 999_700, true),
        (chr1, "chr1", 1_900_000, 2_000_000, false),
        (chr2, "chr2", 4_000, 6_500, true),
    ];
    for (rid, name, beg, end, expect_reads) in cases {
        let mut got: Vec<String> = fetch(&mut indexed, &index, rid, beg, end)
            .unwrap()
            .into_iter()
            .map(|r| String::from_utf8(r.name).unwrap())
            .collect();
        got.sort();
        let mut expect = brute_force(&sam, name, beg, end);
        expect.sort();
        assert_eq!(got, expect, "region {name}:{beg}-{end}");
        assert_eq!(!got.is_empty(), expect_reads, "region {name}:{beg}-{end}");
    }
}

#[test]
fn unsorted_input_is_rejected() {
    let sam = "@SQ\tSN:chr1\tLN:1000\nr1\t0\tchr1\t500\t60\t10M\t*\t0\t0\tACGTACGTAC\t*\nr2\t0\tchr1\t100\t60\t10M\t*\t0\t0\tACGTACGTAC\t*\n";
    let bam = write_bam(sam);
    let mut reader = BamReader::new(Cursor::new(&bam)).unwrap();
    assert!(build_index(&mut reader).is_err());
}
