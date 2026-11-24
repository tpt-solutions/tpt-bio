//! Checks BAM and BAI against files written by samtools 1.22.
//!
//! `tests/data/samtools/sorted.bam` is the 250-record alignment set in coordinate order,
//! written by `samtools view -b`. `samtools_sorted.bai` is `samtools index` output for it.
//! `expected_bam_records.txt` is `samtools view` of the file, sorted, and
//! `expected_bam_region_names.txt` lists the read names that `samtools view chr1:1000-3000`
//! returns. The tests need no samtools installation.

use std::fs::File;
use std::io::{BufReader, Cursor};
use std::path::PathBuf;

use tpt_bio_seq::bai::{Index, build_index, fetch};
use tpt_bio_seq::{BamReader, SamWriter};

fn fixture(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/data/samtools").join(name)
}

fn expected_lines(name: &str) -> Vec<String> {
    std::fs::read_to_string(fixture(name)).unwrap().lines().map(str::to_owned).collect()
}

/// Sorts the tags of a SAM line, so the comparison does not depend on tag order.
fn sort_tags(line: &str) -> String {
    let fields: Vec<&str> = line.split('\t').collect();
    let mut tags: Vec<&str> = fields[11..].to_vec();
    tags.sort_unstable();
    let mut out: Vec<&str> = fields[..11].to_vec();
    out.extend(tags);
    out.join("\t")
}

#[test]
fn bai_built_from_samtools_bam_matches_samtools_index_bytes() {
    let mut reader =
        BamReader::new(BufReader::new(File::open(fixture("sorted.bam")).unwrap())).unwrap();
    let index = build_index(&mut reader).unwrap();
    let mut ours = Vec::new();
    index.write_to(&mut ours).unwrap();
    assert_eq!(ours, std::fs::read(fixture("samtools_sorted.bai")).unwrap());
}

#[test]
fn decodes_samtools_bam_to_samtools_records() {
    let mut reader =
        BamReader::new(BufReader::new(File::open(fixture("sorted.bam")).unwrap())).unwrap();
    let mut sam = SamWriter::new(Vec::new(), reader.header().clone()).unwrap();
    while let Some(record) = reader.next_record().unwrap() {
        sam.write_record(&record).unwrap();
    }
    let text = String::from_utf8(sam.finish().unwrap()).unwrap();
    let mut got: Vec<String> =
        text.lines().filter(|l| !l.starts_with('@')).map(sort_tags).collect();
    got.sort();
    let mut want: Vec<String> =
        expected_lines("expected_bam_records.txt").iter().map(|l| sort_tags(l)).collect();
    want.sort();
    assert_eq!(got, want);
}

#[test]
fn region_query_with_samtools_index_matches_samtools_view() {
    let index = Index::read_from(File::open(fixture("samtools_sorted.bai")).unwrap()).unwrap();
    let bam = std::fs::read(fixture("sorted.bam")).unwrap();
    let mut reader = BamReader::new(Cursor::new(bam)).unwrap();
    // `chr1:1000-3000` in 1-based inclusive coordinates is [999, 3000) in 0-based.
    let records = fetch(&mut reader, &index, 0, 999, 3000).unwrap();
    let mut names: Vec<String> =
        records.iter().map(|r| String::from_utf8_lossy(&r.name).into_owned()).collect();
    names.sort();
    assert_eq!(names, expected_lines("expected_bam_region_names.txt"));
}
