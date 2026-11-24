//! Sequence quickstart: write a sorted BAM, build its BAI index, run a region query
//! and export the records as CSV. Run with:
//!
//! ```text
//! cargo run -p tpt-bio-seq --example seq_quickstart [output-dir]
//! ```
//!
//! The input is generated here, so no data files are needed.

use std::fs::{self, File};
use std::io::{BufReader, BufWriter};
use std::path::{Path, PathBuf};

use tpt_bio_seq::bai::{self, Index};
use tpt_bio_seq::table::records_to_csv;
use tpt_bio_seq::{BamReader, BamWriter, CigarKind, CigarOp, Header, Record, Reference};

type Failure = Box<dyn std::error::Error>;

/// Builds a coordinate-sorted set of 50-base reads on two references.
fn sample_records() -> Vec<Record> {
    // A small deterministic generator, so the output is the same on every run.
    let mut state: u64 = 7;
    let mut next = move || {
        state =
            state.wrapping_mul(6_364_136_223_846_793_005).wrapping_add(1_442_695_040_888_963_407);
        (state >> 33) as i64
    };
    let mut records: Vec<Record> = (0..500)
        .map(|i| Record {
            name: format!("read{i:04}").into_bytes(),
            ref_id: 0,
            pos: next().rem_euclid(99_000),
            mapq: 60,
            cigar: vec![CigarOp { len: 50, kind: CigarKind::Match }],
            seq: b"ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTAC".to_vec(),
            qual: vec![40; 50],
            ..Record::default()
        })
        .collect();
    records.push(Record {
        name: b"read9999".to_vec(),
        ref_id: 1,
        pos: 299,
        mapq: 30,
        cigar: vec![CigarOp { len: 10, kind: CigarKind::Match }],
        seq: b"GGCCGGCCGG".to_vec(),
        qual: vec![40; 10],
        ..Record::default()
    });
    records.sort_by_key(|r| (r.ref_id, r.pos));
    records
}

fn write_bam(path: &Path, header: &Header, records: &[Record]) -> Result<(), Failure> {
    let mut writer = BamWriter::new(BufWriter::new(File::create(path)?), header.clone())?;
    for rec in records {
        writer.write_record(rec)?;
    }
    writer.finish()?;
    Ok(())
}

fn open_bam(path: &Path) -> Result<BamReader<BufReader<File>>, Failure> {
    Ok(BamReader::new(BufReader::new(File::open(path)?))?)
}

fn main() -> Result<(), Failure> {
    let dir: PathBuf = std::env::args_os()
        .nth(1)
        .map_or_else(|| std::env::temp_dir().join("tptbio-seq-quickstart"), PathBuf::from);
    fs::create_dir_all(&dir)?;

    // 1. Write a sorted BAM with a two-reference header.
    let header = Header::from_references(&[
        Reference { name: "chr1".into(), len: 100_000 },
        Reference { name: "chr2".into(), len: 50_000 },
    ])?;
    let records = sample_records();
    let bam_path = dir.join("reads.bam");
    write_bam(&bam_path, &header, &records)?;

    // 2. Build the BAI index from the BAM and write it next to it.
    let index = bai::build_index(&mut open_bam(&bam_path)?)?;
    index.write_to(BufWriter::new(File::create(dir.join("reads.bam.bai"))?))?;

    // 3. Query chr1:1000-2000 (1-based, inclusive) through the index.
    // 0-based half-open coordinates: chr1:1000-2000 is [999, 2000).
    let index = Index::read_from(BufReader::new(File::open(dir.join("reads.bam.bai"))?))?;
    let mut reader = open_bam(&bam_path)?;
    let hits = bai::fetch(&mut reader, &index, 0, 999, 2000)?;
    println!("{} records overlap chr1:1000-2000", hits.len());

    // 4. Export the reads that were found as CSV.
    let csv = records_to_csv(&hits, reader.header())?;
    let csv_path = dir.join("region.csv");
    fs::write(&csv_path, csv)?;
    println!("wrote {}", csv_path.display());
    Ok(())
}
