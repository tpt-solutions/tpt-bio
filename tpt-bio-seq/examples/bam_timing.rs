//! Times BAM decoding and SAM formatting separately, for comparing with samtools.
//!
//! ```text
//! cargo run --release -p tpt-bio-seq --example bam_timing -- in.bam
//! ```

use std::fs::File;
use std::io::BufReader;
use std::time::Instant;

use tpt_bio_seq::{BamReader, SamWriter};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let path = std::env::args().nth(1).ok_or("usage: bam_timing <in.bam>")?;

    let start = Instant::now();
    let mut viewer = BamReader::new(BufReader::with_capacity(1 << 16, File::open(&path)?))?;
    let mut n_views = 0u64;
    while viewer.next_view()?.is_some() {
        n_views += 1;
    }
    let view_only = start.elapsed();

    let start = Instant::now();
    let mut reader = BamReader::new(BufReader::with_capacity(1 << 16, File::open(&path)?))?;
    let mut records = Vec::new();
    while let Some(rec) = reader.next_record()? {
        records.push(rec);
    }
    let decode = start.elapsed();

    let header = reader.header().clone();
    let start = Instant::now();
    let mut sam = SamWriter::new(std::io::sink(), header)?;
    for rec in &records {
        sam.write_record(rec)?;
    }
    sam.finish()?;
    let format = start.elapsed();

    println!("records {n_views}");
    println!("views   {:.1} ms (framing and decompression only)", view_only.as_secs_f64() * 1e3);
    println!("decode  {:.1} ms", decode.as_secs_f64() * 1e3);
    println!("format  {:.1} ms", format.as_secs_f64() * 1e3);
    Ok(())
}
