//! Times CRAM record decoding and SAM formatting separately.
//!
//! ```text
//! cargo run --release -p tpt-bio-seq --example cram_timing -- in.cram
//! ```

use std::fs::File;
use std::io::BufReader;
use std::time::Instant;

use tpt_bio_seq::{CramReader, SamWriter};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let path = std::env::args().nth(1).ok_or("usage: cram_timing <in.cram>")?;
    let start = Instant::now();
    let mut reader = CramReader::new(BufReader::with_capacity(1 << 16, File::open(&path)?))?;
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
    println!("records {}", records.len());
    println!("decode  {:.1} ms", decode.as_secs_f64() * 1e3);
    println!("format  {:.1} ms", start.elapsed().as_secs_f64() * 1e3);
    Ok(())
}
