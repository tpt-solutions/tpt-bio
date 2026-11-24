//! Times BCF record decoding (to VCF text) separately from writing the text.
//!
//! ```text
//! cargo run --release -p tpt-bio-seq --example bcf_timing -- in.bcf
//! ```

use std::fs::File;
use std::io::{BufReader, Write};
use std::time::Instant;

use tpt_bio_seq::BcfReader;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let path = std::env::args().nth(1).ok_or("usage: bcf_timing <in.bcf>")?;
    let start = Instant::now();
    let mut reader = BcfReader::new_bgzf(BufReader::with_capacity(1 << 16, File::open(&path)?))?;
    let mut n = 0u64;
    let mut bytes = 0u64;
    let mut sink = std::io::sink();
    while let Some(rec) = reader.next_record()? {
        n += 1;
        bytes += rec.line().len() as u64;
        sink.write_all(rec.line().as_bytes())?;
    }
    println!("records {n} bytes {bytes} total {:.1} ms", start.elapsed().as_secs_f64() * 1e3);
    Ok(())
}
