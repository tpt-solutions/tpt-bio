//! Times BGZF decompression of a whole file, for comparing with raw inflate cost.
//!
//! ```text
//! cargo run --release -p tpt-bio-seq --example bgzf_timing -- in.bam
//! ```

use std::fs::File;
use std::io::{BufReader, Read};
use std::time::Instant;

use tpt_bio_seq::BgzfReader;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let path = std::env::args().nth(1).ok_or("usage: bgzf_timing <in.bgz>")?;
    let start = Instant::now();
    let mut reader = BgzfReader::new(BufReader::with_capacity(1 << 16, File::open(&path)?));
    let mut out = Vec::new();
    reader.read_to_end(&mut out)?;
    println!("bytes out {} in {:.1} ms", out.len(), start.elapsed().as_secs_f64() * 1e3);
    Ok(())
}
