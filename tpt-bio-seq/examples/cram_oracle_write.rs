//! Converts a SAM file to CRAM 3.0 with `tpt-bio-seq`, for checking with `samtools`.
//!
//! ```text
//! cargo run -p tpt-bio-seq --example cram_oracle_write -- in.sam ref.fa out.cram
//! ```

use std::fs::File;
use std::io::{BufReader, BufWriter};
use std::process::ExitCode;

use tpt_bio_seq::{CramReference, CramWriter, SamReader};

fn main() -> ExitCode {
    let args: Vec<String> = std::env::args().collect();
    let [_, sam_path, ref_path, out_path] = args.as_slice() else {
        eprintln!("usage: cram_oracle_write <in.sam> <ref.fa> <out.cram>");
        return ExitCode::from(2);
    };
    match run(sam_path, ref_path, out_path) {
        Ok(n) => {
            println!("wrote {n} records to {out_path}");
            ExitCode::SUCCESS
        }
        Err(e) => {
            eprintln!("error: {e}");
            ExitCode::FAILURE
        }
    }
}

fn run(sam_path: &str, ref_path: &str, out_path: &str) -> Result<u64, Box<dyn std::error::Error>> {
    let mut reader = SamReader::new(BufReader::new(File::open(sam_path)?))?;
    let header = reader.header().clone();
    let reference = CramReference::from_path(ref_path)?;
    let out = BufWriter::new(File::create(out_path)?);
    let mut writer = CramWriter::new(out, header, Some(reference))?;
    let mut count = 0u64;
    while let Some(record) = reader.next_record()? {
        writer.write_record(&record)?;
        count += 1;
    }
    writer.finish()?;
    Ok(count)
}
