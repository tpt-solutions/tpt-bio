//! Converts a SAM file to BAM with `tpt-bio-seq`, for checking with `samtools`.
//!
//! ```text
//! cargo run -p tpt-bio-seq --example bam_oracle_write -- in.sam out.bam
//! ```

use std::fs::File;
use std::io::{BufReader, BufWriter};
use std::process::ExitCode;

use tpt_bio_seq::{BamWriter, SamReader};

fn main() -> ExitCode {
    let args: Vec<String> = std::env::args().collect();
    let [_, sam_path, bam_path] = args.as_slice() else {
        eprintln!("usage: bam_oracle_write <in.sam> <out.bam>");
        return ExitCode::from(2);
    };
    let result = (|| -> Result<u64, Box<dyn std::error::Error>> {
        let mut reader = SamReader::new(BufReader::new(File::open(sam_path)?))?;
        let mut writer =
            BamWriter::new(BufWriter::new(File::create(bam_path)?), reader.header().clone())?;
        let mut count = 0u64;
        while let Some(record) = reader.next_record()? {
            writer.write_record(&record)?;
            count += 1;
        }
        writer.finish()?;
        Ok(count)
    })();
    match result {
        Ok(n) => {
            println!("wrote {n} records to {bam_path}");
            ExitCode::SUCCESS
        }
        Err(e) => {
            eprintln!("error: {e}");
            ExitCode::FAILURE
        }
    }
}
