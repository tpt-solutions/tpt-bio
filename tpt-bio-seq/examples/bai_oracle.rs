//! Builds a BAI index for a coordinate-sorted BAM with `tpt-bio-seq`, for byte-level
//! comparison with the index `samtools index` writes.
//!
//! ```text
//! cargo run -p tpt-bio-seq --example bai_oracle -- in.bam out.bai
//! ```

use std::fs::File;
use std::io::{BufReader, BufWriter};
use std::process::ExitCode;

use tpt_bio_seq::{BamReader, bai::build_index};

fn main() -> ExitCode {
    let args: Vec<String> = std::env::args().collect();
    let [_, bam_path, bai_path] = args.as_slice() else {
        eprintln!("usage: bai_oracle <in.bam> <out.bai>");
        return ExitCode::from(2);
    };
    let result = (|| -> Result<(), Box<dyn std::error::Error>> {
        let mut reader = BamReader::new(BufReader::new(File::open(bam_path)?))?;
        let index = build_index(&mut reader)?;
        index.write_to(BufWriter::new(File::create(bai_path)?))?;
        Ok(())
    })();
    match result {
        Ok(()) => ExitCode::SUCCESS,
        Err(e) => {
            eprintln!("error: {e}");
            ExitCode::FAILURE
        }
    }
}
