//! Decodes a CRAM file with `tpt-bio-seq` and prints it as SAM on stdout.
//!
//! Used to compare against `samtools view` output when checking CRAM conformance.
//! The reference FASTA is needed for reference-based files. Pass `-` in its place to
//! open the file without a reference.
//!
//! ```text
//! cargo run -p tpt-bio-seq --example cram_oracle_dump -- in.cram ref.fa > ours.sam
//! cargo run -p tpt-bio-seq --example cram_oracle_dump -- in.cram - > ours.sam
//! ```

use std::fs::File;
use std::io::{BufReader, BufWriter, Write};
use std::process::ExitCode;

use tpt_bio_seq::{CramReader, CramReference, SamWriter};

fn main() -> ExitCode {
    let args: Vec<String> = std::env::args().collect();
    let [_, cram_path, ref_path] = args.as_slice() else {
        eprintln!("usage: cram_oracle_dump <in.cram> <ref.fa>");
        return ExitCode::from(2);
    };
    match run(cram_path, ref_path) {
        Ok(()) => ExitCode::SUCCESS,
        Err(e) => {
            eprintln!("error: {e}");
            ExitCode::FAILURE
        }
    }
}

fn run(cram_path: &str, ref_path: &str) -> Result<(), Box<dyn std::error::Error>> {
    let file = BufReader::new(File::open(cram_path)?);
    let mut reader = if ref_path == "-" {
        CramReader::new(file)?
    } else {
        CramReader::with_reference(file, CramReference::from_path(ref_path)?)?
    };
    let header = reader.header().clone();
    let mut out = SamWriter::new(BufWriter::new(std::io::stdout().lock()), header)?;
    while let Some(record) = reader.next_record()? {
        out.write_record(&record)?;
    }
    let mut inner = out.finish()?;
    inner.flush()?;
    Ok(())
}
