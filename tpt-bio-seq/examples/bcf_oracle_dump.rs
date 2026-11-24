//! Decodes a BCF file with `tpt-bio-seq` and prints it as VCF text on stdout.
//!
//! Used to compare against `bcftools view` output when checking BCF conformance.
//!
//! ```text
//! cargo run -p tpt-bio-seq --example bcf_oracle_dump -- in.bcf > ours.vcf
//! ```

use std::fs::File;
use std::io::{BufReader, BufWriter, Write};
use std::process::ExitCode;

use tpt_bio_seq::BcfReader;

fn main() -> ExitCode {
    let args: Vec<String> = std::env::args().collect();
    let [_, path] = args.as_slice() else {
        eprintln!("usage: bcf_oracle_dump <in.bcf>");
        return ExitCode::from(2);
    };
    match run(path) {
        Ok(()) => ExitCode::SUCCESS,
        Err(e) => {
            eprintln!("error: {e}");
            ExitCode::FAILURE
        }
    }
}

fn run(path: &str) -> Result<(), Box<dyn std::error::Error>> {
    let mut reader = BcfReader::new_bgzf(BufReader::new(File::open(path)?))?;
    let mut out = BufWriter::new(std::io::stdout().lock());
    out.write_all(reader.header().vcf().text().as_bytes())?;
    while let Some(record) = reader.next_record()? {
        writeln!(out, "{}", record.line())?;
    }
    out.flush()?;
    Ok(())
}
