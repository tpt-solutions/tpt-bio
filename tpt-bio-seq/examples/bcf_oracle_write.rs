//! Converts a VCF file to BCF with `tpt-bio-seq`, for checking with `bcftools`.
//!
//! ```text
//! cargo run -p tpt-bio-seq --example bcf_oracle_write -- in.vcf out.bcf
//! ```

use std::fs::File;
use std::io::{BufReader, BufWriter};
use std::process::ExitCode;

use tpt_bio_seq::{BcfHeader, BcfWriter, VcfReader};

fn main() -> ExitCode {
    let args: Vec<String> = std::env::args().collect();
    let [_, in_path, out_path] = args.as_slice() else {
        eprintln!("usage: bcf_oracle_write <in.vcf> <out.bcf>");
        return ExitCode::from(2);
    };
    match run(in_path, out_path) {
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

fn run(in_path: &str, out_path: &str) -> Result<u64, Box<dyn std::error::Error>> {
    let mut reader = VcfReader::new(BufReader::new(File::open(in_path)?))?;
    let header = BcfHeader::from_vcf(reader.header().clone())?;
    let out = BufWriter::new(File::create(out_path)?);
    let mut writer = BcfWriter::new_bgzf(out, &header)?;
    let mut count = 0u64;
    while let Some(record) = reader.next_record()? {
        writer.write_record(&record)?;
        count += 1;
    }
    writer.finish()?;
    Ok(count)
}
