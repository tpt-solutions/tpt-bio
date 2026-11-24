//! Builds a TBI index for a bgzipped VCF, and a CSI index for a BCF file, with
//! `tpt-bio-seq`, for byte-level comparison with the indexes `tabix` and `bcftools` write.
//!
//! ```text
//! cargo run -p tpt-bio-seq --example index_oracle -- in.vcf.gz out.tbi in.bcf out.csi
//! ```

use std::fs::File;
use std::io::{BufReader, BufWriter};
use std::process::ExitCode;

use tpt_bio_seq::{IndexKind, TabixConfig, TabixIndex, bcf::build_bcf_index};

fn main() -> ExitCode {
    let args: Vec<String> = std::env::args().collect();
    let [_, vcf_gz, tbi_out, bcf, csi_out] = args.as_slice() else {
        eprintln!("usage: index_oracle <in.vcf.gz> <out.tbi> <in.bcf> <out.csi>");
        return ExitCode::from(2);
    };
    match run(vcf_gz, tbi_out, bcf, csi_out) {
        Ok(()) => ExitCode::SUCCESS,
        Err(e) => {
            eprintln!("error: {e}");
            ExitCode::FAILURE
        }
    }
}

fn run(
    vcf_gz: &str,
    tbi_out: &str,
    bcf: &str,
    csi_out: &str,
) -> Result<(), Box<dyn std::error::Error>> {
    let tbi = TabixIndex::build_text(
        BufReader::new(File::open(vcf_gz)?),
        TabixConfig::vcf(),
        IndexKind::Tbi,
    )?;
    tbi.write_to(BufWriter::new(File::create(tbi_out)?))?;
    // `-` skips the CSI step. Depth 0 only covers contigs up to 16 kb.
    if csi_out != "-" {
        let csi = build_bcf_index(
            BufReader::new(File::open(bcf)?),
            IndexKind::Csi { min_shift: 14, depth: 0 },
        )?;
        csi.write_to(BufWriter::new(File::create(csi_out)?))?;
    }
    Ok(())
}
