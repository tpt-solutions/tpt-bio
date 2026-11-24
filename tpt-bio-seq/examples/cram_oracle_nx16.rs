//! Decodes one rANS Nx16 block payload and writes the raw bytes.
//!
//! Used to compare block-level output against the htscodecs library that samtools
//! links to. The payload is the block data as stored in the file, starting with the
//! flags byte. The raw size is the block's uncompressed size.
//!
//! ```text
//! cargo run -p tpt-bio-seq --example cram_oracle_nx16 -- block.bin 85398 ours.bin
//! ```

use std::process::ExitCode;

use tpt_bio_seq::cram::rans_nx16;

fn main() -> ExitCode {
    let args: Vec<String> = std::env::args().collect();
    let [_, input, raw_size, output] = args.as_slice() else {
        eprintln!("usage: cram_oracle_nx16 <block.bin> <raw size> <out.bin>");
        return ExitCode::from(2);
    };
    let Ok(raw_size) = raw_size.parse::<usize>() else {
        eprintln!("raw size must be an integer");
        return ExitCode::from(2);
    };
    let payload = match std::fs::read(input) {
        Ok(p) => p,
        Err(e) => {
            eprintln!("error: {input}: {e}");
            return ExitCode::FAILURE;
        }
    };
    match rans_nx16::decode(&payload, raw_size) {
        Ok(bytes) => match std::fs::write(output, bytes) {
            Ok(()) => ExitCode::SUCCESS,
            Err(e) => {
                eprintln!("error: {output}: {e}");
                ExitCode::FAILURE
            }
        },
        Err(e) => {
            eprintln!("error: {e}");
            ExitCode::FAILURE
        }
    }
}
