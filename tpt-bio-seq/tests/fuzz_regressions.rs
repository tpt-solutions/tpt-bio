//! Regression tests for inputs that fuzzing found to exhaust memory.
//!
//! Each test reads a minimised input from `tests/data/fuzz_regressions/` and asserts that
//! the parser returns (with `Ok` or `Err`) without panicking.
//!
//! Limitation: these tests do not detect the allocation problem. The inputs declare an
//! ISIZE of about 4 GB, and an ordinary test run reserves that much address space without
//! failing (Windows and Linux overcommit), so the tests pass on the unfixed code. The
//! fuzzer detects the problem through libFuzzer's allocation and RSS limits. Keep the
//! tests as cheap guards against panics; a fix needs its own unit test, for example one
//! that asserts `bgzf` rejects ISIZE above `MAX_BLOCK_SIZE`.

use std::path::PathBuf;

fn regression_input(name: &str) -> Vec<u8> {
    let path: PathBuf =
        [env!("CARGO_MANIFEST_DIR"), "tests", "data", "fuzz_regressions", name].iter().collect();
    std::fs::read(&path).unwrap_or_else(|e| panic!("cannot read {}: {e}", path.display()))
}

/// Found by the `bgzf_decompress` fuzz target on 2026-11-24 (Windows, MSVC ASan):
/// the footer's ISIZE field (0xBE4074F0, about 3 GB) reaches `Vec::with_capacity` in
/// `bgzf::inflate_block` before the block is decompressed or checked against
/// `MAX_BLOCK_SIZE`. libFuzzer reports `out-of-memory (malloc(3191895280))`.
#[test]
fn bgzf_isize_oom_returns_without_panic() {
    let data = regression_input("bgzf_decompress_isize_oom.bin");
    let _ = tpt_bio_seq::bgzf::decompress_parallel(&data);
}

/// Found by the `bcf_read` fuzz target on 2026-11-24: the BGZF layer under
/// `BcfReader::new_bgzf` hits the same unbounded ISIZE allocation (`bgzf::inflate_block`,
/// reached from `bgzf::BgzfReader::load_next_block`). libFuzzer reports
/// `out-of-memory (malloc(4294902015))`.
#[test]
fn bcf_bgzf_isize_oom_returns_without_panic() {
    let data = regression_input("bcf_read_isize_oom.bin");
    if let Ok(mut reader) = tpt_bio_seq::BcfReader::new_bgzf(data.as_slice()) {
        while let Ok(Some(rec)) = reader.next_record() {
            let _ = rec.pos();
            let _ = rec.samples().count();
        }
    }
}

/// Found by the `bam_read` fuzz target on 2026-11-24: a BAM stream whose first BGZF block
/// declares an ISIZE of about 4 GB. It fails in the same `bgzf::inflate_block` allocation
/// (`malloc(4294967295)`).
#[test]
fn bam_bgzf_isize_oom_returns_without_panic() {
    let data = regression_input("bam_read_isize_oom.bin");
    if let Ok(mut reader) = tpt_bio_seq::BamReader::new(data.as_slice()) {
        while let Ok(Some(_)) = reader.next_record() {}
    }
}
