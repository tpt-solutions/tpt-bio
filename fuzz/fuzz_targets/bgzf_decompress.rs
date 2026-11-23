#![no_main]

use libfuzzer_sys::fuzz_target;

fuzz_target!(|data: &[u8]| {
    // Must return an error, never panic or allocate without bound.
    let _ = tpt_bio_seq::bgzf::decompress_parallel(data);
});
