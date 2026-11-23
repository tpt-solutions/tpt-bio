#![no_main]

use libfuzzer_sys::fuzz_target;

fuzz_target!(|data: &[u8]| {
    let _ = tpt_bio_seq::Index::read_from(data);
});
