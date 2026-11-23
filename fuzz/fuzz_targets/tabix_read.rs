#![no_main]

use libfuzzer_sys::fuzz_target;

fuzz_target!(|data: &[u8]| {
    if let Ok(index) = tpt_bio_seq::TabixIndex::read_from(data) {
        let _ = index.names();
        if index.n_refs() > 0 {
            let _ = index.query_by_id(0, 0, 1 << 20);
        }
    }
});
