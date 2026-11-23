#![no_main]

use libfuzzer_sys::fuzz_target;

fuzz_target!(|data: &[u8]| {
    if let Ok(mut reader) = tpt_bio_seq::SamReader::new(data) {
        while let Ok(Some(_)) = reader.next_record() {}
    }
});
