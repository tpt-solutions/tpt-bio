#![no_main]

use libfuzzer_sys::fuzz_target;

fuzz_target!(|data: &[u8]| {
    for record in tpt_bio_seq::FastqParser::new(data) {
        if record.is_err() {
            break;
        }
    }
});
