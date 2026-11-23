#![no_main]

use libfuzzer_sys::fuzz_target;

fuzz_target!(|data: &[u8]| {
    // Plain BCF, then the usual BGZF-wrapped form.
    if let Ok(mut reader) = tpt_bio_seq::BcfReader::new(data) {
        while let Ok(Some(rec)) = reader.next_record() {
            let _ = rec.pos();
            let _ = rec.samples().count();
        }
    }
    if let Ok(mut reader) = tpt_bio_seq::BcfReader::new_bgzf(data) {
        while let Ok(Some(rec)) = reader.next_record() {
            let _ = rec.pos();
            let _ = rec.samples().count();
        }
    }
});
