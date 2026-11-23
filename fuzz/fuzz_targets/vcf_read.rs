#![no_main]

use libfuzzer_sys::fuzz_target;

fuzz_target!(|data: &[u8]| {
    if let Ok(mut reader) = tpt_bio_seq::VcfReader::new(data) {
        let _ = reader.header().contigs().count();
        while let Ok(Some(rec)) = reader.next_record() {
            let _ = rec.pos();
            let _ = rec.alts().count();
            let _ = rec.info().count();
            let _ = rec.samples().count();
        }
    }
});
