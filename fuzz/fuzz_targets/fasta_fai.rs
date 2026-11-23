#![no_main]

use libfuzzer_sys::fuzz_target;

fuzz_target!(|data: &[u8]| {
    // Streaming FASTA reader.
    let mut reader = tpt_bio_seq::fasta::FastaReader::new(data);
    while let Ok(Some(_)) = reader.next_record() {}

    // Index built from the FASTA bytes, then used for a bounded fetch.
    if let Ok(fai) = tpt_bio_seq::Fai::build(data) {
        for entry in fai.entries() {
            let end = entry.length.min(1 << 16);
            let _ = fai.fetch(data, &entry.name, 0, end);
        }
    }

    // Untrusted .fai text.
    if let Ok(text) = std::str::from_utf8(data) {
        let _ = tpt_bio_seq::Fai::read_from(text);
    }
});
