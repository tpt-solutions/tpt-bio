#![no_main]

use libfuzzer_sys::fuzz_target;
use tpt_bio_seq::cram::{block, bytes::Cursor, rans, rans_nx16};

fuzz_target!(|data: &[u8]| {
    // Whole-file CRAM reading without a reference.
    if let Ok(mut reader) = tpt_bio_seq::CramReader::new(data) {
        while let Ok(Some(_)) = reader.next_record() {}
    }

    // Entropy codecs directly, on raw bytes.
    let _ = rans::decode(data);
    let _ = rans_nx16::decode(data, data.len());

    // Block framing with every method code the input can select.
    let mut cur = Cursor::new(data);
    for _ in 0..64 {
        if cur.is_empty() || block::read_block(&mut cur).is_err() {
            break;
        }
    }
    for method in [block::Method::Rans4x8, block::Method::RansNx16] {
        let _ = block::decompress(method, data, data.len() * 4);
    }
});
