//! Property tests for the nucleotide, k-mer and coordinate APIs.

use proptest::prelude::*;
use tpt_bio_core::coords::{Interval, Region};
use tpt_bio_core::kmer::{self, Kmers};
use tpt_bio_core::nucleotide::{Base, Iupac, reverse_complement, reverse_complement_in_place};

// U is excluded: it is read as T on input, so it does not round-trip byte for byte.
const IUPAC_BYTES: &[u8] = b"ACGTNRYKMSWBDHVacgtnrykmswbdhv-";
const ACGT: &[u8] = b"ACGT";

fn iupac_seq(max_len: usize) -> impl Strategy<Value = Vec<u8>> {
    prop::collection::vec(prop::sample::select(IUPAC_BYTES), 0..max_len)
}

fn acgt_seq(min_len: usize, max_len: usize) -> impl Strategy<Value = Vec<u8>> {
    prop::collection::vec(prop::sample::select(ACGT), min_len..max_len)
}

proptest! {
    #[test]
    fn reverse_complement_is_an_involution(seq in iupac_seq(256)) {
        // Bytes outside the IUPAC alphabet are rejected, so only test valid input.
        if let Ok(rc) = reverse_complement(&seq) {
            prop_assert_eq!(reverse_complement(&rc).unwrap(), seq);
        }
    }

    #[test]
    fn reverse_complement_in_place_matches_allocating_version(seq in iupac_seq(256)) {
        let mut buf = seq.clone();
        if let Ok(expected) = reverse_complement(&seq) {
            reverse_complement_in_place(&mut buf).unwrap();
            prop_assert_eq!(buf, expected);
        } else {
            // The in-place version must leave the slice untouched on error.
            prop_assert!(reverse_complement_in_place(&mut buf).is_err());
            prop_assert_eq!(buf, seq);
        }
    }

    #[test]
    fn base_complement_is_an_involution(i in 0u8..4) {
        let base = Base::from_2bit(i).unwrap();
        prop_assert_eq!(base.complement().complement(), base);
        prop_assert_ne!(base.complement(), base);
    }

    #[test]
    fn iupac_ascii_round_trip(b in prop::sample::select(IUPAC_BYTES)) {
        if let Some(code) = Iupac::from_ascii(b) {
            prop_assert_eq!(code.to_ascii(), b.to_ascii_uppercase());
            prop_assert_eq!(Iupac::from_nibble(code.to_nibble()), Some(code));
        }
    }

    #[test]
    fn kmer_reverse_complement_is_an_involution(
        k in 1usize..=kmer::MAX_K,
        raw in any::<u64>(),
    ) {
        let mask = if k == kmer::MAX_K { u64::MAX } else { (1u64 << (2 * k)) - 1 };
        let x = raw & mask;
        prop_assert_eq!(kmer::reverse_complement(kmer::reverse_complement(x, k), k), x);
        prop_assert_eq!(kmer::canonical(x, k), kmer::canonical(kmer::reverse_complement(x, k), k));
    }

    #[test]
    fn kmer_count_is_n_minus_k_plus_one(
        seq in acgt_seq(0, 300),
        k in 1usize..=kmer::MAX_K,
    ) {
        let count = Kmers::new(&seq, k).unwrap().count();
        let expected = (seq.len() + 1).saturating_sub(k);
        prop_assert_eq!(count, expected);
        let canon_count = Kmers::new(&seq, k).unwrap().canonical().count();
        prop_assert_eq!(canon_count, expected);
    }

    #[test]
    fn kmer_count_skips_windows_containing_n(
        seq in iupac_seq(200),
        k in 1usize..=8,
    ) {
        // Expected count: the sum over maximal ACGT runs of max(0, run - k + 1).
        let mut expected = 0usize;
        let mut run = 0usize;
        for &b in &seq {
            if Base::from_ascii(b).is_some() {
                run += 1;
            } else {
                expected += (run + 1).saturating_sub(k);
                run = 0;
            }
        }
        expected += (run + 1).saturating_sub(k);
        let got = Kmers::new(&seq, k).unwrap().count();
        prop_assert_eq!(got, expected);
    }

    #[test]
    fn first_kmer_equals_encode(seq in acgt_seq(1, 33)) {
        let k = seq.len();
        let first = Kmers::new(&seq, k).unwrap().next().unwrap();
        prop_assert_eq!(first, (0, kmer::encode(&seq).unwrap()));
    }

    #[test]
    fn interval_from_one_based_round_trips(start in 1u64..1_000_000, len in 0u64..1_000_000) {
        let end = start + len;
        let iv = Interval::from_one_based_inclusive(start, end).unwrap();
        prop_assert_eq!(iv.start() + 1, start);
        prop_assert_eq!(iv.end(), end);
        prop_assert_eq!(iv.len(), len + 1);
    }

    #[test]
    fn region_text_round_trips(
        contig in "[a-z][a-z0-9_]{0,7}",
        start in 1u64..1_000_000,
        len in 0u64..1_000_000,
    ) {
        let end = start + len;
        let text = format!("{contig}:{start}-{end}");
        let region = Region::parse(&text).unwrap();
        prop_assert_eq!(region.contig(), contig.as_str());
        prop_assert_eq!(region.interval(), Interval::from_one_based_inclusive(start, end));
    }

    #[test]
    fn interval_intersection_is_symmetric_and_bounded(
        a0 in 0u64..1000, al in 0u64..500,
        b0 in 0u64..1000, bl in 0u64..500,
    ) {
        let a = Interval::new(a0, a0 + al).unwrap();
        let b = Interval::new(b0, b0 + bl).unwrap();
        let ab = a.intersect(&b);
        prop_assert_eq!(ab, b.intersect(&a));
        prop_assert_eq!(a.overlaps(&b), ab.is_some());
        if let Some(i) = ab {
            prop_assert!(i.len() <= a.len().min(b.len()));
        }
    }
}
