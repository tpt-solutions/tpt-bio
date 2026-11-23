//! Round-trip tests: records written by [`CramWriter`] are read back by
//! [`CramReader`], and BAM records survive a trip through CRAM.

use std::io::Cursor;

use crate::bam::{BamReader, BamWriter};
use crate::cram::crai::CraiEntry;
use crate::cram::decode::{substitute, substitution_code};
use crate::cram::reader::CramReader;
use crate::cram::reference::CramReference;
use crate::cram::writer::CramWriter;
use crate::error::SeqError;
use crate::header::Header;
use crate::record::{CigarKind, CigarOp, Flags, Record, Tag, TagValue};

const HEADER: &str = "@HD\tVN:1.6\tSO:unsorted\n@SQ\tSN:chr1\tLN:1000\n@SQ\tSN:chr2\tLN:500\n";

/// Deterministic pseudo-random bases.
fn bases(seed: u64, len: usize) -> Vec<u8> {
    let mut s = seed;
    (0..len)
        .map(|_| {
            s = s.wrapping_mul(6_364_136_223_846_793_005).wrapping_add(1_442_695_040_888_963_407);
            b"ACGT"[(s >> 33) as usize % 4]
        })
        .collect()
}

fn fasta(chr1_seed: u64, chr2_seed: u64) -> Vec<u8> {
    let mut out = Vec::new();
    for (name, seq) in [("chr1", bases(chr1_seed, 1000)), ("chr2", bases(chr2_seed, 500))] {
        out.extend_from_slice(format!(">{name}\n").as_bytes());
        for line in seq.chunks(60) {
            out.extend_from_slice(line);
            out.push(b'\n');
        }
    }
    out
}

fn reference() -> CramReference {
    CramReference::from_fasta_bytes(fasta(1, 2)).unwrap()
}

fn header() -> Header {
    Header::parse(HEADER).unwrap()
}

fn quals(n: usize) -> Vec<u8> {
    (0..n).map(|i| (i * 7 % 40) as u8).collect()
}

#[allow(clippy::too_many_arguments)]
fn read(
    name: &str,
    flags: u16,
    ref_id: i32,
    pos: i64,
    mapq: u8,
    cigar: &[(u32, CigarKind)],
    seq: Vec<u8>,
    tags: Vec<Tag>,
) -> Record {
    Record {
        name: name.as_bytes().to_vec(),
        flags: Flags(flags),
        ref_id,
        pos,
        mapq,
        cigar: cigar.iter().map(|&(len, kind)| CigarOp { len, kind }).collect(),
        next_ref_id: -1,
        next_pos: -1,
        tlen: 0,
        qual: quals(seq.len()),
        seq,
        tags,
    }
}

fn slice_of(chr: &[u8], beg: usize, len: usize) -> Vec<u8> {
    chr[beg..beg + len].to_vec()
}

/// A mix covering mapped, reverse, clipped, indel, skip, unmapped and tagged reads.
fn sample_records(chr1: &[u8], chr2: &[u8]) -> Vec<Record> {
    use CigarKind::{Deletion, HardClip, Insertion, Match, Skip, SoftClip};
    let mut out = Vec::new();
    out.push(read("perfect", 0, 0, 99, 60, &[(50, Match)], slice_of(chr1, 99, 50), vec![]));
    let mut mismatch = slice_of(chr1, 200, 20);
    mismatch[5] = if mismatch[5] == b'A' { b'C' } else { b'A' };
    out.push(read("mismatch", 0, 0, 200, 30, &[(20, Match)], mismatch, vec![]));
    out.push(read(
        "reverse",
        Flags::REVERSE,
        0,
        300,
        45,
        &[(30, Match)],
        slice_of(chr1, 300, 30),
        vec![],
    ));
    let mut clipped = bases(9, 5);
    clipped.extend(slice_of(chr1, 400, 20));
    out.push(read("softclip", 0, 0, 400, 20, &[(5, SoftClip), (20, Match)], clipped, vec![]));
    let mut ins = slice_of(chr1, 500, 10);
    ins.extend_from_slice(b"ACG");
    ins.extend(slice_of(chr1, 510, 10));
    out.push(read(
        "insertion",
        0,
        0,
        500,
        50,
        &[(10, Match), (3, Insertion), (10, Match)],
        ins,
        vec![],
    ));
    let mut del = slice_of(chr1, 600, 10);
    del.extend(slice_of(chr1, 614, 10));
    out.push(read(
        "deletion",
        0,
        0,
        600,
        50,
        &[(10, Match), (4, Deletion), (10, Match)],
        del,
        vec![],
    ));
    let mut skip = slice_of(chr1, 700, 10);
    skip.extend(slice_of(chr1, 810, 10));
    out.push(read(
        "refskip",
        0,
        0,
        700,
        40,
        &[(10, Match), (100, Skip), (10, Match)],
        skip,
        vec![],
    ));
    out.push(read(
        "hardclip",
        0,
        0,
        800,
        40,
        &[(4, HardClip), (10, Match)],
        slice_of(chr1, 800, 10),
        vec![],
    ));
    // Read base N against reference base: a substitution to N.
    let mut with_n = slice_of(chr1, 900, 12);
    with_n[3] = b'N';
    out.push(read("n_base", 0, 0, 900, 33, &[(12, Match)], with_n, vec![]));
    // Mate information and a template.
    let mut mate = read(
        "paired",
        Flags::PAIRED | Flags::READ1 | Flags::PROPER_PAIR | Flags::MATE_REVERSE,
        0,
        10,
        60,
        &[(25, Match)],
        slice_of(chr1, 10, 25),
        vec![],
    );
    mate.next_ref_id = 0;
    mate.next_pos = 1000 - 1;
    mate.tlen = -350;
    out.push(mate);
    out.push(read("chr2read", 0, 1, 77, 12, &[(40, Match)], slice_of(chr2, 77, 40), vec![]));
    // Unmapped, and placed-but-unmapped.
    out.push(read("unmapped", Flags::UNMAPPED, -1, -1, 0, &[], b"ACGTNACGT".to_vec(), vec![]));
    out.push(read("placed", Flags::UNMAPPED, 1, 42, 0, &[], b"GGCC".to_vec(), vec![]));
    // Tags of several types.
    let tags = vec![
        Tag { key: *b"NM", value: TagValue::Int(1) },
        Tag { key: *b"RG", value: TagValue::String(b"grp1".to_vec()) },
        Tag { key: *b"XA", value: TagValue::Char(b'x') },
        Tag { key: *b"XF", value: TagValue::Float(1.5) },
        Tag { key: *b"XB", value: TagValue::IntArray(vec![1, -2, 3]) },
        Tag { key: *b"XH", value: TagValue::Hex(b"1AE3".to_vec()) },
    ];
    out.push(read("tagged", 0, 0, 950, 37, &[(8, Match)], slice_of(chr1, 950, 8), tags));
    out
}

fn write_all(
    records: &[Record],
    reference: Option<CramReference>,
    per_slice: usize,
    per_container: usize,
) -> Vec<u8> {
    let mut w =
        CramWriter::with_layout(Vec::new(), header(), reference, per_slice, per_container).unwrap();
    for r in records {
        w.write_record(r).unwrap();
    }
    w.finish().unwrap()
}

fn read_all<R: std::io::Read>(mut r: CramReader<R>) -> Result<Vec<Record>, SeqError> {
    let mut out = Vec::new();
    while let Some(rec) = r.next_record()? {
        out.push(rec);
    }
    Ok(out)
}

#[test]
fn substitution_codes_invert() {
    let m = [0x1B; 5];
    for &r in b"ACGTN" {
        for &alt in b"ACGTN" {
            if alt == r {
                assert_eq!(substitution_code(&m, r, alt), None);
            } else {
                let code = substitution_code(&m, r, alt).expect("alternative has a code");
                assert_eq!(substitute(&m, r, code).unwrap(), alt, "ref {r} alt {alt}");
            }
        }
    }
}

#[test]
fn reference_roundtrip_with_multiple_slices_and_containers() {
    let reference = reference();
    let records = sample_records(
        &reference.fetch("chr1", 0, 1000).unwrap(),
        &reference.fetch("chr2", 0, 500).unwrap(),
    );
    let bytes = write_all(&records, Some(reference.clone()), 3, 2);
    let got = read_all(CramReader::with_reference(Cursor::new(bytes), reference).unwrap()).unwrap();
    assert_eq!(got, records);
}

#[test]
fn reference_free_roundtrip() {
    let reference = reference();
    let records = sample_records(
        &reference.fetch("chr1", 0, 1000).unwrap(),
        &reference.fetch("chr2", 0, 500).unwrap(),
    );
    let bytes = write_all(&records, None, 4, 1);
    let got = read_all(CramReader::new(Cursor::new(bytes)).unwrap()).unwrap();
    assert_eq!(got, records);
}

#[test]
fn reference_needed_without_reference_is_unsupported() {
    let reference = reference();
    let records = sample_records(
        &reference.fetch("chr1", 0, 1000).unwrap(),
        &reference.fetch("chr2", 0, 500).unwrap(),
    );
    let bytes = write_all(&records, Some(reference), 100, 1);
    let mut r = CramReader::new(Cursor::new(bytes)).unwrap();
    let err = loop {
        match r.next_record() {
            Ok(Some(_)) => {}
            Ok(None) => panic!("expected an Unsupported error, got end of file"),
            Err(e) => break e,
        }
    };
    assert!(matches!(err, SeqError::Unsupported(_)), "got {err}");
}

#[test]
fn wrong_reference_fails_md5_check() {
    let reference = reference();
    let records = sample_records(
        &reference.fetch("chr1", 0, 1000).unwrap(),
        &reference.fetch("chr2", 0, 500).unwrap(),
    );
    let bytes = write_all(&records, Some(reference), 100, 1);
    let other = CramReference::from_fasta_bytes(fasta(3, 2)).unwrap();
    let err = read_all(CramReader::with_reference(Cursor::new(bytes), other).unwrap()).unwrap_err();
    assert!(matches!(err, SeqError::Malformed { .. }), "got {err}");
}

#[test]
fn bam_cram_bam_is_lossless() {
    let reference = reference();
    let records = sample_records(
        &reference.fetch("chr1", 0, 1000).unwrap(),
        &reference.fetch("chr2", 0, 500).unwrap(),
    );
    // BAM -> records -> CRAM -> records -> BAM -> records must reproduce the first set.
    let mut bw = BamWriter::new(Vec::new(), header()).unwrap();
    for r in &records {
        bw.write_record(r).unwrap();
    }
    let bam = bw.finish().unwrap();
    let mut br = BamReader::new(Cursor::new(bam)).unwrap();
    let mut from_bam = Vec::new();
    while let Some(r) = br.next_record().unwrap() {
        from_bam.push(r);
    }

    let cram = write_all(&from_bam, Some(reference.clone()), 5, 2);
    let from_cram =
        read_all(CramReader::with_reference(Cursor::new(cram), reference).unwrap()).unwrap();

    let mut bw = BamWriter::new(Vec::new(), header()).unwrap();
    for r in &from_cram {
        bw.write_record(r).unwrap();
    }
    let bam2 = bw.finish().unwrap();
    let mut br = BamReader::new(Cursor::new(bam2)).unwrap();
    let mut back = Vec::new();
    while let Some(r) = br.next_record().unwrap() {
        back.push(r);
    }
    assert_eq!(back, from_bam);
}

#[test]
fn empty_file_round_trips() {
    let bytes = write_all(&[], None, 10, 1);
    let mut r = CramReader::new(Cursor::new(bytes)).unwrap();
    assert_eq!(r.header().text(), HEADER);
    assert!(r.next_record().unwrap().is_none());
}

#[test]
fn truncated_file_is_an_error_not_a_panic() {
    let reference = reference();
    let records = sample_records(
        &reference.fetch("chr1", 0, 1000).unwrap(),
        &reference.fetch("chr2", 0, 500).unwrap(),
    );
    let bytes = write_all(&records, None, 4, 2);
    for cut in [0, 10, 30, bytes.len() / 2, bytes.len() - 1] {
        let part = bytes[..cut].to_vec();
        let res = CramReader::new(Cursor::new(part)).and_then(read_all);
        if cut <= bytes.len() / 2 {
            assert!(res.is_err(), "cut at {cut} should fail");
        }
    }
}

/// `n` reads of 30 bases on chr1 at scattered positions.
fn sorted_reads(chr1: &[u8], n: usize) -> Vec<Record> {
    (0..n)
        .map(|i| {
            let pos = (i * 13) % 900;
            read(
                &format!("r{i}"),
                0,
                0,
                pos as i64,
                60,
                &[(30, CigarKind::Match)],
                slice_of(chr1, pos, 30),
                vec![],
            )
        })
        .collect()
}

#[test]
fn crai_region_query_matches_brute_force() {
    let reference = reference();
    let chr1 = reference.fetch("chr1", 0, 1000).unwrap();
    let mut records = sorted_reads(&chr1, 60);
    records.sort_by_key(|r| r.pos);
    let per_slice = 5;
    let bytes = write_all(&records, Some(reference), per_slice, 1);
    let mut idx_reader = CramReader::new(Cursor::new(bytes)).unwrap();
    let crai = idx_reader.build_index().unwrap();
    let entries = crai.entries().to_vec();
    assert_eq!(entries.len(), records.len().div_ceil(per_slice));

    // Brute force from the records: slice k holds records[k*5..(k+1)*5].
    let slice_span: Vec<(i64, i64)> = records
        .chunks(per_slice)
        .map(|c| {
            let beg = c.iter().map(|r| r.pos).min().unwrap();
            let end = c.iter().map(|r| r.pos + 30).max().unwrap();
            (beg, end)
        })
        .collect();

    for (beg, end) in [(0i64, 10i64), (100, 150), (250, 251), (899, 1000), (500, 500), (0, 1000)] {
        let got: Vec<CraiEntry> = crai.query(0, beg, end);
        let expected: Vec<CraiEntry> = entries
            .iter()
            .zip(&slice_span)
            .filter(|(_, (b, e))| *b < end && beg < *e)
            .map(|(entry, _)| *entry)
            .collect();
        assert_eq!(got, expected, "region [{beg}, {end})");
        // No false negatives: every slice with a record overlapping the region is returned.
        for (k, chunk) in records.chunks(per_slice).enumerate() {
            let hit = chunk.iter().any(|r| r.pos < end && r.pos + 30 > beg);
            if hit {
                assert!(got.contains(&entries[k]), "slice {k} missing for [{beg}, {end})");
            }
        }
    }
}

#[test]
fn crai_round_trips_through_bytes() {
    let reference = reference();
    let chr1 = reference.fetch("chr1", 0, 1000).unwrap();
    let records = sorted_reads(&chr1, 12);
    let bytes = write_all(&records, Some(reference), 4, 1);
    let mut r = CramReader::new(Cursor::new(bytes)).unwrap();
    let crai = r.build_index().unwrap();
    let encoded = crai.to_bytes().unwrap();
    let decoded = crate::cram::crai::Crai::from_bytes(&encoded).unwrap();
    assert_eq!(decoded, crai);
}
