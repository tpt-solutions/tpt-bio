//! Self-comparison benchmarks for the sequence parsers.
//!
//! These track this crate's own throughput over time. They make no comparison with other
//! tools: samtools and bcftools are not part of this benchmark run.

use std::fmt::Write as _;

use tpt_benchmark::tpt_bench_group;
use tpt_bio_seq::BamWriter;
use tpt_bio_seq::sam::SamReader;

const READS: u64 = 20_000;
const READ_LEN: usize = 100;

/// Deterministic pseudo-random bases, so the benchmark needs no fixtures.
fn bases(seed: u64, len: usize) -> String {
    let mut x = seed.wrapping_mul(0x9E37_79B9_7F4A_7C15) | 1;
    (0..len)
        .map(|_| {
            x ^= x << 13;
            x ^= x >> 7;
            x ^= x << 17;
            b"ACGT"[(x % 4) as usize] as char
        })
        .collect()
}

fn fastq_text(n: u64) -> Vec<u8> {
    let mut out = String::new();
    for i in 0..n {
        let seq = bases(i, READ_LEN);
        let _ = write!(out, "@read{i} lane=1\n{seq}\n+\n{}\n", "I".repeat(READ_LEN));
    }
    out.into_bytes()
}

fn sam_text(n: u64) -> String {
    let mut s = String::from("@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:100000000\n");
    for i in 0..n {
        let seq = bases(i + 7, READ_LEN);
        let pos = 1 + i * 50;
        let _ = writeln!(
            s,
            "read{i}\t0\tchr1\t{pos}\t60\t{READ_LEN}M\t*\t0\t0\t{seq}\t{}\tNM:i:0\tRG:Z:bench\n",
            "I".repeat(READ_LEN)
        );
    }
    s
}

fn vcf_lines(n: u64) -> Vec<String> {
    (0..n)
        .map(|i| {
            format!(
                "chr1\t{}\trs{i}\tA\tG\t50\tPASS\tDP={};AF=0.5\tGT:DP\t0/1:12\t1/1:9",
                100 + i * 10,
                i % 40
            )
        })
        .collect()
}

fn bam_bytes(sam: &str) -> Vec<u8> {
    let mut reader = SamReader::new(sam.as_bytes()).expect("valid SAM header");
    let header = reader.header().clone();
    let mut out = Vec::new();
    let mut w = BamWriter::new(&mut out, header).expect("BAM header");
    while let Some(rec) = reader.next_record().expect("valid SAM record") {
        w.write_record(&rec).expect("BAM record");
    }
    w.finish().expect("BAM finish");
    out
}

#[tpt_bench_group]
mod bgzf_decompress {
    use std::io::Write;

    use tpt_benchmark::{Bencher, Throughput, black_box, tpt_bench};
    use tpt_bio_seq::bgzf::{BgzfReader, BgzfWriter};

    use super::{READS, fastq_text};

    fn compressed_payload() -> (u64, Vec<u8>) {
        let payload = fastq_text(READS);
        let mut w = BgzfWriter::new(Vec::new());
        w.write_all(&payload).expect("in-memory write");
        let compressed = w.finish().expect("in-memory finish");
        (payload.len() as u64, compressed)
    }

    #[tpt_bench]
    fn parallel(b: &mut Bencher) {
        let (len, compressed) = compressed_payload();
        b.throughput(Throughput::Bytes(len));
        b.iter(|| black_box(tpt_bio_seq::bgzf::decompress_parallel(black_box(&compressed))));
    }

    #[tpt_bench]
    fn streaming(b: &mut Bencher) {
        let (len, compressed) = compressed_payload();
        b.throughput(Throughput::Bytes(len));
        b.iter(|| black_box(BgzfReader::new(black_box(&compressed[..])).read_to_vec()));
    }
}

#[tpt_bench_group]
mod fastq_parse {
    use tpt_benchmark::{Bencher, Throughput, black_box, tpt_bench};
    use tpt_bio_seq::FastqParser;

    use super::{READS, fastq_text};

    #[tpt_bench]
    fn records(b: &mut Bencher) {
        let data = fastq_text(READS);
        b.throughput(Throughput::Elements(READS));
        b.iter(|| {
            let n = FastqParser::new(black_box(&data)).flatten().count();
            black_box(n)
        });
    }
}

#[tpt_bench_group]
mod sam_parse {
    use tpt_benchmark::{Bencher, Throughput, black_box, tpt_bench};
    use tpt_bio_seq::sam::SamReader;

    use super::{READS, sam_text};

    #[tpt_bench]
    fn records(b: &mut Bencher) {
        let text = sam_text(READS);
        b.throughput(Throughput::Elements(READS));
        b.iter(|| {
            let mut reader = SamReader::new(black_box(text.as_bytes())).expect("SAM header");
            let mut n = 0u64;
            while let Some(rec) = reader.next_record().expect("SAM record") {
                black_box(&rec);
                n += 1;
            }
            n
        });
    }
}

#[tpt_bench_group]
mod vcf_record_parse {
    use tpt_benchmark::{Bencher, Throughput, black_box, tpt_bench};
    use tpt_bio_seq::VcfRecord;

    use super::{READS, vcf_lines};

    #[tpt_bench]
    fn records_2_samples(b: &mut Bencher) {
        let lines = vcf_lines(READS);
        b.throughput(Throughput::Elements(READS));
        b.iter(|| {
            let mut n = 0u64;
            for (i, line) in lines.iter().enumerate() {
                let rec = VcfRecord::parse(black_box(line), 2, i as u64 + 1).expect("VCF record");
                black_box(rec.pos());
                n += 1;
            }
            n
        });
    }
}

#[tpt_bench_group]
mod bam_record_view {
    use tpt_benchmark::{Bencher, Throughput, black_box, tpt_bench};
    use tpt_bio_seq::BamReader;

    use super::{READS, bam_bytes, sam_text};

    #[tpt_bench]
    fn views(b: &mut Bencher) {
        let bam = bam_bytes(&sam_text(READS));
        b.throughput(Throughput::Elements(READS));
        b.iter(|| {
            let mut reader = BamReader::new(black_box(&bam[..])).expect("BAM header");
            let mut n = 0u64;
            while let Some(view) = reader.next_view().expect("BAM record") {
                black_box(view.pos());
                n += 1;
            }
            n
        });
    }
}

tpt_benchmark::main!(bgzf_decompress, fastq_parse, sam_parse, vcf_record_parse, bam_record_view);
