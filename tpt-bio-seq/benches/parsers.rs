//! Self-comparison benchmarks for the sequence parsers.
//!
//! These track this crate's own throughput over time. They make no comparison with other
//! tools: samtools and bcftools are not part of this benchmark run.

use std::fmt::Write as _;
use std::hint::black_box;
use std::io::Write;

use criterion::{BenchmarkId, Criterion, Throughput, criterion_group, criterion_main};
use tpt_bio_seq::bgzf::{BgzfReader, BgzfWriter};
use tpt_bio_seq::sam::SamReader;
use tpt_bio_seq::{BamReader, BamWriter, FastqParser, VcfRecord};

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

fn bench_bgzf_decompress(c: &mut Criterion) {
    let payload = fastq_text(READS);
    let mut w = BgzfWriter::new(Vec::new());
    w.write_all(&payload).expect("in-memory write");
    let compressed = w.finish().expect("in-memory finish");

    let mut group = c.benchmark_group("bgzf_decompress");
    group.throughput(Throughput::Bytes(payload.len() as u64));
    group.bench_function("parallel", |b| {
        b.iter(|| black_box(tpt_bio_seq::bgzf::decompress_parallel(black_box(&compressed))));
    });
    group.bench_function("streaming", |b| {
        b.iter(|| black_box(BgzfReader::new(black_box(&compressed[..])).read_to_vec()));
    });
    group.finish();
}

fn bench_fastq_parse(c: &mut Criterion) {
    let data = fastq_text(READS);
    let mut group = c.benchmark_group("fastq_parse");
    group.throughput(Throughput::Elements(READS));
    group.bench_function("records", |b| {
        b.iter(|| {
            let n = FastqParser::new(black_box(&data)).flatten().count();
            black_box(n)
        });
    });
    group.finish();
}

fn bench_sam_parse(c: &mut Criterion) {
    let text = sam_text(READS);
    let mut group = c.benchmark_group("sam_parse");
    group.throughput(Throughput::Elements(READS));
    group.bench_function("records", |b| {
        b.iter(|| {
            let mut reader = SamReader::new(black_box(text.as_bytes())).expect("SAM header");
            let mut n = 0u64;
            while let Some(rec) = reader.next_record().expect("SAM record") {
                black_box(&rec);
                n += 1;
            }
            n
        });
    });
    group.finish();
}

fn bench_vcf_record_parse(c: &mut Criterion) {
    let lines = vcf_lines(READS);
    let mut group = c.benchmark_group("vcf_record_parse");
    group.throughput(Throughput::Elements(READS));
    group.bench_function(BenchmarkId::new("records", "2 samples"), |b| {
        b.iter(|| {
            let mut n = 0u64;
            for (i, line) in lines.iter().enumerate() {
                let rec = VcfRecord::parse(black_box(line), 2, i as u64 + 1).expect("VCF record");
                black_box(rec.pos());
                n += 1;
            }
            n
        });
    });
    group.finish();
}

fn bench_bam_record_view(c: &mut Criterion) {
    let bam = bam_bytes(&sam_text(READS));
    let mut group = c.benchmark_group("bam_record_view");
    group.throughput(Throughput::Elements(READS));
    group.bench_function("views", |b| {
        b.iter(|| {
            let mut reader = BamReader::new(black_box(&bam[..])).expect("BAM header");
            let mut n = 0u64;
            while let Some(view) = reader.next_view().expect("BAM record") {
                black_box(view.pos());
                n += 1;
            }
            n
        });
    });
    group.finish();
}

criterion_group!(
    benches,
    bench_bgzf_decompress,
    bench_fastq_parse,
    bench_sam_parse,
    bench_vcf_record_parse,
    bench_bam_record_view
);
criterion_main!(benches);
