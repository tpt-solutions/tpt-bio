# Benchmarks

## Baseline (self-measured, Windows x86_64)

These numbers are a self-measured baseline on one Windows 11 x86_64 machine (16 GB RAM),
rustc 1.95.0, release profile. They measure this crate's own throughput over time.
They are **not** compared with samtools, bcftools or any other tool. Other machines,
other OSes and other background load will give different numbers, so only compare runs
made on the same machine.

Command:

```
cargo bench -p tpt-bio-seq --bench parsers -- --warm-up-time 1 --measurement-time 3
```

Criterion reports the low, median and high estimate of each measurement (95% confidence
interval). The table gives the median.

| Benchmark | Median time | Median throughput |
| --- | --- | --- |
| `bgzf_decompress/parallel` | 2.145 ms | 1.92 GiB/s |
| `bgzf_decompress/streaming` | 9.217 ms | 458 MiB/s |
| `fastq_parse/records` | 1.384 ms | 14.4 Melem/s |
| `sam_parse/records` | 14.54 ms | 1.38 Melem/s |
| `vcf_record_parse/records` (2 samples) | 6.953 ms | 2.88 Melem/s |
| `bam_record_view/views` | 7.114 ms | 2.81 Melem/s |

The short measurement time (3 s) keeps the run quick. Expect wider intervals than a
default `cargo bench` run. The `bgzf_decompress/parallel` and `streaming` rows use the
same input, so their ratio shows the parallel speed-up on this machine.

## Rerunning

- Run `cargo bench -p tpt-bio-seq --bench parsers` for the default criterion settings
  (longer, more precise).
- Criterion writes raw results to `target/criterion/`. Compare runs by reading the
  `change:` lines that criterion prints against the previous run in that directory.
- Close other heavy programs first. Do not run benchmarks while `cargo fuzz` or a build
  is running, since they share the CPU.
- The benchmark source is `tpt-bio-seq/benches/parsers.rs`. Its inputs are generated in
  code, so no data files are needed.

## Scheduled run

`.github/workflows/perf.yml` runs the same benchmark weekly on a hosted Linux runner and
keeps the output as a workflow artifact. It is informational: it does not fail on a
regression, and hosted runners are noisy, so read the trend across runs, not one value.

## Known issues found while benchmarking and fuzzing

- `bgzf::inflate_block` (`tpt-bio-seq/src/bgzf.rs`) reserves `ISIZE` bytes before checking
  it, and ISIZE is a 32-bit footer field, so one 100-byte block can request about 4 GB.
  Fuzzing found this through the `bgzf_decompress`, `bam_read` and `bcf_read` targets.
  The fix is to reject ISIZE above `MAX_BLOCK_SIZE` before allocating. The regression
  inputs are in `tpt-bio-seq/tests/data/fuzz_regressions/`.
