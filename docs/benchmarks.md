# Benchmarks

## Baseline (self-measured, Windows x86_64)

The benchmarks use [`tpt-benchmark`](https://github.com/tpt-solutions/tpt-benchmark),
pinned in the workspace `Cargo.toml`. The earlier `criterion` figures were removed: their
benchmark names and numbers are not comparable with `tpt-benchmark`, so the table below
is the first `tpt-benchmark` baseline. It is a single `--quick` run. Re-measure with a full
run before drawing conclusions. The ± value is the tool's relative error for that run.

| Benchmark | Time per iteration |
| --- | --- |
| `bgzf_decompress::parallel` | 2.83 ms ±4.6% |
| `bgzf_decompress::streaming` | 10.42 ms ±2.0% |
| `fastq_parse::records` | 1.46 ms ±1.6% |
| `sam_parse::records` | 16.41 ms ±13.9% (noisy) |
| `vcf_record_parse::records_2_samples` | 4.82 ms ±13.7% (noisy) |
| `bam_record_view::views` | 7.57 ms ±1.8% |

That run was on one Windows 11 x86_64 machine (16 GB RAM), release profile, `--quick`.
`--quick` is short and imprecise. Other machines, OSes and background load give different
numbers, so only compare runs made on the same machine.

## Rerunning

- Run `cargo bench -p tpt-bio-seq --bench parsers` for the default settings (longer,
  more precise). Use `-- --quick` for a shorter run.
- Results are written to `target/tpt-benchmark/`: `tpt-benchmark.json` for the run,
  `index.html` for the report. Save a named baseline with
  `-- --save-baseline <name>` and compare with `-- --baseline <name>`. A significant
  regression above the threshold (default 5%) exits with code 1.
- Run `cargo bench -p tpt-bio-seq --bench parsers -- --list` to list the benchmark names.
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
