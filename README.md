# tpt-bio

A pure Rust engine for bioinformatics sequence data and biological image analysis.
It is intended as a zero-FFI replacement for HTSlib (sequence I/O) and CellProfiler
(image analysis).

> **Status: early development (0.1.0-dev).** Sequence I/O is usable. Image analysis,
> VCF/BCF, CRAM, the pipeline engine and Python bindings are not implemented yet.
> See [todo.md](todo.md) for the phased roadmap.

## What works today

- **BGZF**: pure Rust blocked-gzip reader and writer, with virtual offsets, seeking,
  EOF-marker detection and parallel block (de)compression.
- **FASTQ**: zero-copy parser over byte slices, streaming reader, paired and
  interleaved input, writer. Gzip and BGZF input are detected by magic bytes.
- **SAM and BAM**: text and binary readers and writers with a shared record model,
  a zero-copy BAM record view, and round-trip fidelity.
- **BAI**: build from a coordinate-sorted BAM, read and write, and region queries.
- **Core types**: nucleotides (DNA and IUPAC), amino acids and codon tables, Phred
  qualities, 0-based half-open intervals and region parsing, k-mers, GC content and
  entropy, and Poisson and Binomial distributions (via `tpt-math`).

```sh
tptbio seq view reads.sam -b -o reads.bam   # SAM -> BAM
tptbio seq index reads.bam                  # writes reads.bam.bai
tptbio seq view reads.bam -r chr1:1000-2000 # region query (1-based, inclusive)
```

## Build

The workspace pins Rust 1.85 through `rust-toolchain.toml`.

```sh
cargo build --workspace
cargo test --workspace
cargo build -p tpt-bio-core --no-default-features   # no_std + alloc build
```

## Design

- **Permissive dependencies only.** Every dependency must be MIT, BSD, ISC, Zlib,
  Unicode-3.0, or a dual `MIT OR Apache-2.0` crate. Apache-2.0-only and copyleft
  crates are rejected. `cargo deny check` enforces this in CI.
- **No C code.** No `-sys` crates and no bindgen. Compression uses `flate2` with
  the `miniz_oxide` backend.
- **Minimal unsafe.** The only `unsafe` code is the memory-mapping wrapper in
  `tpt-bio-seq`, documented in [docs/unsafe-policy.md](docs/unsafe-policy.md).
- **`no_std` core.** `tpt-bio-core` builds without `std` for embedded and
  edge-device use.

Full design notes are in [spec.txt](spec.txt).

## Crates

| Crate | Purpose |
| --- | --- |
| `tpt-bio-core` | Shared types, statistics, bio-math (`no_std` + `alloc`) |
| `tpt-bio-seq` | FASTQ, BGZF, SAM, BAM, BAI |
| `tpt-bio-vision` | Image I/O and analysis (not started) |
| `tpt-bio-pipeline` | Workflow engine (not started) |
| `tpt-bio-cli` | The `tptbio` command-line tool |
| `bindings/tpt-bio-py` | Python bindings (not started) |

## License

Dual-licensed under either of

- [MIT license](LICENSE-MIT)
- [Apache License, Version 2.0](LICENSE-APACHE)

at your option. Copyright © 2026 TPT Solutions.

Unless you explicitly state otherwise, any contribution intentionally submitted for
inclusion in this work, as defined in the Apache-2.0 license, shall be dual licensed
as above, without any additional terms or conditions. See
[CONTRIBUTING.md](CONTRIBUTING.md).
