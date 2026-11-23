# tpt-bio-seq

Pure Rust sequencing I/O for tpt-bio: FASTQ, BGZF, SAM, BAM, BAI, VCF, BCF,
tabix, FASTA and CRAM.

No C code is linked. Compression uses `flate2` with the `miniz_oxide` backend.
The only `unsafe` code is the memory-mapping wrapper in the `mmap` module, which
is enabled by the default `mmap` feature.

## Supported formats

| Format | Module | Read | Write | Notes |
| --- | --- | --- | --- | --- |
| BGZF | `bgzf` | yes | yes | Blocked gzip used by BAM, tabix and VCF.gz |
| FASTQ | `fastq` | yes | yes | Streaming record parser and `FastqWriter` |
| FASTA | `fasta` | yes | yes | Includes `.fai` index support and `FastaWriter` |
| SAM | `sam` | yes | yes | Text alignments |
| BAM | `bam` | yes | yes | Binary alignments over BGZF |
| BAI | `bai` | yes | yes | Region queries on sorted BAM; index files can be written |
| VCF | `vcf` | yes | yes | Header, records and metadata lines |
| BCF | `bcf` | yes | yes | Binary VCF |
| Tabix | `tabix` | yes | yes | `.tbi` index for BGZF text files; index files can be written |
| CRAM | `cram` | yes | yes | Containers, slices, rANS and Huffman codecs, reference handling |
| Tables | `table` | — | — | Export of alignment records through `tpt-strata` |

Supporting modules: `header` (shared SAM/BAM header model), `record` (common
record types), `ops` (coordinate and name sorting, merging sorted streams,
flagstat counts), `error`.

## Installation

```toml
[dependencies]
tpt-bio-seq = "0.1"
```

Disable memory mapping if you do not need it:

```toml
[dependencies]
tpt-bio-seq = { version = "0.1", default-features = false }
```

## Features

| Feature | Default | Effect |
| --- | --- | --- |
| `mmap` | yes | Memory-mapped file access through `memmap2`. Only module that uses `unsafe`. |

## Benchmarks and tests

```sh
cargo test -p tpt-bio-seq
cargo bench -p tpt-bio-seq --bench parsers
```

- `tests/malformed.rs` checks that truncated and corrupt inputs return errors
  rather than panicking.
- `tests/random_inputs.rs` feeds generated inputs to the parsers.
- `benches/parsers.rs` measures parser throughput with `criterion`.
- The fuzz targets in the top-level `fuzz/` directory cover BGZF, FASTQ, SAM,
  BAM, BAI, VCF, BCF, tabix and FASTA index parsing.

## Related crates

- `tpt-bio-core` supplies nucleotide, coordinate and quality types.
- `tpt-bio-vision` and `tpt-bio-pipeline` consume these readers.
- The `tptbio seq` subcommands in `tpt-bio-cli` expose view, index, flagstat,
  sort and export.

## Minimum supported Rust version

Rust 1.95, set in the workspace manifest.

## Licence

Dual-licensed under MIT OR Apache-2.0, at your option.
