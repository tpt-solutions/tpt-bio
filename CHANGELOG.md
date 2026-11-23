# Changelog

All notable changes are recorded here. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and the project uses
[Semantic Versioning](https://semver.org/) once 0.1.0 is released.

## [Unreleased]

### Added
- Workspace scaffold: `tpt-bio-core`, `tpt-bio-seq`, `tpt-bio-vision`,
  `tpt-bio-pipeline`, `tpt-bio-cli` (binary `tptbio`), `bindings/tpt-bio-py`.
- `tpt-bio-core`: nucleotide types (DNA and IUPAC, BAM 4-bit codes), amino acids and
  the standard codon table, Phred qualities, 0-based half-open intervals and region
  parsing, k-mer encoding and iteration, GC content and entropy, Poisson and Binomial
  distributions via `tpt-math`. `no_std` + `alloc`, with optional rayon helpers.
- `tpt-bio-seq`: pure Rust BGZF reader and writer with virtual offsets and parallel
  block (de)compression; FASTQ parsing, streaming reader, paired and interleaved
  readers, writer; SAM and BAM readers and writers with a zero-copy BAM record view;
  BAI build, serialization, and region fetch; memory-mapped file access.
- `tptbio seq view` (SAM/BAM, region queries, header-only, BAM output) and
  `tptbio seq index`.
- Project files: licenses, README, contributing guide, `deny.toml`, CI workflow,
  Dependabot configuration, fuzz scaffolding.
- `tpt-bio-seq`: FASTA reader and writer with the samtools `.fai` index; in-memory
  coordinate and name sorting; sorted-stream merge; flagstat counts; streaming batch
  processing of FASTQ with a bounded channel.
- `tpt-bio-seq`: record tables and CSV and JSON Lines export. `tpt-strata` is the
  in-process columnar table, the Arrow replacement. `tptbio seq flagstat`,
  `tptbio seq sort` and `tptbio seq export --format csv|jsonl` expose these.
- `tpt-bio-pipeline`: DAG model with cycle rejection, deterministic topological order,
  sequential and level-parallel execution, error propagation and cancellation.
- `tpt-bio-vision`: separable 2D convolution and Gaussian blur on `tpt-simd`, and an
  AWS SigV4 request signer on `tpt-crypto-hash` (checked against the AWS ListUsers
  example).
- `bindings/tpt-bio-py`: PyO3 module with stubs and pytest tests. It is in the main
  workspace, and `target-lexicon` comes from the MIT-licensed `target-lexicon-mit`
  patch.

### Changed
- Minimum supported Rust version raised to 1.95 for `tpt-simd` and `target-lexicon-mit`.
- The Poisson and Binomial functions require the `std` feature, because they use
  `tpt-math-prob-dist`, which needs `std`.

### Not yet implemented
- `tptbio vision`, `tptbio pipeline`, `tptbio stats`, `tptbio view` (they exit with
  status 2).
- VCF, BCF, CRAM, and CSI/TBI indexes.
- Image I/O, segmentation, feature extraction, the streaming pipeline executor, and
  Parquet export (out of scope; the project uses CSV and JSON Lines instead).
