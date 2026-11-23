# Changelog

All notable changes to `tpt-bio-seq` are recorded here. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and the crate uses
[Semantic Versioning](https://semver.org/).

## [Unreleased]

### Added

- VCF and BCF readers and writers.
- Tabix index reading and writing.
- CRAM codec: containers, slices, rANS and Huffman decoding, reference handling,
  and a writer with round-trip tests.
- FASTA `.fai` index support and `FastaWriter`.
- Malformed-input and random-input test suites.
- `criterion` benchmarks for the parsers.
- Explicit schema for table export, so all-NULL columns (such as RNEXT) keep a type.

### Changed

- BAM, BGZF, FASTA, header and memory-mapping paths updated to match the new readers.

## [0.1.0]

Initial workspace version. Not yet published.
