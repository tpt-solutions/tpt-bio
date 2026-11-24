# tpt-bio

A pure Rust engine for bioinformatics sequence data and biological image analysis.
It is intended as a zero-FFI replacement for HTSlib (sequence I/O) and CellProfiler
(image analysis) in the library crates. The Python bindings use PyO3.

> **Status: early development (0.1.0-dev).** Sequence I/O, image analysis, pipelines,
> the `tptbio` command line and Python bindings all exist, but several parts are
> partial. Nothing is published to crates.io or PyPI yet. See [todo.md](todo.md) for the
> roadmap and the status of each item.

## What works

- **Sequence** (`tpt-bio-seq`): BGZF; FASTQ, streaming, paired and interleaved, gzip and
  BGZF input; SAM and BAM with a zero-copy BAM view; BAI index build and region query;
  VCF (plain and BGZF), BCF, tabix and CSI indexes with region queries; FASTA with `.fai`;
  coordinate and name sort, merge, flagstat; CSV and JSON Lines export.
- **Vision** (`tpt-bio-vision`): PNG, TIFF and BigTIFF, OME-TIFF, Zarr v2 and v3 reading,
  and OME-NGFF plate metadata; global and tile-based Otsu thresholding; 2D and 3D
  labelling; marker-controlled watershed; LoG seeds; morphology; Gaussian, median and
  mean filters; background subtraction; Sobel, Scharr and Canny edges; per-object shape,
  intensity, texture and colocalization features; IoU and F1 matching.
- **3D rendering** (`tpt-bio-vision::render3d`): a CPU reference renderer with ray casting,
  label overlay, marching-tetrahedra isosurfaces and a flat-shaded rasteriser. See
  [docs/viewer.md](docs/viewer.md).
- **Core** (`tpt-bio-core`): nucleotides and IUPAC codes, amino acids and codon tables,
  Phred scores, intervals and region parsing, k-mers, GC content and entropy, Poisson and
  Binomial distributions, normalization, PCA, k-means, quality metrics.
- **Pipelines** (`tpt-bio-pipeline`): a DAG with typed ports and cycle checks, TOML
  definitions, built-in sequence and 2D vision nodes, parallel and checkpointed execution.
- **Command line** (`tptbio`): `seq` (view, index, flagstat, sort, export), `vision`
  (segment, features), `pipeline run`, `stats` (normalize, pca), `view` (image summary),
  shell completions and man pages.

```sh
tptbio seq sort reads.sam -o reads.bam        # SAM -> coordinate-sorted BAM
tptbio seq index reads.bam                    # writes reads.bam.bai
tptbio seq view reads.bam -r chr1:1000-2000   # region query (1-based, inclusive)
```

## Partial

- **CRAM**: CRAM 3.0 reader and writer. rANS 4x8, Huffman and gzip blocks work; rANS
  order 1 and Nx16, bzip2 and LZMA are not implemented. Not yet compared with samtools.
- **VCF, BCF, tabix and CSI**: tested on generated data. No comparison with bcftools.
- **Zarr**: reads v2 and v3 with no compression or gzip; blosc and zstd are not supported.
  Writes v2 only.
- **Remote reads**: HTTP range reads and an AWS SigV4 signer exist. Only plain `http://`
  is supported, because the project has no TLS.
- **HCS plates**: plate and well metadata work. Field and channel layout is not verified.
- **3D viewer**: CPU reference only. No chunked or streamed volumes.
- **Features**: a few dozen per-object features. The 500+ target is not met.
- **Statistics**: t-SNE and UMAP are simplified and not checked against references.
- **CLI**: `vision segment` uses the global Otsu threshold only. `tptbio view` prints a
  summary and does not display images.

## Not implemented

- A GPU (wgpu) viewer, a WebAssembly or WebGPU build, and a desktop app (spec ADR-003).
- Conformance tests against samtools, bcftools, scikit-learn or CellProfiler, and
  benchmarks against samtools.
- Parquet. It needs Arrow, which is excluded by the licence policy.
- HTTPS, and cloud backends other than S3 (GCS, Azure).
- 3D steps, and OME-TIFF series and plane selection, in the TOML pipeline format.
  Its `load_image` node reads PNG and TIFF pages only.
- Prebuilt binaries and published packages.

## Tutorials

- [Sequence quickstart](docs/tutorials/seq-quickstart.md): SAM, sort, index, region query, export.
- [Vision quickstart](docs/tutorials/vision-quickstart.md): PNG, threshold, label, filter, measure.
- [Pipelines](docs/tutorials/pipeline.md): TOML pipelines from the command line and from Rust.
- [Python](docs/tutorials/python.md): build a wheel, load an image, measure objects, run a pipeline.
- [3D rendering](docs/viewer.md): render a synthetic multi-channel volume on the CPU.

Other documents: [API review](docs/api-review.md),
[3D viewer architecture](docs/architecture-3d-viewer.md),
[unsafe policy](docs/unsafe-policy.md), [dependency audit](docs/dependency-audit.md),
[test data](docs/test-data.md), [spec](spec.txt).

## Build

The workspace uses Rust 1.95, the minimum supported version (`rust-version` in
`Cargo.toml`). `rust-toolchain.toml` pins the same toolchain.

```sh
cargo build                 # the default members: core, seq, vision, pipeline and the CLI
cargo test --workspace
cargo build -p tpt-bio-core --no-default-features   # no_std + alloc build
```

The first build needs network access. Two dependencies come from git: `tpt-crypto-hash`,
pinned to a commit, and a MIT-only `target-lexicon` replacement that is patched in the
root `Cargo.toml`. Python bindings are built separately with maturin; see the Python
tutorial.

## Design

- **Permissive dependencies only.** Every dependency must be MIT, BSD, ISC, Zlib,
  Unicode-3.0, or dual MIT/Apache-2.0. Copyleft and Apache-2.0-only crates are rejected.
  `cargo deny check` enforces this in CI.
- **No C code in the library crates.** There are no `-sys` crates and no bindgen.
  Compression uses `flate2` with the `miniz_oxide` backend.
- **Minimal unsafe.** The only `unsafe` code is the memory-mapping wrapper
  `MappedFile::open` in `tpt-bio-seq`, documented in [docs/unsafe-policy.md](docs/unsafe-policy.md).
- **`no_std` core.** `tpt-bio-core` builds without `std` for embedded and edge-device use.

Full design notes are in [spec.txt](spec.txt).

## Crates

| Crate | Purpose |
| --- | --- |
| `tpt-bio-core` | Shared types, statistics and bio-math (`no_std` + `alloc`) |
| `tpt-bio-seq` | FASTQ, BGZF, SAM, BAM, BAI, VCF, BCF, tabix, CRAM, FASTA |
| `tpt-bio-vision` | Image I/O, segmentation, features and the CPU 3D renderer |
| `tpt-bio-pipeline` | DAG model, TOML pipeline definitions and built-in nodes |
| `tpt-bio-cli` | The `tptbio` command-line tool |
| `bindings/tpt-bio-py` | Python bindings (PyO3), built locally with maturin |

## License

Dual-licensed under either of

- [MIT license](LICENSE-MIT)
- [Apache License, Version 2.0](LICENSE-APACHE)

at your option. Copyright © 2026 TPT Solutions.

Unless you explicitly state otherwise, any contribution intentionally submitted for
inclusion in this work, as defined in the Apache-2.0 license, shall be dual licensed
as above, without any additional terms or conditions. See
[CONTRIBUTING.md](CONTRIBUTING.md).
