# tpt-bio TODO

Dual-licensed MIT OR Apache-2.0 · © 2026 TPT Solutions · Pure Rust, zero-FFI, no Apache-2.0-only dependencies.

**Decisions locked in** (override `spec.txt` where they differ)
- cargo-deny allows only MIT, BSD-2/3, ISC, Zlib, Unicode-3.0 and `MIT OR Apache-2.0` expressions (no standalone Apache-2.0)
- Compression is flate2 (miniz_oxide) only. No zstd, since it wraps C
- tpt-math and tpt-simd come from crates.io
- Sequence parsers are fully custom (memmap2 + flate2); noodles is used only as a test oracle
- Remote I/O is HTTP range reads through a permissive client plus our own SigV4 signer. No aws-sdk or object_store
- Crypto comes from tpt-crypto, not `ring`
- CellProfiler parity is functional only (no .cppipe import)
- CI is GitHub Actions; Python wheels are built with maturin

---

## Phase 0: Project Foundations

### Repo & licensing
- [ ] `git init`, `.gitignore`, default branch `master`
- [ ] `LICENSE-MIT` (© 2026 TPT Solutions)
- [ ] `LICENSE-APACHE` (© 2026 TPT Solutions)
- [ ] License headers / `license = "MIT OR Apache-2.0"` in every Cargo.toml
- [ ] `README.md` (overview, license section, "MIT-first" dependency policy)
- [ ] `CONTRIBUTING.md` (dual-license contribution clause, dependency policy)
- [ ] `SECURITY.md`, `CODE_OF_CONDUCT.md`
- [ ] `CHANGELOG.md`

### Workspace scaffold
- [ ] Root `Cargo.toml` workspace (shared `[workspace.package]`, `[workspace.dependencies]`, lints)
- [ ] `rust-toolchain.toml` and pinned MSRV
- [ ] `rustfmt.toml`, `clippy.toml`
- [ ] Scaffold `tpt-bio-core` (`#![no_std]` + `alloc`)
- [ ] Scaffold `tpt-bio-seq`
- [ ] Scaffold `tpt-bio-vision`
- [ ] Scaffold `tpt-bio-pipeline`
- [ ] Scaffold `tpt-bio-cli` (binary `tptbio`)
- [ ] Scaffold `bindings/tpt-bio-py` (PyO3 + maturin `pyproject.toml`)
- [ ] Confirm tpt-math and tpt-simd versions on crates.io and their licenses/features (including `no_std` support)
- [ ] `unsafe` policy: `#![forbid(unsafe_code)]` where possible; document and isolate any exceptions (mmap, SIMD)

### Dependency license audit
- [ ] Verify license expression of each planned dependency: tpt-math, tpt-simd, tpt-crypto, num-traits, flate2 (+ miniz_oxide), memmap2, rayon, crossbeam, image, tiff, ndarray, clap, pyo3, wgpu, ureq, serde/toml
- [ ] Reject or replace anything Apache-2.0-only, MPL, or copyleft (including transitive)
- [ ] Use `tpt-crypto` instead of `ring` for all crypto needs (SigV4 HMAC/SHA-256, HTTPS/TLS client, checksums); drop the `ring` clarify block from `deny.toml`
- [ ] Verify tpt-crypto's license, `no_std` support and crates.io version; make sure the HTTP client's TLS stack can use it (or document the gap)
- [ ] Record the audit result in `docs/dependency-audit.md`

### cargo-deny & CI
- [ ] `deny.toml` with the corrected policy (MIT-option expressions only, copyleft denied, advisories, bans, sources)
- [ ] Verify `cargo deny check` fails on a deliberately Apache-only test crate, then remove it
- [ ] GitHub Actions: fmt, clippy `-D warnings`, test on Linux/macOS/Windows
- [ ] CI: `cargo deny check`
- [ ] CI: `no_std` build of `tpt-bio-core` (e.g. `thumbv7em-none-eabihf`)
- [ ] CI: MSRV job
- [ ] CI: docs build (`cargo doc`, deny warnings)
- [ ] CI: dependency caching and nightly scheduled run
- [ ] Dependabot / Renovate config

### Test data & tooling
- [ ] Choose permissively-licensed test data (small BAM/CRAM/VCF/FASTQ, OME-TIFF, Zarr)
- [ ] Set up `criterion` (check license) or a custom benchmark harness
- [ ] Set up `cargo-fuzz` infrastructure
- [ ] Install samtools/bcftools locally for oracle comparisons (test-only, never linked)

---

## Phase 1: Sequence Core (`tpt-bio-core` + `tpt-bio-seq` basics)

### tpt-bio-core (`no_std` + alloc)
- [ ] Nucleotide types (DNA/RNA, IUPAC codes, 2-bit and 4-bit encodings)
- [ ] Amino acid type and codon table
- [ ] Phred quality score type and conversions
- [ ] Genomic coordinate / region types (0-based vs 1-based clarity)
- [ ] Reverse complement, k-mer iteration
- [ ] Poisson distribution (PMF, CDF, quantile) via tpt-math
- [ ] Binomial distribution (PMF, CDF, quantile) via tpt-math
- [ ] Basic bio-math (GC content, entropy, etc.)
- [ ] Parallel iterator helpers (`std` feature, rayon-based)
- [ ] `no_std` and `std` feature matrix, tested in CI
- [ ] Unit and property tests; rustdoc with examples

### BGZF (flate2, pure Rust)
- [ ] BGZF block reader (header, extra field, CRC32, ISIZE validation)
- [ ] EOF marker handling
- [ ] Parallel block decompression (rayon/crossbeam)
- [ ] BGZF writer with parallel block compression
- [ ] Virtual file offsets (block offset + uncompressed offset)
- [ ] Seek by virtual offset
- [ ] Mmap-backed reader (memmap2)
- [ ] Corruption and truncation tests, fuzz target

### FASTQ
- [ ] Zero-copy FASTQ record parser
- [ ] Gzip / BGZF FASTQ input
- [ ] Paired-end / interleaved handling
- [ ] FASTQ writer
- [ ] Streaming parallel batch reader
- [ ] Tests and fuzz target

### SAM / BAM
- [ ] SAM header parser (@HD, @SQ, @RG, @PG, @CO)
- [ ] SAM record parser (text) and writer
- [ ] BAM header and reference dictionary
- [ ] Zero-copy BAM record view (cigar, seq, qual, aux tags)
- [ ] CIGAR operations and alignment-end calculation
- [ ] Auxiliary tag parsing (all types, arrays)
- [ ] BAM writer
- [ ] SAM <-> BAM round-trip tests
- [ ] Conformance tests vs samtools and noodles outputs

### BAI index & random access
- [ ] BAI parser (bins, chunks, linear index)
- [ ] Region -> chunk query (O(log N))
- [ ] BAI builder (`index` from a coordinate-sorted BAM)
- [ ] Region query API (`reader.query("chr1:1000-2000")`)
- [ ] Tests vs samtools-generated indexes

### CLI (minimal) & milestone
- [ ] `tptbio seq view` (BAM/SAM, region, header-only, format output)
- [ ] `tptbio seq index`
- [ ] Benchmark suite vs samtools (view, index, region query) on a human-genome-scale BAM
- [ ] **Milestone:** `tptbio seq view` reads and indexes a standard human genome BAM faster than samtools

---

## Phase 2: Sequence Complete (VCF, BCF, CRAM)

### VCF / BCF
- [ ] VCF header parser (INFO, FORMAT, FILTER, contig, samples)
- [ ] Zero-copy VCF record parser and writer
- [ ] BGZF-compressed VCF
- [ ] BCF2 header and typed value decoding
- [ ] BCF record parser and writer
- [ ] VCF <-> BCF round-trip tests
- [ ] Tabix (TBI) and CSI index parse + build
- [ ] Region queries for VCF/BCF
- [ ] `tptbio seq view` for VCF/BCF
- [ ] Conformance tests vs bcftools

### CRAM / CRAI
- [ ] CRAM container and slice structure parser
- [ ] CRAM header and compression header decoding
- [ ] Core codecs: rANS 4x8, rANS Nx16, Huffman, beta, gamma, subexp, external
- [ ] Gzip, bzip2 and LZMA block compression: find permissive pure-Rust crates for bzip2/LZMA or document unsupported blocks
- [ ] Reference sequence handling (local FASTA + `.fai`, MD5 lookup)
- [ ] Reference-free / embedded reference mode
- [ ] Record decoding (read features -> alignment records)
- [ ] CRAM writer (at least CRAM 3.0)
- [ ] CRAI parse + build, region queries
- [ ] BAM <-> CRAM round-trip tests
- [ ] Conformance tests vs samtools

### FASTA / misc
- [ ] FASTA reader/writer and `.fai` index
- [ ] Sorting utilities (coordinate / name sort)
- [ ] Merge of sorted BAMs
- [ ] Basic stats (flagstat-like)

### Hardening
- [ ] Fuzz targets for BGZF, BAM, VCF, BCF, CRAM, indexes
- [ ] Malformed-input corpus tests (no panics, bounded memory)
- [ ] Benchmarks for VCF/BCF/CRAM vs bcftools/samtools

---

## Phase 3: Vision Core (`tpt-bio-vision`)

### Storage / I/O
- [ ] `ByteSource` trait (local file, mmap, in-memory, HTTP)
- [ ] HTTP range-read backend (permissive client, retry, caching)
- [ ] Own S3 SigV4 request signer (no AWS SDK)
- [ ] Presigned-URL and public-bucket support
- [ ] Block cache for remote reads

### Formats
- [ ] 2D image loading (TIFF/PNG via `image` / `tiff`; 8/16-bit, float)
- [ ] OME-TIFF reader (OME-XML metadata, multi-series, Z/C/T planes, BigTIFF, tiled and striped)
- [ ] Memory-mapped plane access
- [ ] Zarr reader (v2, then v3; chunked, compressed with permissive codecs)
- [ ] HCS plate layout (wells, fields, channels; OME-NGFF plate metadata)
- [ ] ndarray integration for N-d image data
- [ ] Image writer (OME-TIFF and/or Zarr) for results and masks

### SIMD image processing (tpt-simd)
- [ ] 2D convolution kernel (separable and non-separable)
- [ ] 3D convolution kernel (ADR-002)
- [ ] Gaussian blur (2D/3D)
- [ ] Background subtraction (rolling ball / morphological)
- [ ] Median / mean filters
- [ ] Edge detection (Sobel, Scharr, Canny)
- [ ] Morphology (erode, dilate, open, close, top-hat) in 2D/3D
- [ ] Runtime feature detection (AVX2 / AVX-512 / NEON) with scalar fallback
- [ ] Benchmarks vs scalar baseline and reference implementations

### Segmentation
- [ ] Otsu thresholding (global + adaptive/local)
- [ ] Connected-component labeling (2D/3D)
- [ ] Distance transform
- [ ] Marker-controlled watershed
- [ ] Seed detection (local maxima / LoG)
- [ ] Graph-based segmentation (sparse matrix graph cuts via tpt-math)
- [ ] Object filtering (size, border touch), label relabeling
- [ ] Validation against known datasets (e.g. BBBC) with IoU/F1 thresholds

### Feature extraction (500+)
- [ ] Feature framework (per-object, per-channel, parallel over objects)
- [ ] Shape and size features (area, perimeter, eccentricity, solidity, Zernike, etc.)
- [ ] Intensity features (mean, std, quantiles, MAD, integrated, edge intensity)
- [ ] Texture features (Haralick/GLCM, granularity, Gabor)
- [ ] Neighbor / spatial features
- [ ] Colocalization and correlation features
- [ ] Reach 500+ features; document the list vs CellProfiler equivalents
- [ ] Numerical parity tests vs CellProfiler outputs (test-only, tolerance-based)
- [ ] Tabular export (CSV, Parquet or Arrow: license check first)

### CLI & milestone
- [ ] `tptbio vision segment`
- [ ] `tptbio vision features`
- [ ] **Milestone:** `tptbio vision segment` processes a 96-well plate microscopy image in seconds

---

## Phase 4: Statistics & 3D Viewer

### Statistical analysis (tpt-bio-core + tpt-math)
- [ ] Feature normalization (z-score, robust, per-plate)
- [ ] PCA
- [ ] Clustering (k-means, hierarchical, optionally UMAP/t-SNE if feasible)
- [ ] Phenotype grouping and per-well aggregation
- [ ] Quality metrics (Z-prime, SSMD)
- [ ] Tests vs scikit-learn reference results

### 3D viewer (wgpu, ADR-003)
- [ ] Architecture note: share the rendering paradigm with tpt-vertex
- [ ] wgpu (check license expression) volume renderer: texture upload, ray-marching
- [ ] Label/mask overlay and per-cell coloring
- [ ] Mesh / isosurface rendering of segmented cells
- [ ] Large-volume handling (chunked/streamed, level of detail)
- [ ] Browser target (WASM + WebGPU) demo
- [ ] Optional Tauri desktop app (license check for Tauri and its dependency tree)
- [ ] Camera controls, channel and clipping controls
- [ ] Viewer documentation and sample dataset

---

## Phase 5: Integration, CLI, Python

### tpt-bio-pipeline
- [ ] DAG model (nodes, typed ports, edges, cycle detection)
- [ ] Streaming executor with bounded channels and backpressure (crossbeam)
- [ ] Parallel node execution (rayon)
- [ ] Built-in nodes: seq readers/writers/filters, vision load/filter/segment/features, stats
- [ ] Pipeline definition format (TOML/YAML) + validation (license check for serde crates)
- [ ] Error handling, cancellation, progress reporting
- [ ] Checkpointing / resume
- [ ] End-to-end example: high-content screening flow (spec section 7, steps 1-6)
- [ ] Integration tests across seq + vision nodes

### tpt-bio-cli (`tptbio`)
- [ ] clap command tree: `seq`, `vision`, `pipeline run`, `stats`, `view`
- [ ] Shell completions and man pages
- [ ] Consistent output formats, exit codes, logging/verbosity
- [ ] Config file support
- [ ] CLI integration tests
- [ ] Prebuilt binaries via CI (Linux/macOS/Windows)

### tpt-bio-py (PyO3 + maturin)
- [ ] Bindings for seq readers (iterators, region queries)
- [ ] Bindings for vision (load, segment, features)
- [ ] NumPy interop (zero-copy where possible; check numpy crate license)
- [ ] Bindings for pipeline execution
- [ ] Type stubs (`.pyi`), docstrings
- [ ] Release the GIL during heavy work
- [ ] Python tests (pytest)
- [ ] Wheels for Linux (manylinux), macOS (x86_64 + arm64), Windows via maturin in CI
- [ ] Jupyter notebook example: results from the spec section 7 flow

---

## Phase 6: Hardening & Release v0.1.0

### Quality
- [ ] API review across all crates (naming, error types, semver surface)
- [ ] Full rustdoc coverage, mdBook or docs site
- [ ] Tutorials: seq quickstart, vision quickstart, pipeline, Python
- [ ] Performance regression checks in CI
- [ ] Extended fuzzing run, fix findings
- [ ] Memory-safety review of every remaining `unsafe` block
- [ ] Cross-platform test pass (Linux/macOS/Windows, x86_64 + aarch64)
- [ ] `no_std` verification on an embedded target

### License & compliance (final)
- [ ] Full `cargo deny check` clean on the release commit
- [ ] Regenerate `docs/dependency-audit.md` and a third-party notices file
- [ ] Confirm zero-FFI: no `-sys` crates in the dependency tree (`cargo tree` check in CI)
- [ ] Confirm headers/Cargo metadata/LICENSE files consistent in all crates and in the Python package

### Release
- [ ] Choose final crate names and check availability on crates.io
- [ ] Publish order: core -> seq / vision -> pipeline -> cli
- [ ] `cargo publish --dry-run` for all crates
- [ ] Publish the Python package to PyPI
- [ ] Tag `v0.1.0`, GitHub release, release notes
- [ ] Announce (docs site, README badges)

---

## Later / Post-v0.1.0 (parking lot)
- [ ] Additional formats (e.g. GFF/GTF, BED, 2bit)
- [ ] Optional zstd support via a pure-Rust implementation, if one reaches parity
- [ ] .cppipe import (currently out of scope)
- [ ] GPU-accelerated segmentation
- [ ] Cloud-native backends beyond S3 (GCS, Azure) through the ByteSource trait
