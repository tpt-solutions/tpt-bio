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

**Open decisions** (what is still open after the 2026-11-24 review)
- **Poisson:** recommendation is to reimplement the PMF and CDF on `libm` (`lgamma`, incomplete gamma and beta) so the distributions work in `no_std`. Keep `tpt-math-prob-dist` as a dev-dependency oracle in tests. Waiting for a decision.
- **Crates.io publishing:** `tpt-crypto` and `target-lexicon-mit` are git pins, and crates.io rejects git dependencies. Phase 6 needs published versions or vendored copies.
- **TLS for remote I/O:** `tpt-crypto` covers SHA-256 and HMAC but not a TLS client. `tpt-crypto-rustls` is unaudited. Decide before Phase 3 HTTPS.

---

## Phase 0: Project Foundations

### Repo & licensing
- [x] `git init`, `.gitignore`, default branch `master`
- [x] `LICENSE-MIT` (© 2026 TPT Solutions)
- [x] `LICENSE-APACHE` (© 2026 TPT Solutions)
- [x] License headers / `license = "MIT OR Apache-2.0"` in every Cargo.toml _(Cargo `license` field set for every crate; per-file headers not added)_
- [x] `README.md` (overview, license section, "MIT-first" dependency policy)
- [x] `CONTRIBUTING.md` (dual-license contribution clause, dependency policy)
- [x] `SECURITY.md`, `CODE_OF_CONDUCT.md`
- [x] `CHANGELOG.md`

### Workspace scaffold
- [x] Root `Cargo.toml` workspace (shared `[workspace.package]`, `[workspace.dependencies]`, lints)
- [x] `rust-toolchain.toml` and pinned MSRV
- [x] `rustfmt.toml`, `clippy.toml`
- [x] Scaffold `tpt-bio-core` (`#![no_std]` + `alloc`)
- [x] Scaffold `tpt-bio-seq`
- [x] Scaffold `tpt-bio-vision` _(empty crate; no code yet)_
- [x] Scaffold `tpt-bio-pipeline` _(contains the DAG model)_
- [x] Scaffold `tpt-bio-cli` (binary `tptbio`)
- [x] Scaffold `bindings/tpt-bio-py` (PyO3 + maturin `pyproject.toml`) _(separate workspace; see open decision on pyo3)_
- [ ] Confirm tpt-math and tpt-simd versions and licenses (including `no_std` support) _(tpt-simd 0.1.1 is adopted and needs Rust 1.95, which is now the MSRV; its `no_std` features are not verified)_
- [x] `unsafe` policy: `#![forbid(unsafe_code)]` where possible; document and isolate any exceptions (mmap, SIMD)

### Dependency license audit
- [x] Verify license expression of each planned dependency: tpt-math, tpt-simd, tpt-crypto, num-traits, flate2 (+ miniz_oxide), memmap2, rayon, crossbeam, image, tiff, ndarray, clap, pyo3, wgpu, ureq, serde/toml _(see docs/dependency-audit.md)_
- [x] Reject or replace anything Apache-2.0-only, MPL, or copyleft (including transitive) _(arrow/parquet rejected; target-lexicon via pyo3 awaits a decision)_
- [ ] Use `tpt-crypto` instead of `ring` for all crypto needs (SigV4 HMAC/SHA-256, HTTPS/TLS client, checksums); drop the `ring` clarify block from `deny.toml` _(SigV4 uses `tpt-crypto-hash`, pinned to a git commit; TLS is still open)_
- [x] Verify tpt-crypto's license, `no_std` support and crates.io version; make sure the HTTP client's TLS stack can use it (or document the gap) _(MIT OR Apache-2.0, `tpt-crypto-hash` is `no_std`; not on crates.io, so it is pinned to a git commit; TLS not covered)_
- [x] Record the audit result in `docs/dependency-audit.md`

### cargo-deny & CI
- [x] `deny.toml` with the corrected policy (MIT-option expressions only, copyleft denied, advisories, bans, sources)
- [x] Verify `cargo deny check` fails on a deliberately Apache-only test crate, then remove it _(verified with a temporary path crate, since removed)_
- [x] GitHub Actions: fmt, clippy `-D warnings`, test on Linux/macOS/Windows _(workflow written; not yet run on GitHub)_
- [x] CI: `cargo deny check`
- [x] CI: `no_std` build of `tpt-bio-core` (e.g. `thumbv7em-none-eabihf`)
- [x] CI: MSRV job
- [x] CI: docs build (`cargo doc`, deny warnings)
- [x] CI: dependency caching and nightly scheduled run
- [x] Dependabot / Renovate config _(Dependabot)_

### Test data & tooling
- [ ] Choose permissively-licensed test data (small BAM/CRAM/VCF/FASTQ, OME-TIFF, Zarr) _(tests are self-contained with generated data; external fixtures not chosen, see docs/test-data.md)_
- [ ] Set up `criterion` (check license) or a custom benchmark harness _(not added yet)_
- [x] Set up `cargo-fuzz` infrastructure _(fuzz/ builds; targets not yet run under libFuzzer)_
- [ ] Install samtools/bcftools locally for oracle comparisons (test-only, never linked) _(BLOCKED: not installed on this machine)_

---

## Phase 1: Sequence Core (`tpt-bio-core` + `tpt-bio-seq` basics)

### tpt-bio-core (`no_std` + alloc)
- [x] Nucleotide types (DNA/RNA, IUPAC codes, 2-bit and 4-bit encodings)
- [x] Amino acid type and codon table
- [x] Phred quality score type and conversions
- [x] Genomic coordinate / region types (0-based vs 1-based clarity)
- [x] Reverse complement, k-mer iteration
- [x] Poisson distribution (PMF, CDF, quantile) via tpt-math _(PMF, CDF, quantile; requires the std feature)_
- [x] Binomial distribution (PMF, CDF, quantile) via tpt-math _(PMF, CDF, quantile (quantile by linear search); requires the std feature)_
- [x] Basic bio-math (GC content, entropy, etc.)
- [x] Parallel iterator helpers (`std` feature, rayon-based)
- [x] `no_std` and `std` feature matrix, tested in CI
- [ ] Unit and property tests; rustdoc with examples _(unit tests and doctests in place; property-based tests not yet added)_

### BGZF (flate2, pure Rust)
- [x] BGZF block reader (header, extra field, CRC32, ISIZE validation)
- [x] EOF marker handling
- [x] Parallel block decompression (rayon/crossbeam)
- [x] BGZF writer with parallel block compression _(parallel via compress_parallel; the streaming writer is sequential)_
- [x] Virtual file offsets (block offset + uncompressed offset)
- [x] Seek by virtual offset
- [x] Mmap-backed reader (memmap2)
- [x] Corruption and truncation tests, fuzz target

### FASTQ
- [x] Zero-copy FASTQ record parser
- [x] Gzip / BGZF FASTQ input
- [x] Paired-end / interleaved handling
- [x] FASTQ writer
- [x] Streaming parallel batch reader _(bounded-channel reader thread; per-batch work via rayon)_
- [x] Tests and fuzz target

### SAM / BAM
- [ ] SAM header parser (@HD, @SQ, @RG, @PG, @CO) _(@SQ is parsed; other header lines are kept as text, without structured accessors)_
- [x] SAM record parser (text) and writer
- [x] BAM header and reference dictionary
- [x] Zero-copy BAM record view (cigar, seq, qual, aux tags)
- [x] CIGAR operations and alignment-end calculation
- [x] Auxiliary tag parsing (all types, arrays)
- [x] BAM writer
- [x] SAM <-> BAM round-trip tests
- [ ] Conformance tests vs samtools and noodles outputs _(BLOCKED: samtools not installed; noodles not used)_

### BAI index & random access
- [x] BAI parser (bins, chunks, linear index)
- [x] Region -> chunk query (O(log N))
- [x] BAI builder (`index` from a coordinate-sorted BAM)
- [x] Region query API (`reader.query("chr1:1000-2000")`) _(fetch() takes a reference ID and coordinates; no string API)_
- [ ] Tests vs samtools-generated indexes _(BLOCKED: samtools not installed)_

### CLI (minimal) & milestone
- [x] `tptbio seq view` (BAM/SAM, region, header-only, format output)
- [x] `tptbio seq index`
- [ ] Benchmark suite vs samtools (view, index, region query) on a human-genome-scale BAM _(BLOCKED: samtools not installed; no benchmark harness yet)_
- [ ] **Milestone:** `tptbio seq view` reads and indexes a standard human genome BAM faster than samtools _(not measured; see benchmark item)_

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
- [x] FASTA reader/writer and `.fai` index
- [x] Sorting utilities (coordinate / name sort)
- [x] Merge of sorted BAMs _(merge over record streams, so any sorted source works)_
- [x] Basic stats (flagstat-like)

### Hardening
- [ ] Fuzz targets for BGZF, BAM, VCF, BCF, CRAM, indexes
- [ ] Malformed-input corpus tests (no panics, bounded memory)
- [ ] Benchmarks for VCF/BCF/CRAM vs bcftools/samtools

---

## Phase 3: Vision Core (`tpt-bio-vision`)

### Storage / I/O
- [ ] `ByteSource` trait (local file, mmap, in-memory, HTTP)
- [ ] HTTP range-read backend (permissive client, retry, caching)
- [x] Own S3 SigV4 request signer (no AWS SDK) _(signer only; checked against the AWS ListUsers example vector)_
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
- [ ] 2D convolution kernel (separable and non-separable) _(separable done via `tpt-simd`; non-separable not started)_
- [ ] 3D convolution kernel (ADR-002)
- [ ] Gaussian blur (2D/3D) _(2D done via the separable kernel; 3D not started)_
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
- [x] Tabular export: CSV and JSON Lines through `tpt-strata` and `tptbio seq export`. Parquet is out of scope, so Arrow is not used.

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
- [ ] DAG model (nodes, typed ports, edges, cycle detection) _(nodes, edges and cycle detection done; typed ports not implemented)_
- [ ] Streaming executor with bounded channels and backpressure (crossbeam) _(not implemented; the executor is in-memory)_
- [x] Parallel node execution (rayon) _(level-parallel execution)_
- [ ] Built-in nodes: seq readers/writers/filters, vision load/filter/segment/features, stats
- [ ] Pipeline definition format (TOML/YAML) + validation (license check for serde crates)
- [ ] Error handling, cancellation, progress reporting _(errors and cancellation done; progress reporting not)_
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
- [x] Type stubs (`.pyi`), docstrings
- [x] Release the GIL during heavy work _(sequence utilities only)_
- [x] Python tests (pytest) _(5 tests, run against a locally built module)_
- [ ] Wheels for Linux (manylinux), macOS (x86_64 + arm64), Windows via maturin in CI _(maturin not run in CI yet)_
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
