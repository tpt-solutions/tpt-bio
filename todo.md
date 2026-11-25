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
- **Poisson:** recommendation is to reimplement the PMF and CDF on `libm` (`lgamma`, incomplete gamma and beta) so the distributions work in `no_std`. Keep `tpt-math-prob-dist` as a dev-dependency oracle in tests. **Decided: recommendation followed.** Implemented on `libm`; no oracle crate is available offline, so tests use hand-computed values.
- **Crates.io publishing:** `tpt-crypto` and `target-lexicon-mit` are git pins, and crates.io rejects git dependencies. Phase 6 needs published versions or vendored copies.
- **wgpu and the Apache-2.0-only rule (decided 2026-11-25):** wgpu is rejected. Its tree has `codespan-reporting`, a non-optional naga dependency that is Apache-2.0 only, so no feature choice removes it. `gl_generator`, `glutin_wgl_sys` (Apache-2.0 only), `hexf-parse` (CC0-1.0) and `paste` (RUSTSEC-2024-0436) also fail the policy. wgpu-based native GPU rendering is dropped. Browser WebGPU goes through `web-sys` and `wasm-bindgen` (MIT OR Apache-2.0) instead. **Native follow-up (2026-11-25): Vulkan through `ash` with an offline-compiled SPIR-V shader, behind the `RenderBackend` trait.** `ash` passed the licence audit; the webview route (`wry`) failed it. `tpt-gpu` was considered and skipped (no Vulkan or D3D12, no WGSL or SPIR-V input, no 3D textures, pre-1.0).
- **TLS for remote I/O:** `tpt-crypto` covers SHA-256 and HMAC but not a TLS client. `tpt-crypto-rustls` is unaudited. **Decided for now: no TLS.** Remote I/O is plain `http://` only, and `https://` is rejected. SigV4 over `http://` sends signatures in cleartext, so this is not safe on untrusted networks.

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
- [x] Confirm tpt-math and tpt-simd versions and licenses (including `no_std` support) _(tpt-simd 0.1.1 is `#![no_std]` with `alloc` and `std` as opt-in features, MIT OR Apache-2.0, and needs Rust 1.95, which is now the MSRV; licenses also pass `cargo deny`)_
- [x] `unsafe` policy: `#![forbid(unsafe_code)]` where possible; document and isolate any exceptions (mmap, SIMD)

### Dependency license audit
- [x] Verify license expression of each planned dependency: tpt-math, tpt-simd, tpt-crypto, num-traits, flate2 (+ miniz_oxide), memmap2, rayon, crossbeam, image, tiff, ndarray, clap, pyo3, wgpu, ureq, serde/toml _(see docs/dependency-audit.md)_
- [x] Reject or replace anything Apache-2.0-only, MPL, or copyleft (including transitive) _(arrow/parquet rejected; target-lexicon via pyo3 awaits a decision)_
- [ ] Use `tpt-crypto` instead of `ring` for all crypto needs (SigV4 HMAC/SHA-256, HTTPS/TLS client, checksums); drop the `ring` clarify block from `deny.toml` _(SigV4 uses `tpt-crypto-hash`, pinned to a git commit; TLS is still open. 2026-11-26: `ring` is not in Cargo.lock, and `deny.toml` only keeps it as a ban. No clarify block is left to drop)_
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
- [ ] Choose permissively-licensed test data (small BAM/CRAM/VCF/FASTQ, OME-TIFF, Zarr) _(decided, see docs/test-data.md: generated in code, or written by reference tools and committed with provenance; no downloads. Still open: the human-genome-scale BAM, and BBBC images, whose host did not answer)_
- [x] Set up the benchmark harness: `tpt-benchmark` (git dev-dependency, pinned to `8384f44`) replaces `criterion`, which is removed _(its Apache-only transitive crates are no longer in the lock; `borsh-derive` is a stale lock entry, see docs/dependency-audit.md)_
- [x] Set up `cargo-fuzz` infrastructure _(fuzz/ builds; targets not yet run under libFuzzer)_
- [x] Install samtools/bcftools locally for oracle comparisons (test-only, never linked) _(installed in WSL Ubuntu: samtools 1.22.1, bcftools 1.22, tabix 1.22.1; not on the Windows PATH, and never linked into the crates)_

---

## Phase 1: Sequence Core (`tpt-bio-core` + `tpt-bio-seq` basics)

### tpt-bio-core (`no_std` + alloc)
- [x] Nucleotide types (DNA/RNA, IUPAC codes, 2-bit and 4-bit encodings)
- [x] Amino acid type and codon table
- [x] Phred quality score type and conversions
- [x] Genomic coordinate / region types (0-based vs 1-based clarity)
- [x] Reverse complement, k-mer iteration
- [x] Poisson distribution (PMF, CDF, quantile) via tpt-math _(PMF, CDF, quantile on libm; works under no_std; hand-computed tests)_
- [x] Binomial distribution (PMF, CDF, quantile) via tpt-math _(PMF, CDF, quantile (quantile by linear search); requires the std feature)_
- [x] Basic bio-math (GC content, entropy, etc.)
- [x] Parallel iterator helpers (`std` feature, rayon-based)
- [x] `no_std` and `std` feature matrix, tested in CI
- [x] Unit and property tests; rustdoc with examples _(proptest properties for nucleotide, k-mer and coordinate code; doctests for public types)_

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
- [x] SAM header parser (@HD, @SQ, @RG, @PG, @CO) _(structured accessors: `hd()`, `read_groups()`, `programs()`, `comments()`; the raw text is still kept verbatim)_
- [x] SAM record parser (text) and writer
- [x] BAM header and reference dictionary
- [x] Zero-copy BAM record view (cigar, seq, qual, aux tags)
- [x] CIGAR operations and alignment-end calculation
- [x] Auxiliary tag parsing (all types, arrays)
- [x] BAM writer
- [x] SAM <-> BAM round-trip tests
- [x] Conformance tests vs samtools and noodles outputs _(samtools 1.22 only: our BAM is read by samtools with identical records; samtools BAM decodes to its own records; noodles not used)_

### BAI index & random access
- [x] BAI parser (bins, chunks, linear index)
- [x] Region -> chunk query (O(log N))
- [x] BAI builder (`index` from a coordinate-sorted BAM)
- [x] Region query API (`reader.query("chr1:1000-2000")`) _(fetch() takes a reference ID and coordinates; no string API)_
- [x] Tests vs samtools-generated indexes _(our BAI matches `samtools index` for a 250-record file byte for byte, and for 20,000 and 300,000 reads once parsed; only bin order in the file differs, as samtools uses hash order; region queries match samtools)_

### CLI (minimal) & milestone
- [x] `tptbio seq view` (BAM/SAM, region, header-only, format output)
- [x] `tptbio seq index`
- [ ] Benchmark suite vs samtools (view, index, region query) on a human-genome-scale BAM _(synthetic 300,000-read BAM: view 304 ms vs samtools 298 ms, index 94 ms vs 273 ms, region 51 ms vs 14 ms; samtools' medians varied from 269 to 461 ms, and the environments differ (WSL vs native); no human genome data, so the milestone stays open; see docs/samtools-bcftools-conformance.md)_ 2026-11-26 update: run on a real human chromosome 20 BAM (318 MB, 2.9 M records, not genome-scale). Output matches samtools record for record. tptbio is faster for index (0.41x) and region query (0.16x) and slower for full view (1.39x). See docs/samtools-bcftools-conformance.md. The whole-genome BAM (15.6 GB) would take about 9.5 h at the measured 0.46 MB/s, so it is not run.
- [ ] **Milestone:** `tptbio seq view` reads and indexes a standard human genome BAM faster than samtools _(not measured; see benchmark item)_ 2026-11-26: chr20 result is mixed (view slower, index and region faster). Still open until a whole-genome BAM is measured.

---

## Phase 2: Sequence Complete (VCF, BCF, CRAM)

### VCF / BCF
- [x] VCF header parser (INFO, FORMAT, FILTER, contig, samples) _(`tpt-bio-seq::vcf`; structured meta lines are split into fields)_
- [x] Zero-copy VCF record parser and writer _(`VcfRecord` borrows its line; `VcfWriter` validates columns; line-numbered errors)_
- [x] BGZF-compressed VCF _(`VcfReader::new_bgzf`, `VcfWriter::new_bgzf` + `finish`)_
- [x] BCF2 header and typed value decoding _(checked against bcftools 1.22: bcftools reads our BCF with identical records, and we read its BCF identically; an empty ID was a bug, now fixed)_
- [x] BCF record parser and writer _(records decode to VCF text, so `next_record` returns a `VcfRecord`)_
- [x] VCF <-> BCF round-trip tests _(both directions checked against bcftools 1.22 on 240 generated records)_
- [x] Tabix (TBI) and CSI index parse + build _(TBI and CSI are byte-identical to tabix and bcftools once decompressed; fixed: the trailing n_no_coor count, and CSI depth 0 was rejected)_
- [x] Region queries for VCF/BCF _(match tabix and bcftools on seven regions, including empty and edge regions)_
- [x] `tptbio seq view` for VCF/BCF _(VCF, bgzipped VCF and BCF, with region queries via tabix or CSI)_
- [x] Conformance tests vs bcftools _(BCF, TBI and CSI checked against bcftools 1.22 and tabix 1.22 on generated data; not exhaustive over every BCF type)_

### CRAM / CRAI
- [x] CRAM container and slice structure parser _(container, block framing and slice header; EOF marker is a self-consistent substitute, not the spec's)_
- [x] CRAM header and compression header decoding _(compression header and encodings decoded; field layouts from memory, not verified against samtools)_
- [ ] Core codecs: rANS 4x8, rANS Nx16, Huffman, beta, gamma, subexp, external _(all decode; rANS Nx16 order 1 and 32 states verified against samtools; RLE and STRIPE only from htscodecs-derived vectors; bzip2 decoded since 2026-11-26 (pure Rust, `cram/bzip2.rs`, verified against samtools' bzip2 output). Still open: order 2, which samtools does not write)_
- [x] Gzip, bzip2 and LZMA block compression _(gzip, rANS and bzip2 blocks decode. bzip2 decoded 2026-11-26 in pure Rust and checked against samtools 1.22 `use_bzip2=1` output. LZMA still returns Unsupported, documented)_
- [x] Reference sequence handling (local FASTA + `.fai`, MD5 lookup) _(FASTA reference with MD5 region check; in-memory FASTA also accepted)_
- [x] Reference-free / embedded reference mode _(reads needing copied bases return Unsupported; unknown-bases flag not handled)_
- [x] Record decoding (read features -> alignment records) _(all standard read features except RG; spec bit meanings not checked against samtools)_
- [x] CRAM writer (at least CRAM 3.0) _(CRAM 3.0 writer; gzip data blocks; EOF marker not the spec's 38 bytes)_
- [x] CRAI parse + build, region queries _(build, parse and query, tested against a brute-force scan)_
- [x] BAM <-> CRAM round-trip tests _(round trips for multi-slice, reference and reference-free files; no samtools comparison)_
- [ ] Conformance tests vs samtools _(CRAM 3.0 and 3.1 decode both ways with samtools 1.22, every column and tag matching; TOK3 name codec verified for rANS, arithmetic and duplicate/zero-delta name shapes; NOP is never produced by samtools; RLE and STRIPE only from reference vectors; TOK3 writing not implemented; see docs/cram-conformance.md)_

### FASTA / misc
- [x] FASTA reader/writer and `.fai` index
- [x] Sorting utilities (coordinate / name sort)
- [x] Merge of sorted BAMs _(merge over record streams, so any sorted source works)_
- [x] Basic stats (flagstat-like)

### Hardening
- [x] Fuzz targets for BGZF, BAM, VCF, BCF, CRAM, indexes _(targets for BGZF, BAM, BAI, FASTQ, SAM, VCF, BCF, tabix, FASTA and CRAM)_
- [x] Malformed-input corpus tests (no panics, bounded memory) _(12 parsers, no panics found; zero-filled over-allocations fixed; bounded memory not asserted by tests)_
- [ ] Benchmarks for VCF/BCF/CRAM vs bcftools/samtools _(300,000-record medians on one machine: VCF.gz to text 114 ms vs bcftools 165 ms; BCF to text 165 ms vs 111 ms; CRAM 3.1 to SAM text 819 ms vs samtools 150 ms (reference-free); BAM from the earlier round; harness is tpt-benchmark pinned to 8384f44, and the timings are not a comparable native-vs-native run)_

---

## Phase 3: Vision Core (`tpt-bio-vision`)

### Storage / I/O
- [x] `ByteSource` trait (local file, mmap, in-memory, HTTP) _(trait with in-memory and local-file sources in `tpt-bio-vision::source`; `ByteSource` is implemented for `tpt_bio_seq::mmap::MappedFile`, so vision adds no `unsafe`; the caller opens the map with the unsafe `MappedFile::open` under its contract. Checked by a scratch binary that mapped a 4 MB file and read ranges through the trait, and by tests on the safe `MappedFile::from_vec`. HTTP is a separate item below)_
- [x] HTTP range-read backend (permissive client, retry, caching) _(plain http only; https is rejected, since no TLS is used)_
- [x] Own S3 SigV4 request signer (no AWS SDK) _(signer only; checked against the AWS ListUsers example vector)_
- [x] Presigned-URL and public-bucket support
- [x] Block cache for remote reads _(LRU, 64 KiB default blocks)_

### Formats
- [x] 2D image loading (TIFF/PNG via `image` / `tiff`; 8/16-bit, float) _(PNG and TIFF; greyscale and per-channel planes)_
- [x] OME-TIFF reader (OME-XML metadata, multi-series, Z/C/T planes, BigTIFF, tiled and striped) _(OME-XML is hand-written in tests, not Bio-Formats output; multi-file OME not supported)_
- [x] Memory-mapped plane access _(MappedTiff is generic over its backing; mapping a file needs unsafe in tpt-bio-seq)_
- [x] Zarr reader (v2, then v3; chunked, compressed with permissive codecs) _(v2 and v3 read with none and gzip; blosc and zstd return Unsupported)_
- [x] HCS plate layout (wells, fields, channels; OME-NGFF plate metadata) _(parses NGFF 0.4 `.zattrs` and 0.5 `zarr.json`; `validate` checks the spec's MUST rules; `ImageChannels` reads each field's channel axis and omero labels, and `check_same_channels` compares fields. Verified against plates written by the reference `ome-zarr` 0.21.0 writer, fixtures in `tests/fixtures/hcs`; 16 tests, 2026-11-25)_
- [x] ndarray integration for N-d image data _(Array2 and Array3 conversion)_
- [x] Image writer (OME-TIFF and/or Zarr) for results and masks _(TIFF, PNG, OME-TIFF and Zarr v2 writers)_

### SIMD image processing (tpt-simd)
- [x] 2D convolution kernel (separable and non-separable) _(separable and general square kernels)_
- [x] 3D convolution kernel (ADR-002) _(separable 3D, with boundary modes)_
- [x] Gaussian blur (2D/3D)
- [x] Background subtraction (rolling ball / morphological)
- [x] Median / mean filters
- [x] Edge detection (Sobel, Scharr, Canny) _(Canny has no built-in smoothing)_
- [x] Morphology (erode, dilate, open, close, top-hat) in 2D/3D
- [ ] Runtime feature detection (AVX2 / AVX-512 / NEON) with scalar fallback _(decided 2026-11-26: deferred past v0.1.0. Runtime dispatch needs `#[target_feature]` functions, which are `unsafe` and clash with `tpt-bio-vision`'s `forbid(unsafe_code)`. Per-ISA builds are the supported route. Reopened 2026-11-25. `tpt-simd` chooses its vector kernels at compile time, from target features, so runtime detection does not change the kernels. `simd::backend()` now reports the compiled kernels, and `detected_features()` reports the CPU. The default build stays baseline x86_64 for portability. Still open: true runtime dispatch, which needs a change in tpt-simd or per-ISA builds)_
- [ ] Benchmarks vs scalar baseline and reference implementations _(`examples/simd_bench.rs`, a 2048x2048 separable Gaussian blur, sigma 2, median of 7. Baseline build: about 117 ms. Built with `RUSTFLAGS="-C target-feature=+avx,+fma"`: about 22 ms, roughly 5 times faster. Still open: a reference implementation to compare against, and the other filters)_

### Segmentation
- [x] Otsu thresholding (global + adaptive/local) _(global, tile-based and sliding-window Otsu)_
- [x] Connected-component labeling (2D/3D)
- [x] Distance transform _(exact Euclidean, 2D and 3D)_
- [x] Marker-controlled watershed
- [x] Seed detection (local maxima / LoG)
- [x] Graph-based segmentation (sparse matrix graph cuts via tpt-math) _(two-label min-cut with own sparse Dinic graph; tpt-math not used)_
- [x] Object filtering (size, border touch), label relabeling
- [ ] Validation against known datasets (e.g. BBBC) with IoU/F1 thresholds _(IoU/F1 matching is implemented and tested on hand-computed cases; no BBBC data was available)_

### Feature extraction (500+)
- [x] Feature framework (per-object, per-channel, parallel over objects) _(serial only, since rayon is not a dependency of tpt-bio-vision; per-channel is not done)_
- [x] Shape and size features (area, perimeter, eccentricity, solidity, Zernike, etc.) _(2D area, crack perimeter, eccentricity, solidity, Zernike; 3D volume, surface count, axes)_
- [x] Intensity features (mean, std, quantiles, MAD, integrated, edge intensity) _(mean, std, quantiles, MAD, integrated and edge intensity)_
- [x] Texture features (Haralick/GLCM, granularity, Gabor) _(Haralick over four directions, granularity spectrum, Gabor bank; not compared with mahotas)_
- [x] Neighbor / spatial features _(touching-neighbour count and nearest-centroid distance)_
- [x] Colocalization and correlation features _(Pearson, Manders M1 and M2 per object)_
- [ ] Reach 500+ features; document the list vs CellProfiler equivalents _(a few dozen features exist; the 500+ target is far off)_
- [ ] Numerical parity tests vs CellProfiler outputs (test-only, tolerance-based)
- [x] Tabular export: CSV and JSON Lines through `tpt-strata` and `tptbio seq export`. Parquet is out of scope, so Arrow is not used.

### CLI & milestone
- [x] `tptbio vision segment` _(Otsu, labelling and size filtering; writes a 16-bit TIFF)_
- [x] `tptbio vision features` _(shape and intensity columns; takes labels and an optional intensity image)_
- [ ] **Milestone:** `tptbio vision segment` processes a 96-well plate microscopy image in seconds

---

## Phase 4: Statistics & 3D Viewer

### Statistical analysis (tpt-bio-core + tpt-math)
- [x] Feature normalization (z-score, robust, per-plate) _(`tpt-bio-core::stats::normalize`)_
- [x] PCA _(Jacobi eigendecomposition, no external crate; `stats::pca`)_
- [x] Clustering (k-means, hierarchical, optionally UMAP/t-SNE if feasible) _(k-means, agglomerative, t-SNE and UMAP; simplified and unverified against references)_
- [x] Phenotype grouping and per-well aggregation _(grouping and mean/median aggregation; no threshold classifier, since none was specified)_
- [x] Quality metrics (Z-prime, SSMD) _(`stats::qc`)_
- [x] Tests vs scikit-learn reference results _(2026-11-26, scikit-learn 1.9.0 now installed. `tpt-bio-core/tests/sklearn_reference.rs` checks z-score, robust scaling, PCA (mean, variance, ratio, axes up to sign, scores) and k-means (partition and inertia) to 1e-9. The generator is `tpt-bio-core/tests/fixtures/sklearn/generate.py`. Not covered: t-SNE, UMAP, agglomerative clustering and the QC metrics)_

### 3D viewer (ADR-003, wgpu superseded)
- [x] Architecture note: share the rendering paradigm with tpt-vertex _(docs/architecture-3d-viewer.md; tpt-vertex not available, so the comparison is unverified)_
- [x] ~~wgpu volume renderer~~ _(wgpu rejected by the licence policy, see the open decisions. Closed as won't-do. The WGSL shader is `tpt-bio-vision/src/render3d/raycast.wgsl`, used by the browser path, and the retired wgpu host code is in `gpu-draft/wgpu_host.rs.txt`)_
- [x] Native render backend seam _(`render3d::backend`: `RenderBackend` trait, `CpuBackend`. `render_volume` and `render_labelled` now run rows in parallel with rayon, output identical to the serial reference)_
- [x] Licence audit of native GPU candidates _(2026-11-25, docs/dependency-audit.md: `ash` 0.38 + `libloading` pass; `wry` fails on `dpi` (Apache-2.0 AND MIT) and `dunce` (MIT-0 not allowed) and has a `*-sys` crate)_
- [ ] Native Vulkan backend: `ash` + offline-compiled SPIR-V ray caster implementing `RenderBackend` _(chosen 2026-11-25. Decided 2026-11-26: `libloading` is accepted. It is dlopen of the system loader, not a `-sys` crate, so it passes the zero-FFI CI rule. The backend goes in a new opt-in crate, `tpt-bio-vulkan`, because `tpt-bio-vision` forbids unsafe code. BLOCKED 2026-11-26: no WGSL-to-SPIR-V compiler is installed, and installing naga-cli was denied. Next step: provide `glslc`, `naga` or a Vulkan SDK, or commit the `.spv` built elsewhere. Then the shared CPU-side prep (rays, LUT, packing) moves out of `webgpu.rs`, and the backend is checked against `CpuBackend` on the Intel GPU here. Windows and Linux; Metal is separate)_
- [x] Label/mask overlay and per-cell coloring _(CPU renderer, first-hit blending)_
- [x] Mesh / isosurface rendering of segmented cells _(marching tetrahedra and a flat-shaded CPU rasteriser)_
- [x] Large-volume handling (chunked/streamed, level of detail) _(block-max LOD with the usual rule: coarsest level whose voxels cover at most N pixels. `render3d::stream` builds levels from a `SlabSource` one z-slab at a time, equal to the in-memory result. `stream_level_for_camera` streams only the level the camera needs. `OmeSlabs` is a file-backed source over OME-TIFF planes, verified on single- and multi-channel files. `ZarrSlabs` reads `[z, y, x]` and `[c, z, y, x]` Zarr arrays through the new `ZarrArray::read_block`, which fetches only the chunks a slab overlaps; a test counts the chunk reads. The interactive viewer is not built, as with the camera controls item)_
- [ ] Browser target (WASM + WebGPU) demo _(2026-11-26: the page loads over HTTP and the WASM module is fetched, but headless Edge reports `requestAdapter returned null`, with and without WebGPU flags and SwiftShader. Still unverified, so it needs a machine with a WebGPU adapter. Earlier notes: built: `tpt-bio-web` compiles for wasm32, `tpt-bio-vision` feature `webgpu` uses `web-sys` WebGPU with the same shader and packing as the CPU reference, and `wasm-bindgen` 0.2.130 generated the bindings. NOT yet run in a browser: the page needs an HTTP server, and starting one was denied, so the GPU output is unverified. Next step: serve `tpt-bio-web/www` and check `run_comparison` reports PASS in Chrome)_
- [x] ~~Optional Tauri desktop app~~ _(won't do, decided 2026-11-26: Tauri depends on `wry`, which failed the licence audit on `dpi`, `dunce` and a `*-sys` crate. See docs/dependency-audit.md)_
- [x] Camera controls, channel and clipping controls _(library API only; no interactive UI)_
- [x] Viewer documentation and sample dataset _(docs/viewer.md and a sample-volume generator example)_

---

## Phase 5: Integration, CLI, Python

### tpt-bio-pipeline
- [x] DAG model (nodes, typed ports, edges, cycle detection) _(`TypedDag` rejects mismatched types at `add_edge`)_
- [x] Streaming executor with bounded channels and backpressure (crossbeam) _(`stream::stream_map`; with more than one worker the reorder buffer is unbounded, and the reported error is the first to arrive)_
- [x] Parallel node execution (rayon) _(level-parallel execution)_
- [x] Built-in nodes: seq readers/writers/filters, vision load/filter/segment/features, stats _(seq read, count and CSV write; vision load, filter, Otsu, label, measure; stats nodes)_
- [x] Pipeline definition format (TOML/YAML) + validation (license check for serde crates) _(TOML via serde and toml; explicit edge table not implemented)_
- [x] Error handling, cancellation, progress reporting _(`Progress` trait with `execute_observed`)_
- [x] Checkpointing / resume _(`execute_checkpointed`; checkpoints are keyed by node name only, with no invalidation when code or data changes)_
- [x] End-to-end example: high-content screening flow _(spec section 7 steps 1 to 6 in 2D, on generated data; no 3D, OME-TIFF or 500+ features)_
- [x] Integration tests across seq + vision nodes _(BAM count and PNG threshold, label and measure pipelines)_

### tpt-bio-cli (`tptbio`)
- [x] clap command tree: `seq`, `vision`, `pipeline run`, `stats`, `view` _(`vision` and `pipeline run` are stubs)_
- [x] Shell completions and man pages _(bash, zsh, fish, powershell, elvish; one man page per command)_
- [ ] Consistent output formats, exit codes, logging/verbosity _(exit codes 0/1/2 are documented in the CLI README and tested; `seq export` now renders `table` as aligned columns, as the other row commands do; logging is still limited to -q/-v and log levels in the config file)_
- [x] Config file support _(tptbio.toml with threads (recorded only), log_level, output_format)_
- [x] CLI integration tests _(integration tests run against the built binary)_
- [x] Prebuilt binaries via CI (Linux/macOS/Windows) _(2026-11-26: `.github/workflows/cli-binaries.yml` builds `tptbio` for x86_64 Linux, x86_64 and aarch64 macOS, and x86_64 Windows, on manual dispatch. It uploads workflow artifacts only, with no release. Not yet run on GitHub. Replayed locally on Windows: release build, and `tptbio --version` prints `tptbio 0.1.0`)_

### tpt-bio-py (PyO3 + maturin)
- [x] Bindings for seq readers (iterators, region queries) _(FASTQ iterator and BAM region queries)_
- [x] Bindings for vision (load, segment, features) _(image loading, threshold, labelling, watershed and per-object features)_
- [x] NumPy interop (zero-copy where possible; check numpy crate license) _(numpy 0.30 (BSD-2-Clause); output zero-copy, input copies once)_
- [x] Bindings for pipeline execution _(Pipeline.from_toml(path).run())_
- [x] Type stubs (`.pyi`), docstrings
- [x] Release the GIL during heavy work _(sequence utilities only)_
- [x] Python tests (pytest) _(5 tests, run against a locally built module)_
- [ ] Wheels for Linux (manylinux), macOS (x86_64 + arm64), Windows via maturin in CI _(wheels.yml builds artifacts on manual dispatch only, not run; a Windows wheel built locally)_
- [x] Jupyter notebook example: results from the spec section 7 flow _(`examples/spec7_flow.ipynb` executes end to end with nbclient on a locally built Windows wheel, 2 s, outputs saved; needs numpy, which the package does not declare)_

---

## Phase 6: Hardening & Release v0.1.0

### Quality
- [x] API review across all crates (naming, error types, semver surface) _(docs/api-review.md with 39 findings; not all fixed)_
- [x] Full rustdoc coverage, mdBook or docs site _(rustdoc passes with -D warnings for four crates; `missing_docs` reports 0 undocumented public items in all five library crates, checked 2026-11-25; no mdBook site, so the "or" branch is the docs-coverage half only)_
- [x] Tutorials: seq quickstart, vision quickstart, pipeline, Python _(docs/tutorials, snippets run as examples; Python via a locally built wheel)_
- [x] Performance regression checks in CI _(.github/workflows/perf.yml runs weekly and on demand, informational only; not run on GitHub)_
- [x] Extended fuzzing run, fix findings _(2026-11-25, nightly cargo-fuzz in WSL: each of the ten parser targets ran for about 500 s, bgzf_decompress 500 s twice. One real bug found and fixed: a BGZF block whose BSIZE was smaller than its header and footer panicked (c338925). Every saved artifact now replays cleanly. Open: a 42-byte CRAM input declares a 67 MB rANS output, within the 256 MiB per-block cap. Release decode takes about 0.2 s, so the fuzz build's 14 s is debug assertions. A limit across all blocks in one file is not yet implemented)_
- [x] Memory-safety review of every remaining `unsafe` block _(one unsafe block in tpt-bio-seq mmap; MappedFile::open is now unsafe)_
- [ ] Cross-platform test pass (Linux/macOS/Windows, x86_64 + aarch64) _(Windows x86_64 passes locally, and `cargo test --workspace` passes with `-D warnings` on 2026-11-26. Linux x86_64 in WSL passed the full default-member workspace test run, 23 suites, 2026-11-25, with a separate target directory. ci.yml's test matrix (ubuntu, macos, windows) and the aarch64-macOS build in cli-binaries.yml exist but have not run on GitHub. macOS and aarch64 are not run here)_
- [x] `no_std` verification on an embedded target _(tpt-bio-core builds for thumbv7em-none-eabihf with --no-default-features)_

### License & compliance (final)
- [ ] Full `cargo deny check` clean on the release commit _(2026-11-26: `cargo deny check` passes on the working tree, with advisories, bans, licences and sources all ok. The release commit still does not exist, so this stays open)_
- [x] Regenerate `docs/dependency-audit.md` and a third-party notices file _(THIRD_PARTY.md generated from cargo metadata, 214 packages; docs/dependency-audit.md updated)_
- [x] Confirm zero-FFI: no `-sys` crates in the dependency tree (`cargo tree` check in CI) _(only windows-sys among -sys crates; pyo3-ffi only in the bindings, which CI now excludes on purpose; docs/zero-ffi.md)_
- [x] Confirm headers/Cargo metadata/LICENSE files consistent in all crates and in the Python package _(READMEs and readme keys are consistent; `LICENSE-MIT` and `LICENSE-APACHE` copied byte-identical into each library and CLI crate and into `bindings/tpt-bio-py`; crates inherit `license` from the workspace)_

### Release
- [x] Choose final crate names and check availability on crates.io _(decided: keep the current names. On 2026-11-25 the crates.io index returned 404 for tpt-bio, tpt-bio-core, tpt-bio-seq, tpt-bio-vision, tpt-bio-pipeline and tpt-bio-cli, and PyPI returned 404 for tpt-bio. Nothing is reserved, since nothing is published)_
- [ ] Publish order: core -> seq / vision -> pipeline -> cli
- [ ] `cargo publish --dry-run` for all crates _(run 2026-11-25, dry-run only, nothing uploaded: tpt-bio-core passes; seq, vision, pipeline and cli fail because tpt-bio-core is not on crates.io yet; tpt-bio-vision also depends on git pins, `tpt-crypto-hash` and `tpt-benchmark`, which crates.io rejects)_
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
