# Dependency license audit

Recorded 2026-11-24. `cargo deny check` enforces the policy in `deny.toml`. Re-run this
audit when adding a dependency or before a release.

Policy: MIT, BSD-2/3, ISC, Zlib, Unicode-3.0, or a dual license that includes MIT.
Apache-2.0-only, GPL/LGPL/AGPL, and MPL are rejected. Git sources are allowed only
for the two repositories listed under `allow-git` in `deny.toml`, and each is pinned
to a commit.

## Toolchain

| Item | Value | Note |
| --- | --- | --- |
| MSRV | Rust 1.95 | Raised from 1.85 for `tpt-simd` 0.1.1 and `target-lexicon-mit` 0.13.6. |

## Substrate crates (tpt-solutions)

| Dependency | Source | License | Use | Status |
| --- | --- | --- | --- | --- |
| `tpt-crypto-hash` | git, `tpt-solutions/tpt-crypto` @ `e08fafb` | MIT OR Apache-2.0 | SHA-256 and HMAC-SHA256 for the SigV4 signer (`tpt-bio-vision`) | In use. Pre-release and unaudited. The `tpt-crypto` workspace says so in its README. Not on crates.io, so it cannot be published to crates.io while pinned to git. |
| `target-lexicon` 0.13.6 (patched) | git, `tpt-solutions/target-lexicon-mit` @ `297f442` | MIT | Replaces the upstream crate through `[patch.crates-io]` in the root `Cargo.toml`. Needed by `pyo3-build-config`. | In use. This is a clean-room rewrite of the subset PyO3 uses, not a relicensed copy of upstream. Upstream is `Apache-2.0 WITH LLVM-exception`, and that license is no longer in the tree. |
| `tpt-simd` 0.1.1 | crates.io | MIT OR Apache-2.0 | SIMD 2D separable convolution (`tpt-bio-vision::filter`) | In use. Needs Rust 1.95. |
| `tpt-strata` 0.1.0 | crates.io | MIT OR Apache-2.0 | In-process columnar table and CSV export (`tpt-bio-seq::table`). This is the Arrow replacement. | In use. The core crate has zero dependencies. |
| `tpt-strata-parquet` 0.1.0 | crates.io | MIT OR Apache-2.0 | Parquet bridge | **Not used.** Parquet is out of scope for this project, and the bridge depends on `arrow` 60 and `parquet` 60 (Apache-2.0 only). |
| `tpt-math-prob-dist` 0.2.0 | crates.io | MIT OR Apache-2.0 | Poisson and Binomial functions (`tpt-bio-core`, `std` feature only) | In use. See the Poisson decision in `todo.md`. |
| `tpt-math` 0.2.0 | crates.io | MIT OR Apache-2.0 | Not a direct dependency. Its `prob` feature enables `std` throughout, which breaks the `no_std` build. | Not used directly. |

## Other planned and in-use dependencies

| Dependency | Version checked | License | Status |
| --- | --- | --- | --- |
| `flate2` | 1.1.11 | MIT OR Apache-2.0 | In use. `rust_backend` (miniz_oxide, Zlib-licensed). |
| `memmap2` | 0.9.12 | MIT OR Apache-2.0 | In use (`tpt-bio-seq`). |
| `rayon` | 1.12.0 | MIT OR Apache-2.0 | In use. |
| `crossbeam` | 0.8.x | MIT OR Apache-2.0 | In use (`tpt-bio-seq` batch channel). |
| `libm` | 0.2.16 | MIT OR Apache-2.0 | In use (`tpt-bio-core`). |
| `clap` | 4.6.9 | MIT OR Apache-2.0 | In use (`tpt-bio-cli`). |
| `pyo3` | 0.30.1 | MIT OR Apache-2.0 | In use (`bindings/tpt-bio-py`). |
| `image` | 0.25.9 | MIT OR Apache-2.0 | Planned, Phase 3. |
| `tiff` | 0.11.3 | MIT | Planned, Phase 3 (OME-TIFF). |
| `ndarray` | 0.17.2 | MIT OR Apache-2.0 | Planned. |
| `numpy` (rust-numpy) | 0.30.0 | BSD-2-Clause | Planned (Python interop). |
| `wgpu` | 26.0.1 | MIT OR Apache-2.0 | Planned, Phase 4. The latest is 30.0.2; check its MSRV before upgrading. |
| `ureq` | 3.5.2 | MIT OR Apache-2.0 | Planned, HTTP range reads. Its TLS stack still needs a decision (see below). |
| `serde`, `toml` | 1.0.229, 1.2.2 | MIT OR Apache-2.0 | Planned, pipeline definitions. |
| `criterion` | 0.7.0 | Apache-2.0 OR MIT | Planned, dev-dependency for benchmarks. |
| `tokio` | 1.54.1 | MIT | Not expected to be needed. |
| `rustls` | 0.24.0-dev | Apache-2.0 OR ISC OR MIT | Not used. Its default crypto providers need review. |

## Dev-only tooling

| Tool | Version | License | Status |
| --- | --- | --- | --- |
| `libfuzzer-sys` | 0.4.13 | `(MIT OR Apache-2.0) AND NCSA` | **Accepted as a dev/tooling dependency.** It lives only in the separate `fuzz/` workspace, which `cargo deny` does not check. NCSA is not on the allow list, and that is intentional. |
| `cargo-deny` | 0.20.2 | MIT OR Apache-2.0 | Runs the policy check. |
| `cargo-fuzz` | 0.13.2 | MIT OR Apache-2.0 | Runs the fuzz targets. |
| `samtools`, `bcftools` | not installed | MIT | Oracles for conformance tests and benchmarks. Never linked. |

## Rejected

| Dependency | Reason | Replacement |
| --- | --- | --- |
| `arrow`, `parquet` (and `tpt-strata-parquet`, which depends on them) | Apache-2.0 only, and Parquet is out of scope | `tpt-strata` for in-process tables. Tabular export is CSV (`tptbio seq export --format csv`) and JSON Lines (`--format jsonl`). |
| `target-lexicon` (upstream) | `Apache-2.0 WITH LLVM-exception` | `target-lexicon-mit` via `[patch.crates-io]`. |
| `ring`, any `*-sys` crate | Zero-FFI policy | `tpt-crypto`. `windows-sys` is excluded from the check because it is Microsoft's OS binding, not a third-party C library. |

## Open questions

- **TLS for remote I/O.** `tpt-crypto` gives SHA-256 and HMAC, but not a TLS client.
  The companion `tpt-crypto-rustls` repository is unaudited. Phase 3 HTTPS needs a decision.
- **Crates.io publishing.** The `tpt-crypto` and `target-lexicon-mit` pins are git
  dependencies. crates.io does not accept git dependencies, so Phase 6 needs a published
  version of each, or a vendored copy.

## Current `cargo deny` result

`cargo deny check` passes on advisories, bans, licenses, and sources. The licenses in the
tree are `MIT`, `Apache-2.0` (always as a member of an `OR` expression that includes MIT),
`BSD-2-Clause`, `Unicode-3.0`, `Zlib`, and `0BSD` (also an `OR` member). The string
`LLVM-exception` no longer appears in the tree.
