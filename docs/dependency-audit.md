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
| `wgpu` | 26.0.1 | MIT OR Apache-2.0 | **Rejected.** Its tree includes `codespan-reporting` (Apache-2.0 only, a non-optional dependency of naga), `gl_generator`, `glutin_wgl_sys` (Apache-2.0 only), `hexf-parse` (CC0-1.0) and `paste` (RUSTSEC-2024-0436). No feature selection removes them. Browser WebGPU will use `web-sys` instead. |
| `ash` | 0.38.0+1.3.281 | MIT OR Apache-2.0 | **Passes** (audited 2026-11-25, see "Native GPU candidates"). Not yet used. Candidate Vulkan backend for `render3d`, fed offline-compiled SPIR-V. |
| `libloading` | 0.8.9 | ISC | **Passes.** Only dependency of `ash` besides `windows-link`. It loads `vulkan-1.dll` at run time, so it is FFI in the dlopen sense (see the zero-FFI note below). **Decided 2026-11-26: accepted.** It is not a `-sys` crate, so the zero-FFI CI rule holds. |
| `wry` | 0.57.3 | Apache-2.0 OR MIT | **Fails the policy on Windows** (audited 2026-11-25). `dpi` and `dunce` are rejected. Not used. See "Native GPU candidates". |
| `ureq` | 3.5.2 | MIT OR Apache-2.0 | Planned, HTTP range reads. Its TLS stack still needs a decision (see below). |
| `serde`, `toml` | 1.0.229, 1.2.2 | MIT OR Apache-2.0 | Planned, pipeline definitions. |
| `tpt-benchmark` (git) | 0.1.0 | MIT OR Apache-2.0 | **Adopted**, dev-dependency for benchmarks (`tpt-bio-seq`). Replaced `criterion`. Pinned to `8384f44` of `tpt-solutions/tpt-benchmark`, which is not yet on crates.io. |
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

## Native GPU candidates (2026-11-25)

wgpu was rejected (above), so native GPU rendering has no backend. Two replacements were
audited in scratch projects outside the repo, with this repository's `deny.toml` (`all-features`
off) and `cargo deny --target x86_64-pc-windows-msvc check licenses bans advisories`.

| Candidate | Packages (Windows) | licenses | bans | advisories |
| --- | --- | --- | --- | --- |
| `ash` 0.38.0 | 3: `ash`, `libloading` 0.8.9 (ISC), `windows-link` 0.2.1 | ok | ok | ok |
| `wry` 0.57.3 (webview shell) | 38 | **FAILED** | ok | ok |

`wry` failures, both transitive:

- `dpi` 0.1.2 is `Apache-2.0 AND MIT`. The `AND` needs Apache-2.0, which the policy rejects.
  This one is a real conflict and no `deny.toml` setting fixes it short of allowing Apache-2.0.
- `dunce` 1.0.5 is `CC0-1.0 OR MIT-0 OR Apache-2.0`. MIT-0 is permissive and MIT-compatible but
  is not on the allow list. Adding `MIT-0` to `[licenses] allow` would clear it.

Other points for `wry`:

- `webview2-com-sys` is a `*-sys` crate. The zero-FFI rule excepts only `windows-sys`, so this
  needs an explicit exception or a decision.
- Only the Windows target was checked. macOS and Linux pull in other webview stacks (WebKit,
  GTK) that were not audited.

Points for `ash`:

- Every package passes with no policy change.
- Vulkan is loaded with `libloading` at run time, so this is FFI. It is not a `*-sys` crate and
  no `deny.toml` ban matches it, but it is a departure from the zero-FFI intent. Needs an
  explicit decision before adoption.
- `ash` takes SPIR-V only. WGSL would be compiled offline (outside the dependency tree) and the
  `.spv` committed. naga is not needed.

Nothing here has been added to a `Cargo.toml`.

## Open questions

- **TLS for remote I/O.** `tpt-crypto` gives SHA-256 and HMAC, but not a TLS client.
  The companion `tpt-crypto-rustls` repository is unaudited. Phase 3 HTTPS needs a decision.
- **Crates.io publishing.** The `tpt-crypto` and `target-lexicon-mit` pins are git
  dependencies. crates.io does not accept git dependencies, so Phase 6 needs a published
  version of each, or a vendored copy.

## Current workspace dependency list (2026-11-24)

Direct dependencies from the workspace manifests, with the locked version and the licence
the crate declares. The full transitive list (214 third-party packages) and the per-package
licence check are in `THIRD_PARTY.md`, generated by `scripts/gen_third_party.py`.

| Crate | Direct dependencies (locked version, licence) |
| --- | --- |
| `tpt-bio-core` | `libm` 0.2.16 (MIT); `rayon` 1.12.0 (MIT OR Apache-2.0), optional; `tpt-math-prob-dist` 0.2.0 (MIT OR Apache-2.0), optional, `std` only |
| `tpt-bio-seq` | `tpt-bio-core`; `flate2` 1.1.11 (MIT OR Apache-2.0); `rayon` 1.12.0; `crossbeam` 0.8.5 (MIT OR Apache-2.0); `tpt-strata` 0.1.0 (MIT OR Apache-2.0); `memmap2` 0.9.12 (MIT OR Apache-2.0), optional, `mmap` feature |
| `tpt-bio-vision` | `tpt-bio-core`, `tpt-bio-seq`; `tpt-crypto-hash` 0.1.0 (git, MIT OR Apache-2.0); `tpt-simd` 0.1.1 (MIT OR Apache-2.0); `image` 0.25.10 (MIT OR Apache-2.0); `tiff` 0.9.1 (MIT); `ndarray` 0.16.1 (MIT OR Apache-2.0); `quick-xml` 0.41.0 (MIT); `serde` 1.0.229, `serde_json` 1.0.154 (MIT OR Apache-2.0); `flate2` 1.1.11 |
| `tpt-bio-pipeline` | `tpt-bio-core`, `tpt-bio-seq`, `tpt-bio-vision`; `crossbeam` 0.8.5; `rayon` 1.12.0; `serde` 1.0.229; `toml` 0.8.23 (MIT OR Apache-2.0) |
| `tpt-bio-cli` | `clap` 4.6.9, `clap_complete` 4.6.12, `clap_mangen` 0.2.33 (all MIT OR Apache-2.0); `serde`, `serde_json`, `toml`; the four library crates |
| `bindings/tpt-bio-py` | `pyo3` 0.30.1 (MIT OR Apache-2.0); `numpy` 0.30.0 (BSD-2-Clause); the four library crates |
| dev-only | `proptest` 1.12.0 (MIT OR Apache-2.0); `tpt-benchmark` 0.1.0 (git, MIT OR Apache-2.0, benches in `tpt-bio-seq`) |

Patched and git-sourced packages:

| Package | Source | Licence | Effect |
| --- | --- | --- | --- |
| `tpt-crypto-hash` 0.1.0 (and `tpt-crypto-core`, `tpt-crypto-ct`) | git, `tpt-solutions/tpt-crypto` @ `e08fafb` | MIT OR Apache-2.0 | **Blocks crates.io publishing.** crates.io rejects git dependencies. |
| `target-lexicon` 0.13.6 | git, `tpt-solutions/target-lexicon-mit` @ `297f442`, via `[patch.crates-io]` | MIT | **Blocks crates.io publishing.** The patch is not applied on crates.io. |

Both pins are commits, not released versions. Phase 6 needs a published release of each, or
a vendored copy, before any crate that depends on them can be published.

## Licence findings (2026-11-24)

- No third-party package lacks a licence field. No package uses MPL-2.0, GPL, LGPL-only, or
  an unlisted licence. `LGPL-2.1-or-later` appears only as an alternative in `r-efi`, which is
  `MIT OR Apache-2.0 OR LGPL-2.1-or-later`, so MIT is available and is the one used.
- Apache-2.0-only packages: `ciborium` 0.2.2, `ciborium-io` 0.2.2 and `ciborium-ll` 0.2.2
  came in through `criterion`, which has been removed, so they are no longer in the lock.
  `borsh-derive` 1.8.1 is still in `Cargo.lock` but not in the active dependency graph
  (`cargo tree` finds no path to it), so `cargo deny` passes. The generator flags it anyway.
  Remove the stale lock entry (`cargo update -p borsh` or a clean re-lock) to clear it.
- `zlib-rs` 0.6.8 is Zlib-only and is on the allow list.
- `libfuzzer-sys` is not in the main workspace lock. It lives only in `fuzz/`.

## Current `cargo deny` result (2026-11-24)

`cargo deny check` (all categories, run in the workspace root):

- **bans**: ok. Warnings only: duplicate versions of `miniz_oxide`, `syn` and `tiff`
  (`tiff` 0.9.1 and 0.11.3).
- **licenses**: ok.
- **sources**: ok.
- **advisories**: ok after the upgrade. The check first failed on two advisories against
  `quick-xml` 0.37.5 (a dependency of `tpt-bio-vision`): a quadratic run time when checking a
  start tag for duplicate attribute names, and unbounded namespace-declaration allocation in
  `NsReader`. The working copy now pins `quick-xml = "0.41"` (lock: 0.41.0), which fixes both.
  That change was uncommitted when this audit ran, and it was made outside this audit, so it
  should be committed with a fresh `cargo deny check` before release.

The earlier statement that the advisory check passed is out of date. The earlier run failed
on those advisories, and this run passes.

The licences in the tree are `MIT`, `Apache-2.0` (always as a member of an `OR` expression
that includes MIT, except for the four packages listed above), `BSD-2-Clause`, `BSD-3-Clause`,
`Unicode-3.0`, `Zlib`, `0BSD` (also an `OR` member), `Unlicense` (also an `OR` member), and
`(MIT OR Apache-2.0) AND Unicode-3.0`. The string `LLVM-exception` no longer appears in the
tree, except in `linux-raw-sys` and `rustix`, where `Apache-2.0 WITH LLVM-exception` is an
`OR` alternative and MIT is also available.
