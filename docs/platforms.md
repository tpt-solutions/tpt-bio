# Platform coverage

This page records which platforms the project has been built and tested on, and what
the CI matrix covers. Update it when either changes.

## Local testing

| Platform | Arch | Tested locally | Notes |
| --- | --- | --- | --- |
| Windows 11 Pro (10.0.26200) | x86_64 | Yes | MSVC toolchain. Build, tests, fuzzing (MSVC ASan) and benchmarks. |
| Linux | any | No | Covered by CI only. |
| macOS | any | No | Covered by CI only. |
| Windows on ARM, other BSDs, 32-bit targets | any | No | Not tested. |

Local testing is Windows x86_64 only. Nothing else has been run by hand.

## CI matrix (`.github/workflows/ci.yml`)

| Job | Platform | What it checks |
| --- | --- | --- |
| test | ubuntu-latest, macos-latest, windows-latest | `cargo test --workspace`, toolchain 1.95.0 |
| fmt | ubuntu-latest | `cargo fmt --check` (stable) |
| clippy | ubuntu-latest | `cargo clippy --workspace --all-targets -D warnings` (1.95.0) |
| no-std | ubuntu-latest, cross-target `thumbv7em-none-eabihf` | `tpt-bio-core` builds and tests with `--no-default-features` |
| msrv | ubuntu-latest | `cargo check --workspace --all-targets` on 1.95.0 |
| docs | ubuntu-latest | `cargo doc` with warnings as errors |
| deny | ubuntu-latest | `cargo deny check` (licences, advisories, sources) |
| python-bindings | ubuntu-latest | clippy on `bindings/tpt-bio-py` (separate workspace) |
| zero-ffi | ubuntu-latest | fails if a C-library `-sys` crate enters the dependency tree |

Other workflows:

- `.github/workflows/perf.yml`: criterion benches on ubuntu-latest, weekly and on manual
  dispatch. Informational only; it does not fail on regressions.
- `.github/workflows/wheels.yml`: Python wheels for linux, macOS and Windows, manual
  dispatch only, uploaded as workflow artifacts. Nothing is published.

## Gaps

- No CI job runs the fuzz targets or `proptest` with extended case counts.
- No aarch64 or big-endian target is tested. The BAM, BCF and BGZF code assumes
  little-endian reads through explicit `from_le_bytes`, but this has not been checked on
  a big-endian machine.
- The wheels workflow has not been run on a hosted runner yet.
