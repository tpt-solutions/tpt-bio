# Zero-FFI status

Checked 2026-11-24 against the workspace at commit `c11a2c4`.

tpt-bio is zero-FFI for its library and CLI crates. No C library is linked into
`tpt-bio-core`, `tpt-bio-seq`, `tpt-bio-vision`, `tpt-bio-pipeline` or `tptbio`.
The one expected exception is the Python extension, `bindings/tpt-bio-py`, which
links the Python C API through `pyo3-ffi`. That link is required for any PyO3
extension module and is not a third-party C library dependency.

## How the check works

The `zero-ffi` job in `.github/workflows/ci.yml` runs:

```sh
cargo tree --workspace --exclude tpt-bio-py --edges normal,build --prefix none \
  | grep -E '^[a-zA-Z0-9_]+(-[a-zA-Z0-9_]+)*-(sys|ffi) ' \
  | grep -v '^windows-sys '
```

It fails if any package whose name ends in `-sys` or `-ffi` appears in the normal or
build dependency graph. Two exclusions are deliberate:

| Excluded | Why |
| --- | --- |
| `tpt-bio-py` (`--exclude`) | The bindings crate pulls in `pyo3-ffi` (links `python`). That is the expected FFI link for a Python extension. Every other workspace crate is still checked. |
| `windows-sys` (`grep -v`) | Microsoft's raw Windows API binding, which std-adjacent crates use on Windows. It wraps no third-party C library. |

The `-ffi` suffix was added to the pattern. The earlier pattern (`-sys` only) would
not have caught `pyo3-ffi`, so it did not express the bindings exception at all.

## Result

Measured with `cargo tree --workspace --exclude tpt-bio-py --edges normal,build --prefix none --target all --offline`:

- No package ending in `-sys` or `-ffi` other than `windows-sys`.
- No `cc`, `cmake`, `pkg-config` or `vcpkg` crate, so no crate in the tree
  compiles C code through a build script.
- No package in the tree declares a `links` key naming a C library.

Measured with `cargo tree --workspace --edges normal,build --prefix none --target all --offline`
(including the bindings crate):

- `pyo3-ffi` 0.30.1 (`links = "python"`), reached through `pyo3` and `numpy` from
  `tpt-bio-py`. This is the expected exception.
- `windows-sys` 0.61.2, excluded as above.

Packages that are not `-sys` crates but are worth knowing about:

| Package | Why it is present | Note |
| --- | --- | --- |
| `libc` 0.2.191 | Used by `rayon-core` and other std-adjacent crates on Unix | Raw bindings to the platform C library supplied by the OS. It is not a `-sys` crate and the policy check does not match it. Its build script does not compile C. |
| `rayon-core` 1.13.0 | Parallel executor | Declares `links = "rayon-core"`. That key only prevents duplicate builds, and nothing is linked. |
| `wasm-bindgen-shared` 0.2.130 | Appears in the lockfile metadata, through `wasm-bindgen` (wasm targets) | Declares `links = "wasm_bindgen"`. Not C. `cargo tree --target all` does not list it in the workspace's normal and build graph. |

`linux-raw-sys` appears in the full metadata because `rustix` depends on it on Linux.
It did not appear in the Linux or macOS trees for these crates when checked with
`--target`, so it is not part of the CI dependency set for the library crates.
If a future change makes it appear, it is pure Rust (generated syscall bindings)
and would need an explicit exclusion in the CI grep with a note.

## Re-running the check

```sh
cargo tree --workspace --exclude tpt-bio-py --edges normal,build --prefix none --target all
cargo tree -i pyo3-ffi --workspace --edges normal,build
```

The first command must print only `windows-sys` lines for `-sys`/`-ffi` names. The
second shows the bindings path to `pyo3-ffi`.
