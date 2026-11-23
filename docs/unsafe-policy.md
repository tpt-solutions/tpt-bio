# Unsafe code policy

The workspace sets `unsafe_code = "deny"` (see the root `Cargo.toml`). Crates that
can be written without `unsafe` also say so with `#![forbid(unsafe_code)]`.

## Inventory

| Location | Why | Contract |
| --- | --- | --- |
| `tpt-bio-seq/src/mmap.rs`, `MappedFile::open` | `memmap2::Mmap::map` is `unsafe` because the OS does not stop another process from changing a mapped file. | The file must not be truncated or modified while the map is alive. The API documents this. The mapping is read-only. |

`tpt-bio-core` is `#![forbid(unsafe_code)]`. The rest of the workspace has no
`unsafe` code.

## Rules for new unsafe code

1. Put it in its own module, and add `#[allow(unsafe_code)]` to that item only.
2. Add a `// SAFETY:` comment to every `unsafe` block, stating the invariant it relies on.
3. Add it to this inventory.
4. Add tests that exercise the boundary. Run them under Miri where the platform
   allows it.

## Planned SIMD exceptions

`tpt-simd` (Phase 3) is expected to provide its own safe interfaces, so tpt-bio
should need no `unsafe` for SIMD. If that changes, the exception will be recorded here.
