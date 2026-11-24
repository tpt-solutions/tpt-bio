# tpt-bio-core

Shared biological types, statistics and bio-math for the tpt-bio workspace.

`tpt-bio-core` is the leaf crate of the workspace. It holds the value types and
numerical routines that the sequencing, imaging and pipeline crates build on. It
has no I/O and no `unsafe` code, and it is `no_std` + `alloc` when the default
`std` feature is turned off, so it can run on embedded targets.

## Contents

| Module | What it provides |
| --- | --- |
| `nucleotide` | DNA/RNA base types, complement and reverse complement, validation |
| `amino` | Amino-acid types and codon translation |
| `coords` | Genomic coordinates and half-open intervals, with 0-based and 1-based conversion |
| `kmer` | Packed k-mer encoding, rolling k-mer iteration, canonical k-mers |
| `quality` | Phred and Phred+33 quality score conversion and error probabilities |
| `stats` | Discrete distributions (Poisson, Binomial), k-means and agglomerative clustering, PCA, t-SNE, UMAP, normalisation, phenotype aggregation and plate QC (Z', SSMD) |
| `par` | Rayon-based parallel helpers (requires `std`) |
| `error` | The shared error type used across modules |

## Installation

```toml
[dependencies]
tpt-bio-core = "0.1"
```

For `no_std` builds, disable default features:

```toml
[dependencies]
tpt-bio-core = { version = "0.1", default-features = false }
```

## Features

| Feature | Default | Effect |
| --- | --- | --- |
| `std` | yes | Enables rayon-based parallel helpers and the Binomial distribution, which depends on `tpt-math-prob-dist`. Poisson works without it. |

## Example

```rust
use tpt_bio_core::nucleotide::reverse_complement;

assert_eq!(reverse_complement(b"ACGTN").unwrap(), b"NACGT".to_vec());
```

## Behaviour notes

- Invalid distribution parameters return `NaN` rather than panicking.
- k-means and t-SNE take an explicit seed, so their results are reproducible for
  a given input and seed.
- The t-SNE and UMAP implementations have not been compared against reference
  implementations such as scikit-learn; treat their output as unvalidated.
- Poisson distribution functions are computed with `libm` (log-gamma and the
  regularised incomplete gamma function).

## Testing

```sh
cargo test -p tpt-bio-core
cargo test -p tpt-bio-core --no-default-features   # no_std build
```

The integration tests in `tests/properties.rs` use `proptest` to check invariants
on the core types.

## Minimum supported Rust version

Rust 1.95, set in the workspace manifest.

## Status

Version 0.1.0, pre-release. The modules in the Contents table are the implemented surface described in the crate documentation. Open items are tracked in `todo.md` in the repository root.

## Licence

Dual-licensed under MIT OR Apache-2.0, at your option.
