# Changelog

All notable changes to `tpt-bio-core` are recorded here. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and the crate uses
[Semantic Versioning](https://semver.org/).

## [Unreleased]

### Added

- `stats` module with discrete distributions, k-means and agglomerative clustering,
  PCA, t-SNE and UMAP.
- Column and group normalisation, phenotype aggregation, and plate QC metrics
  (Z' factor and SSMD).
- Seeded random number generator for reproducible clustering and embeddings.
- Property tests under `tests/properties.rs`.

### Changed

- Poisson distribution functions now run without the `std` feature, using `libm`.
  Only the Binomial distribution still requires `std`.

## [0.1.0]

Initial workspace version. Not yet published.
