# Changelog

All notable changes to `tpt-bio-cli` are recorded here. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and the crate uses
[Semantic Versioning](https://semver.org/).

## [Unreleased]

### Added

- `tptbio` binary with `seq`, `vision`, `pipeline`, `stats` and `view` command groups.
- `seq view`, `seq index`, `seq flagstat`, `seq sort` and `seq export`.
- `vision segment` and `vision features`.
- `pipeline run` for TOML-defined pipelines.
- `stats normalize` and `stats pca`.
- `completions` and `man` subcommands.
- `--config` file support and a global `--output-format` option.
- Integration tests for the CLI and the vision pipeline.

## [0.1.0]

Initial workspace version. Not yet published.
