# Changelog

All notable changes to `tpt-bio-pipeline` are recorded here. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and the crate uses
[Semantic Versioning](https://semver.org/).

## [Unreleased]

### Added

- `Dag` model with cycle checking on edge insertion.
- Sequential and level-parallel execution (`execute`, `execute_parallel`).
- `Progress` callbacks for node start and finish.
- Typed ports (`typed` module).
- Checkpointing so an interrupted run can resume (`checkpoint` module).
- Bounded, backpressured streaming map (`stream` module).
- TOML pipeline definitions (`definition` module).
- Built-in sequence and vision nodes (`seq_nodes`, `vision_nodes`).
- Worked screening example (`examples/spec7_screen.toml`) and integration tests.

## [0.1.0]

Initial workspace version. Not yet published.
