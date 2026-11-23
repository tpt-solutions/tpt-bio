# tpt-bio-pipeline

DAG-based streaming workflow engine for tpt-bio.

A pipeline is a directed acyclic graph of nodes. Each node takes the outputs of
its predecessors and produces a value. The crate checks for cycles when edges
are added, so every graph is acyclic by construction. Nodes can run sequentially
or level-by-level in parallel, can stream data through bounded queues, and can
checkpoint their outputs so an interrupted run resumes where it stopped.

## Concepts

| Item | Module | Purpose |
| --- | --- | --- |
| `Dag` | crate root | Nodes with payloads and edges; `execute` runs in dependency order, `execute_parallel` runs each independent level on the rayon pool |
| `Progress` | crate root | Callback notified when each node starts and finishes |
| Typed ports | `typed` | Port types so edges can be checked for compatible data |
| Checkpoints | `checkpoint` | Persist node outputs to disk so a rerun skips completed nodes |
| Streams | `stream` | Bounded, backpressured map over an iterator |
| Definitions | `definition` | Parse TOML pipeline files into a graph |
| Built-in nodes | `execute` | Statistics, sequence I/O and 2D vision steps |
| Sequence nodes | `seq_nodes` | Pipeline nodes backed by `tpt-bio-seq` |
| Vision nodes | `vision_nodes` | Pipeline nodes backed by `tpt-bio-vision` |

## Installation

```toml
[dependencies]
tpt-bio-pipeline = "0.1"
```

## Defining a pipeline in TOML

Pipelines can be described in a TOML file and loaded through `definition`. A
worked example lives at `examples/spec7_screen.toml`, and the `spec7_example`
integration test runs it end to end.

## Running a pipeline

Command-line access is through `tptbio pipeline run <definition>` in
`tpt-bio-cli`. The worked example is run by the test suite:

```sh
cargo test -p tpt-bio-pipeline
```

## Design notes

- `Dag::execute` runs nodes in dependency order. `Dag::execute_parallel` runs
  each level of independent nodes on the rayon thread pool.
- Edges are checked for cycles when they are added, so a built `Dag` is acyclic.
- `stream` provides a bounded map over an iterator, which limits how far a
  producer can run ahead of its consumer.

## Dependencies on the workspace

- `tpt-bio-core`, `tpt-bio-seq` and `tpt-bio-vision` supply the built-in nodes.
- `toml` and `serde` parse definitions.
- `crossbeam` and `rayon` provide channels and the parallel executor.

## Testing

```sh
cargo test -p tpt-bio-pipeline
```

Integration tests cover the built-in nodes (`tests/builtin_nodes.rs`) and the
worked example (`tests/spec7_example.rs`).

## Minimum supported Rust version

Rust 1.95, set in the workspace manifest.

## Licence

Dual-licensed under MIT OR Apache-2.0, at your option.
