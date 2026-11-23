# tpt-bio-cli

The `tptbio` command-line tool: sequencing, imaging, statistics and pipeline
commands from one binary.

## Installation

From the workspace root:

```sh
cargo install --path tpt-bio-cli
```

This installs a binary named `tptbio`.

## Command overview

```text
tptbio seq          sequencing files: view, index, flagstat, sort, export
tptbio vision       image work: segment, features
tptbio pipeline     run a TOML-defined pipeline: run
tptbio stats        statistics: normalize, pca
tptbio view         inspect an image file
tptbio completions  print a shell completion script
tptbio man          write man pages to a directory
```

Run `tptbio <command> --help` for the full options of each command.

## Global options

| Option | Effect |
| --- | --- |
| `--config <file>` | Load defaults from a TOML file |
| `--output-format <fmt>` | Output format for commands that print tables; command-line values override the config file |

The config file's `threads` key is validated (it must be at least 1) and
recorded, but the CLI does not currently create a thread pool from it.

## Sequencing commands

| Command | Purpose |
| --- | --- |
| `tptbio seq view` | Print records from SAM, BAM, VCF, bgzipped VCF (`.vcf.gz`) or BCF. The format is detected from the file. Region queries need a BAM input and its index. |
| `tptbio seq index` | Build an index for a sorted BAM |
| `tptbio seq flagstat` | Summarise alignment flags |
| `tptbio seq sort` | Sort alignments by coordinate |
| `tptbio seq export` | Convert a file to another supported format |

FASTQ and FASTA are read by the library crates but are not accepted by `seq view`.

## Vision commands

| Command | Purpose |
| --- | --- |
| `tptbio vision segment` | Threshold an image with the chosen method (for example Otsu) and write the result |
| `tptbio vision features` | Measure per-object shape and intensity features and write them as a table |

## Pipeline command

| Command | Purpose |
| --- | --- |
| `tptbio pipeline run <definition>` | Load a TOML pipeline definition and execute it |

## Statistics commands

| Command | Purpose |
| --- | --- |
| `tptbio stats normalize` | Normalise columns. The method is z-score or robust (median and MAD based). |
| `tptbio stats pca` | Principal component analysis |

## Image inspection

`tptbio view <input>` opens an image file and reports its contents.

## Completions and man pages

```sh
tptbio completions bash > tptbio.bash
tptbio man ./man            # writes one page per command into ./man
```

## Exit codes

- `0`: success.
- Non-zero: a usage error, a sequencing error or a runtime error. Each is reported
  on standard error with a short message.

## Source layout

| Module | Role |
| --- | --- |
| `main.rs` | Argument parsing and command dispatch |
| `config.rs` | Loading and validating the TOML config file |
| `output.rs` | Table and CSV formatting |
| `variants.rs` | Format detection and viewing for VCF and BCF |

## Testing

```sh
cargo test -p tpt-bio-cli
```

`tests/cli.rs` exercises the binary's commands, and `tests/vision_pipeline.rs`
runs an end-to-end vision pipeline.

## Minimum supported Rust version

Rust 1.95, set in the workspace manifest.

## Licence

Dual-licensed under MIT OR Apache-2.0, at your option.
