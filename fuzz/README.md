# Fuzzing

This directory is a separate Cargo workspace (see `fuzz/Cargo.toml`), so
`libfuzzer-sys` does not enter the main dependency graph.

Install the tool once, then run a target:

```sh
cargo install cargo-fuzz
cd fuzz
cargo fuzz run bgzf_decompress
cargo fuzz run fastq_parse
cargo fuzz run sam_parse
cargo fuzz run bam_read
cargo fuzz run bai_read
```

Crashes are written to `fuzz/artifacts/`, and the corpus goes to `fuzz/corpus/`.
Both are git-ignored. A crash found by fuzzing becomes a regression test in the
crate it affects, with the input stored inline.

`libfuzzer-sys` is licensed `(MIT OR Apache-2.0) AND NCSA`. See
`docs/dependency-audit.md` for why it sits outside the `cargo deny` policy.
