# Test data

The test suite is self-contained. It generates its fixtures in code, so `cargo test`
needs no downloads.

- `tpt-bio-core`: unit tests use literal sequences and known constants.
- `tpt-bio-seq`: BGZF, FASTQ, SAM and BAM tests build their inputs in memory, and
  `tests/bam_index.rs` generates a 2000-read BAM with a deterministic PRNG.

## External fixtures (not yet chosen)

Larger real-world fixtures are needed for conformance tests against samtools and
for benchmarks. Each one must have its license recorded here before it is added.

| Data | Needed for | Source | License | Status |
| --- | --- | --- | --- | --- |
| Small BAM with index | Phase 1 conformance | TBD | TBD | Not chosen |
| Human-genome-scale BAM | Phase 1 milestone benchmark | TBD | TBD | Not chosen |
| VCF, BCF, CRAM samples | Phase 2 conformance | TBD | TBD | Not chosen |
| OME-TIFF, Zarr plates | Phase 3 | TBD | TBD | Not chosen |
| BBBC reference images | Segmentation validation | TBD | TBD | Not chosen |

Any external data goes in `data/`, which `.gitignore` excludes from git.
