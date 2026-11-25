# Test data

The test suite is self-contained. It generates its fixtures in code, so `cargo test`
needs no downloads.

- `tpt-bio-core`: unit tests use literal sequences and known constants.
- `tpt-bio-seq`: BGZF, FASTQ, SAM and BAM tests build their inputs in memory, and
  `tests/bam_index.rs` generates a 2000-read BAM with a deterministic PRNG.

## Decisions

No third-party dataset is downloaded for tests. Each fixture is one of:

1. generated in code by the test itself;
2. written by a permissively licensed reference tool (samtools and bcftools 1.22, MIT
   and BSD; ome-zarr 0.21.0, BSD-3-Clause), run once, and committed only when small. Its
   provenance is recorded next to it, for example in `tpt-bio-vision/tests/fixtures/hcs/README.md`.

| Data | Needed for | Decision | Status |
| --- | --- | --- | --- |
| Small BAM with index | Phase 1 conformance | Generated in code, and checked against samtools 1.22 output run in WSL | Chosen: generated |
| Human-genome-scale BAM | Phase 1 milestone benchmark | No public dataset is committed. The benchmark uses a generated 300,000-read BAM, and the milestone stays open until real data is measured | Open |
| VCF, BCF, CRAM samples | Phase 2 conformance | Generated in code, and checked against bcftools and samtools | Chosen: generated |
| OME-TIFF, Zarr plates | Phase 3 | OME-TIFF from the crate's own encoder. HCS plate metadata from ome-zarr 0.21.0 | Chosen: generated |
| BBBC reference images | Segmentation validation | The BBBC host did not answer from this environment, so it can't be licence-checked. Validation uses synthetic images with exact masks, through the IoU and F1 matching code | Open: BBBC when reachable |

Any external data goes in `data/`, which `.gitignore` excludes from git. Nothing is in
`data/` yet.
