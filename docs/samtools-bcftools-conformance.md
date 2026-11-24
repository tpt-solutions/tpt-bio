# Conformance with samtools, bcftools and tabix

Checked against samtools 1.22.1, bcftools 1.22 and tabix 1.22.1, run in WSL Ubuntu.
Nothing here links against those tools. The committed tests use the fixtures under
`tpt-bio-seq/tests/data/` and run without them. CRAM is covered in
[cram-conformance.md](cram-conformance.md).

## Results

| Format | Direction | Result |
|---|---|---|
| BAM | samtools reads our BAM | Identical records, `quickcheck` passes |
| BAM | we read samtools' BAM | Identical records |
| BAI | our index from a samtools BAM | Byte-identical for the 250-record fixture. For 20,000 and 300,000 reads, the bins, chunks, linear index and metadata are identical once parsed; only the order of the bins in the file differs (see below). |
| BAI | region queries on a 300,000-read BAM | Six regions return the same reads as `samtools view`. |
| BAI | region query `chr1:1000-3000` | Same read names as `samtools view` |
| BCF | we read bcftools' BCF | Identical records (240 generated) |
| BCF | bcftools reads our BCF | Identical records to its own file |
| TBI | our index for a bgzipped VCF | Identical to `tabix` once decompressed |
| CSI | our index for a BCF | Identical to `bcftools index` once decompressed |
| VCF.gz and BCF regions | seven regions, including empty and contig edges | Same records as `tabix` and `bcftools view` |

## Bugs this found and fixed

- **BCF empty ID.** An empty string was written as type 0 (null). bcftools rejects it.
  An empty ID is now a zero-length `char`.
- **TBI and CSI trailer.** Both formats end with an 8-byte count of records without a
  coordinate. We omitted it.
- **CSI depth 0.** bcftools writes depth 0 when every contig fits in one bin. We rejected
  it, so the reader could not open bcftools' own CSI files.
- **BAI pseudo-bin.** samtools stores, for each reference, a pseudo-bin 37450 holding
  the record span and the mapped and unmapped counts. We left it out, so the index bytes
  differed.
- **BAI bins at scale.** htslib groups consecutive records in one bin into a single
  chunk, then moves small bins into their parent and merges chunks that share a BGZF
  block. We merged overlapping chunks as records arrived, so the bins differed on large
  files. Our builder now follows htslib's rules.
- **BAI linear index.** htslib fills an empty window from the next window with an
  offset, not from the previous one. We filled forward, which gave earlier offsets in
  about 10% of windows at 20,000 reads. The builder (and the TBI builder) now fills
  backwards.
- **BAI file order.** samtools writes bins in hash-table order, we write them in
  ascending order. Readers do not depend on this. The byte comparison in the tests
  therefore only holds where the two orders agree.

## Not checked

- **BCF types.** The generated data covers integers, floats, strings, flags, missing
  values, multi-allelic sites and genotypes. It does not cover every integer width, the
  `Number=R` and `Number=G` cases in every position, or large contigs (which give a
  deeper CSI).
- **Large BAM.** The 300,000-read timing file is synthetic and not committed. Deep
  indexes and the CSI format at depth above 0 are not compared on large data.
- **Compressed-file bytes.** Indexes are compared after decompression. The BGZF block
  boundaries of our output are not the same as tabix's.

## Reproducing

The fixtures were made with samtools, bcftools and tabix in the WSL distribution:

```sh
samtools view -b sorted.sam -o sorted.bam && samtools index sorted.bam
bcftools view -O b -o oracle.bcf in.vcf && bcftools index oracle.bcf
bgzip -c in.vcf > oracle.vcf.gz && tabix -p vcf oracle.vcf.gz
```

The examples in `tpt-bio-seq/examples/` run our side of each comparison:
`bai_oracle`, `bam_oracle_write`, `bcf_oracle_dump`, `bcf_oracle_write` and
`index_oracle`.

## Timings

Median of five runs, in milliseconds, on this machine, for a synthetic BAM of 300,000
100-base reads over three contigs (11 MB). samtools ran in WSL Ubuntu and tptbio ran
natively on Windows (release build). The two environments are not identical, and each
time includes process start-up. Treat these as indicative only.

| Operation | samtools 1.22.1 | tptbio |
|---|---|---|
| `view` to SAM text (all records) | 341 | 521 |
| `view` region `chr2:1000000-3000000` (indexed) | 30 | 48 |
| index build | 290 | 176 |

Our records and region results were checked against samtools on the same file. Human
genome data was not available, so the milestone for a human-genome BAM remains open.

