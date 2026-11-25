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
time includes process start-up. samtools' own medians ranged from 269 to 461 ms for
`view` across runs, so treat these as indicative only.

| Operation | samtools 1.22.1 | tptbio |
|---|---|---|
| `view` to SAM text (all records) | 298 | 304 |
| `view` region `chr2:1000000-3000000` (indexed) | 14 | 51 |
| index build | 273 | 94 |

Changes made to reach these numbers: SAM lines are written into one reused buffer,
BAM sequences decode through a lookup table, and BGZF blocks are inflated in parallel
in batches (read-ahead starts at 64 KiB after each seek, and grows to 4 MiB).

Our records and region results were checked against samtools on the same file. Human
genome data was not available, so the milestone for a human-genome BAM remains open.

## Human BAM benchmark (chromosome 20, 2026-11-26)

Input: `HG00096.chrom20.ILLUMINA.bwa.GBR.low_coverage.20120522.bam` from the 1000 Genomes
phase 3 FTP (318,238,352 bytes, 2,931,500 records). This is real human chromosome 20, not
the whole genome. The whole-genome BAM for the same sample is 15.6 GB, and the link measured
about 0.46 MB/s, so it was not downloaded.

One block in the first download was damaged in transit. `tptbio` rejected it with an ISIZE
mismatch, and samtools failed to inflate it too. The block was re-fetched from the server,
and all 16,284 BGZF blocks then pass their CRC and length checks.

Environments differ. tptbio is a release build on native Windows, reading the file from NTFS.
samtools 1.22.1 runs in WSL Ubuntu and reads a copy on the Linux filesystem. Each operation
ran five times, and the median is shown.

| Operation | tptbio | samtools | tptbio / samtools |
|---|---|---|---|
| view (full file, SAM text to null) | 4.438 s (4.390-4.521) | 3.190 s (3.108-4.355) | 1.39 |
| index (BAI) | 0.965 s (0.920-1.001) | 2.360 s (1.728-2.389) | 0.41 |
| region `20:10000000-10100000` | 0.355 s (0.340-0.397) | 2.156 s (2.001-2.195) | 0.16 |

The region timings include WSL start-up for samtools. Output agreement, with header lines
excluded: the full view has 2,931,500 records with the same SHA-256 from both tools, and the
region has 4,813 records with the same SHA-256.

So tptbio builds an index faster and answers region queries faster than samtools here, but
full-file SAM output is slower. The genome-scale milestone is still open.

## View after parallel output (2026-11-26)

The first `tptbio seq view` run on this file used one thread to write output. Two changes
follow, and the output is unchanged:

- The BGZF writer compresses blocks in batches of 64 on the rayon pool, writing them in order.
  `tptbio seq view -b` went from 35 s to about 5 s on the same file.
- `seq view` reads records in batches of 4096 raw bodies, and parses and formats each batch on
  the rayon pool. SAM output went from about 4.0 s to about 2.2 to 2.7 s, against samtools'
  3.2 s median.

Checks: the SAM record stream has the same SHA-256 as before (2,931,500 records). The BAM
output passes `samtools quickcheck`, and samtools reads 2,931,500 records from it with the
same record hash. A new unit test checks that batched output matches record-by-record output.
