# Sequence quickstart

This tutorial reads alignments from a SAM file, sorts them into a BAM, builds an
index, queries a genomic region and exports the records. Every input is generated
below, so you need no downloaded data.

Everything here is covered by the `tptbio` command-line tool and by the library
crate `tpt-bio-seq`. The commands assume `tptbio` is on your `PATH`. Build it with:

```sh
cargo build -p tpt-bio-cli
# the binary is target/debug/tptbio (target\debug\tptbio.exe on Windows)
```

## 1. Generate a small SAM file

SAM is the text form of alignments. The file below has a header with two
references and 200 reads on `chr1`. The reads are not sorted, which is what
real aligners often produce.

```sh
printf '@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:100000\n@SQ\tSN:chr2\tLN:50000\n' > reads.sam
for i in $(seq 1 200); do
  pos=$(( (i * 487) % 99000 + 1 ))
  printf 'read%03d\t0\tchr1\t%d\t60\t50M\t*\t0\t0\tACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTAC\tIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII\n' "$i" "$pos" >> reads.sam
done
tptbio seq view reads.sam -H
```

The last command prints the header only (`-H`). Each record line has the eleven
mandatory SAM columns: name, FLAG, reference, 1-based position, MAPQ, CIGAR, mate
reference, mate position, template length, bases and qualities.

## 2. Sort into a BAM

An index needs coordinate-sorted input. `seq sort` sorts in memory and writes BAM.
`seq view -b` only converts the format and does not sort, so use it only for input that is already sorted.

```sh
tptbio seq sort reads.sam -o sorted.bam
```

## 3. Build the index

```sh
tptbio seq index sorted.bam
```

This writes `sorted.bam.bai`. `seq index` fails with `input is not coordinate-sorted`
if the BAM is out of order, so run step 2 first.

## 4. Query a region

Regions use 1-based, inclusive coordinates, as in `samtools`:

```sh
tptbio seq view sorted.bam -r chr1:1000-2000
```

The output starts with the header, then every record that overlaps the region.
Use `-i` to point at an index with another name. Querying a reference that is not
in the header is an error (`not found: reference chr3`).

## 5. Summarize and export

```sh
tptbio seq flagstat sorted.bam
tptbio seq export sorted.bam -f csv -o reads.csv
tptbio seq export sorted.bam -f jsonl -o reads.jsonl
head -n 2 reads.csv
```

`seq export` writes CSV or JSON Lines. Parquet is not supported, because it needs
Arrow. VCF and BCF are read by `seq view` and by the library, but `seq export`
handles alignment records only.

## Library version

The same steps are in `tpt-bio-seq/examples/seq_quickstart.rs`. It writes a BAM,
builds the index in memory, runs `bai::fetch` for `chr1:1000-2000` and exports
the result with `records_to_csv`:

```sh
cargo run -p tpt-bio-seq --example seq_quickstart
```

The core calls are:

```rust
let index = bai::build_index(&mut open_bam(&bam_path)?)?;
let index = Index::read_from(BufReader::new(File::open(dir.join("reads.bam.bai"))?))?;
let mut reader = open_bam(&bam_path)?;
// 0-based half-open coordinates: chr1:1000-2000 is [999, 2000).
let hits = bai::fetch(&mut reader, &index, 0, 999, 2000)?;
let csv = records_to_csv(&hits, reader.header())?;
```

## What is not covered

- CRAM, VCF and BCF are not covered here. The library has readers for them
  (`tpt_bio_seq::cram`, `vcf`, `bcf`), and the CLI reads VCF and BCF with `seq view`.
- Coordinate sorting is in memory only. There is no external merge sort yet.
