# CRAM conformance against samtools

This records how the CRAM reader and writer were checked against samtools 1.22.1
(htslib 1.22.1, which links libhtscodecs 2.1.6). samtools ran in WSL Ubuntu
(`/usr/bin/samtools`). Nothing in the crate links against it.

## What was checked

| Check | Result |
|---|---|
| samtools decodes our CRAM 3.0 file | Yes. 250 records. Identical to samtools' own decode of the same data in every column, RG, NM and MD. |
| Our reader decodes a samtools CRAM 3.0 file | Yes. 250 records. Identical in every column and tag, including NM and MD, which are recomputed from the reference when the file omits them. |
| Our writer ends the file with samtools' EOF container | Yes, byte for byte. Before this check the writer used a substitute marker, which samtools reported as missing. |
| We decode a samtools CRAM 3.1 file (default codecs) | Yes. 250 records, identical to samtools in every column and tag. rANS Nx16 with PACK and CAT, and order-0 blocks. |
| We decode samtools CRAM 3.1 with 32-state rANS Nx16 (N32) | Yes. 3000 records, same digest as samtools' decode. Blocks with flags 0x04 (N32, order 0) and 0x05 (N32, order 1). |
| We decode samtools CRAM 3.1 with PACK, RLE, N32 and order 1 | Yes. 3000 records, same digest. Blocks with flags 0x80, 0x85 (PACK, N32, order 1), 0xA0 and 0xE0 (PACK, RLE, CAT). |
| We decode a reference-free samtools CRAM 3.1 file | Yes. Opened with `CramReader::new` (no reference). 250 records, same digest as samtools. |
| We decode a samtools CRAM 3.1 file whose read names use TOK3 (block method 8) | Yes. `samtools_v31_tok3.cram`, samtools default codecs (`version=3.1`), Illumina-style names. 250 records, identical to samtools in every column and tag. The file has two method 8 blocks, one per slice. |
| TOK3 decoder against htscodecs' `tok3_decode_names` | Ported from `tokenise_name3.c`. Hand-built unit tests cover duplicates, zero-padded numbers, constant type streams, truncation, stored arithmetic streams and bad flags. Both forms are decoded: rANS Nx16 (flag 0) and arithmetic (flag 1, see below). |
| We decode a samtools CRAM 3.1 file whose TOK3 names use the arithmetic flag (flag 1) | Yes. `samtools_v31_tok3_arith.cram`, written with `--output-fmt-option use_arith=1`. 250 records; the projected, sorted records equal samtools' own decode. The two name blocks have flag 1. The 23 external data-series blocks are method 6 (arithmetic) and decode as well. |
| Arithmetic coder (htscodecs `arith_dynamic.c`) for TOK3 streams and block method 6 | Ported in `tpt-bio-seq/src/cram/arith.rs`: order 0 and 1 models, run-length, bit-packing, stored, and striped streams. Checked end to end by the arithmetic fixture, which uses order 0 and 1, packed and packed run-length streams, stored streams, and striped streams. Unit tests cover packing, striping, stored streams and truncation. |
| Block level: every rANS Nx16 block of four samtools files | Yes. 156 payloads. Our output is byte-identical to the libhtscodecs decode that samtools uses. |
| We decode a samtools CRAM 3.1 file whose blocks are bzip2 (method 2) | Yes. `samtools_v31_bzip2.cram`, written with `--output-fmt-option use_bzip2=1`. 250 records, identical to samtools in every column and tag. The file holds bzip2 streams. The decoder is in `tpt-bio-seq/src/cram/bzip2.rs`, checked on libbzip2 fixtures (multi-block, concatenated, corrupt and truncated input) in `tests/fixtures/bzip2`. |
| `cargo fuzz run cram_read` for 60 seconds | Ran with no crash (8606 runs). |

## Transforms samtools produced

Block method 5 (rANS Nx16) flags seen in samtools output, across about 40 generated files
(levels 1, 2 to 9, use_tok on and off, use_lzma, use_bzip2, use_rans, use_tok, lossy_names,
and inputs of 250 to 10000 records):

| Flags | Meaning | Seen |
|---|---|---|
| 0x00 | order 0 | yes |
| 0x04 | N32, order 0 | yes (large blocks only) |
| 0x05 | N32, order 1 | yes (same) |
| 0x80 | PACK | yes |
| 0x85 | PACK, N32, order 1 | yes |
| 0xA0 | PACK, CAT | yes |
| 0xE0 | PACK, RLE, CAT | yes |
| 0x01 | order 1 without N32 | no |
| 0x08 | STRIPE | no |
| 0x40 | RLE with rANS data (no CAT) | no |

Order 1 was only produced together with N32, and only in large blocks (tens of KB). The
0x85 block in `samtools_v31_pack_n32.cram` has a compressed order-1 table (header byte 0xa1).
Every order-1 block seen uses table shift 10 (header high nibble 0xa). Shift 12 is accepted
but no sample uses it. STRIPE and RLE without CAT were never produced, so those are tested
only with streams from the htscodecs encoder (see the unit tests below).

## What was not checked

- **STRIPE, RLE with rANS data, and plain order 1 (without N32)** from samtools. Verified only
  against streams that the htscodecs library produced (unit tests).
- **Writing CRAM 3.1 and writing order 1 or N32**: our writer produces 3.0 only. There is no
  order-1 or N32 encoder in the crate.
- **Arithmetic streams with order 2**: reported as `Unsupported`. samtools does not write
  order 2. The bzip2 (external) form of arithmetic streams is decoded, but only checked against
  libbzip2 fixtures, not a samtools file. Plain run-length without packing, unpacked stripes and the
  no-size flag are ported but not covered by a samtools fixture.
- **The `archive` profile** (`-O cram,version=3.1,archive`) wrote the same block methods as
  `use_arith=1` for this file (2 flag-1 name blocks, 23 method 6 blocks). Its decode also matches
  samtools (250 records, checked by hand; it is not kept as a fixture).
- **Writing TOK3**: not implemented. The writer never emits method 8.
- **Which names samtools tokenizes**: TOK3 appeared for Illumina-style names
  (`HWI-ST1234:8:1101:...`). It did not appear for `r00005`-style names, nor for a set of
  mixed names with mate duplicates, where samtools chose rANS Nx16. Those names were probed
  with samtools 1.22.1 default settings, and the choice is samtools' own.
- **Other TOK3 token kinds**: the only samtools-written TOK3 data checked is the Illumina-style
  fixture. The CHAR, NOP and match-by-digits paths are ported from the reference but no test
  exercises them.
- **BAM to CRAM to BAM equality through samtools**: not run.

## Bugs found by this check and fixed

1. **rANS frequency tables summing to 4095.** htscodecs accepts these and fills the last
   slot with the preceding symbol. We rejected them.
2. **Tag map keys are ITF8 integers.** The key packs the two-letter tag name and its
   type into one integer. We read three raw bytes, which shifted every later field.
3. **RG was dropped.** The read-group series is now mapped to the header's `@RG` IDs.
4. **EOF container.** The writer now emits the canonical 38-byte CRAM 3.0 EOF container.
5. **MD and NM.** samtools recomputes them on decode, so the decoder now does too.
6. **rANS Nx16.** The module was rebuilt from the htscodecs reference. The varints are big-endian, the alphabet comes before the frequencies, the states are read after the table, and the PACK, RLE and STRIPE layers are applied.
7. **Order 1 and N32 were unsupported.** They are now decoded, including order-1 tables
   stored compressed with order 0.

## Reproducing the fixtures

Fixtures live in `tpt-bio-seq/tests/data/samtools/`. The 250 records come from
`sorted.bam` (`samtools view -h sorted.bam > in.sam`). The reference is `ref.fa`.

```sh
# CRAM 3.0 (earlier check)
samtools view -C -T ref.fa --output-fmt-option version=3.0 in.sam -o samtools_v30.cram

# CRAM 3.1, default codecs
samtools view -C -T ref.fa --output-fmt-option version=3.1 in.sam -o samtools_v31.cram

# CRAM 3.1, default codecs, with read names that samtools stores as TOK3 (block method 8).
# names_illum.sam is in.sam with each record renamed (records only, header kept):
#   awk 'BEGIN{OFS="\t"} {$1=sprintf("HWI-ST1234:8:%d:%d:%d:%d#0/%d", 1101+int(NR/7), 1000+NR*37%5000, 2000+NR*91%9000, NR, 1+NR%2); print}'
samtools view -C -T ref.fa --output-fmt-option version=3.1 names_illum.sam -o samtools_v31_tok3.cram
# The expected records are samtools' decode of that file, projected and sorted like the others.
samtools view -T ref.fa samtools_v31_tok3.cram   # -> expected_tok3_records.txt after projection

# 3000 records: the 250 records, 12 times, each name suffixed with _1 .. _12.
# use_tok=0 keeps these fixtures on the rANS path, so their digests are unchanged by TOK3 support.
# Arithmetic TOK3: use_arith=1 writes flag-1 name blocks and arithmetic (method 6) data series.
samtools view -C -T ref.fa --output-fmt-option version=3.1 --output-fmt-option use_arith=1 names_illum.sam -o samtools_v31_tok3_arith.cram
samtools view -T ref.fa samtools_v31_tok3_arith.cram   # -> expected_tok3_arith_records.txt after projection
(grep '^@' in.sam; for i in $(seq 1 12); do grep -v '^@' in.sam | awk -v i=$i 'BEGIN{OFS="\t"}{$1=$1"_"i; print}'; done) > c12.sam
samtools view -C -T ref.fa --output-fmt-option version=3.1 --output-fmt-option level=1 --output-fmt-option use_tok=0 c12.sam -o samtools_v31_n32.cram
samtools view -C -T ref.fa --output-fmt-option version=3.1 --output-fmt-option level=5 --output-fmt-option use_tok=0 c12.sam -o samtools_v31_pack_n32.cram

# Reference-free (no -T)
samtools view -C --output-fmt-option version=3.1 --output-fmt-option no_ref=1 in.sam -o samtools_v31_noref.cram
```

The committed tests compare record digests rather than full text, to keep the repository
small. The digest is FNV-1a 64 over the records, each projected to its 11 mandatory columns
plus the tags in sorted order, then sorted, each followed by a newline. The digest was
computed from `samtools view -h` of each file. Both 3000-record files have digest
`de3d334c10dfaf0d`, and the reference-free file has `b71334f4f9810ce9`.

The rANS Nx16 unit tests use streams from the htscodecs library in WSL, built with
`gcc` against `/usr/lib/x86_64-linux-gnu/libhtscodecs.so.2`. The order-0 and order-1
32-state bodies come from its exported `rans_compress_O0_32x16` and `rans_compress_O1_32x16`,
because the public `rans_compress_to_4x16` drops the 32-state flag for inputs of 1000 bytes
or less. The 4-way order-1, RLE and STRIPE vectors come from `rans_compress_to_4x16`.

The example programs in `tpt-bio-seq/examples/` decode and encode CRAM for manual checks:

```sh
cargo run -p tpt-bio-seq --example cram_oracle_dump -- in.cram ref.fa > ours.sam
cargo run -p tpt-bio-seq --example cram_oracle_dump -- in.cram - > ours.sam   # no reference
cargo run -p tpt-bio-seq --example cram_oracle_nx16 -- block.bin RAW_SIZE out.bin
cargo run -p tpt-bio-seq --example cram_oracle_write -- in.sam ref.fa out.cram
```

`cram_oracle_nx16` decodes one block payload, starting at its flags byte, so block-level
output can be compared with the htscodecs library.

## Tests

`tpt-bio-seq/tests/cram_samtools.rs`:

- `decodes_samtools_cram_30_to_samtools_records`
- `decodes_samtools_cram_31_to_samtools_records`
- `decodes_samtools_cram_31_n32_blocks`
- `decodes_samtools_cram_31_pack_n32_order1_blocks`
- `decodes_reference_free_samtools_cram_31`
- `decodes_samtools_cram_31_tok3_names_to_samtools_records`
- `decodes_samtools_cram_31_arithmetic_tok3_names`
- `writer_ends_with_samtools_eof_container`

`tpt-bio-seq/src/cram/tok3.rs` (unit tests, hand-built blobs; expected values derived from the
htscodecs `decode_name` logic):

- `decodes_zero_padded_numeric_tokens`, `decodes_duplicate_name`, `decodes_constant_type_stream`
- `rejects_length_mismatch_and_truncation`, `decodes_stored_arithmetic_streams`, `rejects_unknown_arithmetic_flag`

`tpt-bio-seq/src/cram/rans_nx16.rs` (unit tests, streams from htscodecs):

- `order1_stream_from_htscodecs`: 4-way order 1, shift 10.
- `order0_32_states_from_htscodecs` and `order1_32_states_from_htscodecs`: N32 bodies.
- `rle_stream_from_htscodecs` and `stripe_stream_from_htscodecs`.
- `order1_shift_other_than_10_or_12_is_malformed`.

All of them run without samtools or htscodecs installed.
