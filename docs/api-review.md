# API review (read-only)

Scope: the public surface of `tpt-bio-core`, `tpt-bio-seq`, `tpt-bio-vision` and
`tpt-bio-pipeline`, plus the Python binding where it differs from Rust. This is a
review for the 0.1.0 release, not a change list. No library code was changed. Line
numbers refer to the tree at commit `c11a2c4` and may drift as other work lands.

Items are grouped by crate. "Breaking" means the change would break existing callers,
so it is best made before 0.1.0.

## Cross-crate

1. **No public enum is `#[non_exhaustive]`.** The only occurrence is a `finish_non_exhaustive`
   call in `tpt-bio-seq/src/ops.rs:78`. Adding a variant to any error or data enum is
   breaking. Affected: `AminoAcid` (`tpt-bio-core/src/amino.rs:7`), `Base`
   (`nucleotide.rs:23`), `Iupac` (`nucleotide.rs:97`), `SeqError` (`tpt-bio-seq/src/error.rs:8`),
   `FormatError` (`tpt-bio-vision/src/io/mod.rs:38`), `RenderError`
   (`tpt-bio-vision/src/render3d/mod.rs:33`), `Value` (`tpt-bio-pipeline/src/execute.rs:48`)
   and `RunError` (`execute.rs:89`). Breaking if left as is: yes.
2. **Eleven error types across four crates, with no shared trait or umbrella type.**
   Core: `RegionParseError`, `KmerError`, `StatsError`. Seq: `SeqError`. Vision:
   `ImageError`, `FormatError`, `FilterError` (`filter.rs:12`), `GraphCutError`
   (`graphcut.rs:19`), `ShapeError` (`shape3d.rs:22`), `ZernikeError` (`zernike.rs:18`),
   `RenderError`. Pipeline: `DagError` (`lib.rs:81`), `ExecError<E>` (`lib.rs:108`),
   `RunError`, `DefinitionError` (`definition.rs:399`), `TypedError` (`typed.rs:51`),
   `CheckpointError<E>` (`checkpoint.rs:66`). A caller using several crates needs one
   conversion per pair. Recommend one error type per crate, with per-module variants,
   before 0.1.0.
3. **Experimental code is public with no marking.** `stats::tsne` and `stats::umap`
   (`tpt-bio-core/src/stats.rs:18-19`) are described in todo.md as simplified and not
   verified against references. `vision::http` is plain HTTP only, and SigV4 signatures
   go in cleartext (`tpt-bio-vision/src/http.rs:1-12`). Mark these `#[doc(hidden)]`, or
   document them as experimental, before they become part of the semver surface.
4. **Third-party types in public signatures.** `toml::Table` and `toml::ser::Error` appear
   in `tpt-bio-pipeline` (`definition.rs:68`, `264`, `306`). Callers then depend on the
   `toml` version too.
5. **Crate-level docs are incomplete.** The vision crate doc (`tpt-bio-vision/src/lib.rs`)
   does not list `cache`, `conv3d`, `graphcut`, `http`, `shape3d`, `spatial_features`,
   `texture` or `zernike`, though all are `pub mod`.

## tpt-bio-core

6. **Three conventions for invalid input.** Distributions return `NaN` for bad parameters
   (`stats.rs:5`). `qc` and `tsne` return `StatsError` (`qc.rs:18`, `36`; `tsne.rs:62`).
   `gc_content` returns `Option` (`stats.rs:68`), while `reverse_complement` returns
   `Result<_, InvalidSymbol>` (`nucleotide.rs:250`). Pick one rule. Breaking: yes.
7. **`kmer::canonical` exists twice.** A free function (`kmer.rs:66`) and a method on the
   k-mer iterator (`kmer.rs:113`). Harmless, but the docs should point from one to the other.

## tpt-bio-seq

8. **Reader method names differ.** `FastqReader::read_record(&mut FastqBuf) -> Result<bool>`
   (`fastq.rs:162`) is the only reader that fills a caller buffer and returns a bool.
   The others use `next_record() -> Result<Option<T>>`: BAM (`bam.rs:550`), SAM
   (`sam.rs:294`), VCF (`vcf.rs:365`), BCF (`bcf.rs:935`), FASTA (`fasta.rs:54`).
   Breaking: yes if renamed.
9. **Two ways to finish VCF and BCF writers.** `VcfWriter` has `into_inner`
   (`vcf.rs:447`) and `finish` (`vcf.rs:458`); `BcfWriter` has the same pair
   (`bcf.rs:996`, `1015`). Every other writer has only `finish` (`bam.rs:613`,
   `sam.rs:345`, `fastq.rs:341`, `fasta.rs:143`, `bgzf.rs:359`). Keep `finish`.
10. **Region queries have several shapes.** `bai::fetch(reader, index, ref_id: usize,
    beg: u64, end: u64)` (`bai.rs:279`) takes an ID and 0-based half-open coordinates.
    `TabixIndex::query_chunks(chrom: &str, ...)` returns `Result` (`tabix.rs:320`), while
    `query_by_id` returns a plain `Vec` (`tabix.rs:330`). The free functions
    `tabix::fetch_text` (`tabix.rs:692`) and `bcf::fetch_lines` (`bcf.rs:1063`) are
    separate again. Rust has no string region API (todo.md, Phase 1), although
    `core::coords::Region` exists. Recommend one region type from core and one fetch per
    format. Breaking: yes.
11. **Index naming.** `bai::Index` is re-exported at the crate root as `Index`
    (`tpt-bio-seq/src/lib.rs`), next to `TabixIndex`. The root docs do not say which format
    `Index` is. Rename to `BaiIndex`.
12. **`reg2bin` and `reg2bins` disagree.** `reg2bin(beg: i64, end: i64) -> u16`
    (`bam.rs:28`) and `reg2bins(beg: u64, end: u64) -> Vec<u32>` (`bai.rs:222`) use
    different integer types and return shapes for one concept.
13. **The `Result` alias shadows `std::result::Result`.** `pub type Result<T, E = SeqError>`
    (`error.rs:70`) is re-exported by `use tpt_bio_seq::*`. The vision crate uses
    `FormatResult` instead. Rename to `SeqResult` or drop the alias.
14. **`FastqReader` cannot open a path.** `fastq::open_reader(path)` is a free function
    (`fastq.rs:353`), while `FastqReader::new` takes any reader (`fastq.rs:142`). Add
    `FastqReader::open`, or document the pattern.
15. **`records_to_jsonl` cannot fail; the other exporters can.** `table.rs:206` returns
    `String`, while `records_to_csv` (`table.rs:132`) and `records_to_table`
    (`table.rs:35`) return `Result`. Make all three the same.
16. **Public fields on core records.** `Flags(pub u16)` (`record.rs:11`) and every field of
    `Record` (`record.rs:210-233`) are public. This is convenient, but any new field
    breaks struct literals. Consider `#[non_exhaustive]` on `Record`, with a constructor.
17. **Stringly-typed errors.** `SeqError::Unsupported(String)` and `NotFound(String)`
    (`error.rs:21`, `23`) hold free text, so callers cannot match on "reference not
    found". Consider a dedicated variant.

## tpt-bio-vision

18. **Two Gaussian blurs with different types.** `filter::gaussian_blur(input: &[f32],
    width, height, sigma: f32)` (`filter.rs:269`) takes raw slices and returns
    `FilterError`. `spatial::gaussian_smooth(&Image<f64>, sigma: f64)` (`spatial.rs:83`)
    returns `ImageError`. Keep one, on `Image`, and deprecate the other before 0.1.0.
19. **Boolean parameters for TIFF flavour.** `save_tiff(path, img, description, big: bool)`
    (`io/write.rs:71`), `encode_tiff(img, description, big)` (`write.rs:53`) and
    `encode_ome_tiff(volume, (z, c, t), name, big: bool)` (`io/ome.rs:560`). An enum with
    `Classic` and `Big` variants is clearer. Breaking: yes.
20. **Readers and writers take different inputs.** Readers take a `ByteSource`
    (`io/png.rs:19`, `io/tiff.rs:25`). Writers take a `&Path` (`io/write.rs:41`, `71`).
    `read_png` returns `Planes` (`io/mod.rs:91`), not an `Image`. Pick one model: both take
    `ByteSource` and `Write`, or both take paths.
21. **TIFF writer accepts only `u8`, `u16` and `f32`** (`io/sample.rs:68-70`). A label image
    with more than 65,535 objects cannot be written as TIFF without clipping. The CLI's
    `vision segment` writes 16-bit labels. Document the cap, or add a 32-bit label path.
22. **Dimension appears in some names only.** `label_2d` returns `Result` (`label.rs:42`),
    while `label_3d` is infallible (`label.rs:55`). `shape_features_2d` (`features.rs:84`)
    names its dimension, while `intensity_features` (`features.rs:192`) does not, although
    it is also 2D only. `Connectivity2` and `Connectivity3` (`label.rs:13`, `22`) use digit
    suffixes; `Connectivity2D` and `Connectivity3D` read better.
23. **`Volume` constructors take different collection types.** `Volume::new(Vec<Image<f32>>,
    spacing)` (`render3d/volume.rs:26`) takes the vector by value; `Volume::from_u16(&[Image<u16>],
    spacing)` (`volume.rs:44`) takes a slice. Add `from_f32(&[Image<f32>])`, or accept
    `impl IntoIterator`.
24. **Background colour is set two ways.** `RenderSettings::background` (`render3d/raycast.rs:23-26`)
    is a field, but `rasterise_mesh(mesh, camera, colour, background)`
    (`render3d/rasterise.rs:24-29`) takes it as an argument. Use one settings struct for both.
25. **An internal argument is exposed.** `Camera::ray(px, py, ortho_back: f64)`
    (`render3d/camera.rs:135`) takes an orthographic back-offset that only the renderer can
    choose sensibly. Make it private, or give a public `ray(px, py)` that uses the default.
26. **Stringly-typed render errors.** `RenderError::InvalidParameter(&'static str)`
    (`render3d/mod.rs`) covers every range check. Callers cannot tell which parameter failed.
    Add a variant per parameter, or a `kind` enum.
27. **Error conversions are incomplete.** `FormatError` converts from `ImageError`
    (`io/mod.rs:82`), but `FilterError` (`filter.rs:12`) has no `From` impl into
    `ImageError` or `FormatError`. Code that mixes the filter and image functions must
    convert by hand. Add the conversions, or unify the types (item 2).
28. **Zarr writing is v2 only.** The reader handles v2 and v3 (`io/zarr.rs:11`), and the
    writers are `write_group_v2`, `write_array_v2` and `write_image_v2`
    (`zarr.rs:684`, `705`, `799`). The `_v2` suffix is accurate, but the asymmetry should
    be stated in the crate docs.

## tpt-bio-pipeline

29. **`Pipeline::run` has no cancel or progress.** `Pipeline::run(&self)` (`execute.rs:120`)
    takes neither a cancel token nor a progress callback. `Dag` has four executors
    (`lib.rs:268`, `279`, `317`, `334`), and `execute_checkpointed` (`checkpoint.rs:161`).
    Add `run_with(cancel, progress)` so TOML users get the same controls. Breaking: no.
30. **Cancellation is a bare `&AtomicBool`** in every `Dag` executor. A small `CancelToken`
    type would document the intent and leave room for change without breaking callers.
31. **Five overlapping `Dag` executors.** `execute`, `execute_observed`, `execute_parallel`
    and `execute_parallel_observed` differ only by two options (parallel, progress). An
    options struct would replace them. Keep the old names as thin wrappers until 0.1.0.
32. **`RunError` and `ExecError` describe the same failure.** `ExecError<E>` is generic;
    `RunError` is concrete, with `String` fields (`execute.rs:89-98`). `Pipeline::run`
    converts one into the other. Use one error for the run, and keep node errors generic.
33. **Stringly-typed node errors.** `parse_values(params) -> Result<_, String>`
    (`definition.rs:537`) and the node evaluators in `execute.rs` return `String`. A
    node-level error enum would let callers tell a bad parameter from bad data.
34. **`Value` is exhaustive** (`execute.rs:48`). Adding an output kind, such as a VCF
    records value, breaks every `match`. See item 1.
35. **Checkpoint keys are node names only** (todo.md, Phase 5). A changed input or parameter
    reuses the old result. This is a semantic issue rather than a signature one, but the docs
    should say so before anyone relies on resume.

## Python binding (`bindings/tpt-bio-py`)

36. **Same name, different meaning.** Rust `Pipeline::from_toml(text: &str)`
    (`definition.rs:288`) parses a string. Python `Pipeline.from_toml(path: str)`
    (`bindings/tpt-bio-py/src/lib.rs:481`) reads a file. Rename the Python method to
    `from_file`, or add `from_toml_file`. Breaking for Python users.
37. **Region strings exist only in Python and the CLI.** `IndexedBam.fetch("chr1:1000-2000")`
    takes a string, while Rust takes IDs (item 10). The parser lives in the binding and the
    CLI. Move it into `core::coords::Region` and use it from both.
38. **Only the sequence utilities release the GIL** (todo.md, Phase 5). The image and pipeline
    calls hold it for their whole run, which blocks other Python threads.

## Found while reading (a bug, not an API issue)

39. **`tptbio seq sort` writes a duplicate `SO` tag** when the input header already has one.
    `tpt-bio-cli/src/main.rs:661-668` appends `SO:coordinate` to the `@HD` line without
    removing an existing `SO:unsorted`. Reproduced with a SAM whose `@HD` line says
    `SO:unsorted`: the output header is `@HD VN:1.6 SO:unsorted SO:coordinate`. The fix
    belongs in the CLI. It was not made here, since this review is read-only.

## Suggested order before 0.1.0

1. Error types (2) and `#[non_exhaustive]` (1), since they touch everything.
2. Breaking renames: items 8, 9, 10, 12, 13, 19, 36.
3. Mark or move experimental items (3, 5).
4. Additive changes: items 29, 30, 37.
