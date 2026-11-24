//! `tptbio`: the unified command-line tool.
//!
//! Implemented: `seq view` (SAM, BAM, VCF, bgzipped VCF and BCF), `seq index`,
//! `seq flagstat`, `seq sort`, `seq export`, `stats normalize` and `stats pca` on CSV
//! tables, `vision segment` and `vision features` on PNG and TIFF images, `pipeline run`
//! on TOML pipeline definitions, `view` (an image data summary), plus `completions` and
//! `man`.
//!
//! Settings come from the command line, then `tptbio.toml` (see [`config`]), then the
//! built-in defaults. Exit status: 0 on success, 1 for runtime errors (bad input, I/O,
//! bad configuration), and 2 for usage errors (clap's own parse errors, unsupported
//! file names, and invalid pipeline definitions). Log messages go to standard error.
//! `-q` shows errors only, and `-v` adds debug output.

mod config;
mod output;
mod variants;

use std::ffi::OsStr;
use std::fmt::{self, Write as _};
use std::fs::{self, File};
use std::io::{self, BufRead, BufReader, BufWriter, Read, Write};
use std::path::{Path, PathBuf};
use std::process::ExitCode;
use std::sync::atomic::{AtomicU8, Ordering};

use clap::{ArgAction, Args, CommandFactory, Parser, Subcommand, ValueEnum};
use clap_complete::Shell;
use config::{LogLevel, OutputFormat};
use tpt_bio_core::Region;
use tpt_bio_core::stats::normalize::{self, Method};
use tpt_bio_core::stats::pca::Pca;
use tpt_bio_core::stats::{Matrix, StatsError};
use tpt_bio_pipeline::{Pipeline, Value};
use tpt_bio_seq::bai::{self, Index};
use tpt_bio_seq::{BamReader, BamWriter, Header, SamReader, SamWriter, SeqError};
use tpt_bio_vision::features::{intensity_features, shape_features_2d};
use tpt_bio_vision::image::{Image, ImageError};
use tpt_bio_vision::io::Planes;
use tpt_bio_vision::io::png::read_png;
use tpt_bio_vision::io::tiff::read_tiff;
use tpt_bio_vision::io::write::save_tiff;
use tpt_bio_vision::label::{Connectivity2, Labeling, label_2d};
use tpt_bio_vision::objects::{ObjectFilter, filter_objects};
use tpt_bio_vision::source::FileSource;
use tpt_bio_vision::threshold::{apply_threshold, otsu_threshold};
use variants::Kind;

/// The active log threshold, stored as a [`LogLevel`] discriminant.
static LOG_LEVEL: AtomicU8 = AtomicU8::new(LogLevel::Info as u8);

fn log(level: LogLevel, msg: impl fmt::Display) {
    if LOG_LEVEL.load(Ordering::Relaxed) >= level as u8 {
        eprintln!("tptbio: {msg}");
    }
}

fn log_warn(msg: impl fmt::Display) {
    log(LogLevel::Warn, msg);
}

fn log_info(msg: impl fmt::Display) {
    log(LogLevel::Info, msg);
}

fn log_debug(msg: impl fmt::Display) {
    log(LogLevel::Debug, msg);
}

#[derive(Debug, Parser)]
#[command(name = "tptbio", version, about = "High-performance bioinformatics and imaging toolkit")]
struct Cli {
    /// Print only errors. Conflicts with --verbose.
    #[arg(short, long, global = true, conflicts_with = "verbose")]
    quiet: bool,
    /// Print debug messages as well. Repeat for no extra effect.
    #[arg(short, long, action = ArgAction::Count, global = true)]
    verbose: u8,
    /// Configuration file. Defaults to `tptbio.toml` in the current directory, if present.
    #[arg(long, global = true, value_name = "PATH")]
    config: Option<PathBuf>,
    /// Layout of tabular output (csv, jsonl or table). Overrides the configuration file.
    #[arg(long, global = true, value_enum, value_name = "FORMAT")]
    output_format: Option<OutputFormat>,
    #[command(subcommand)]
    command: Command,
}

#[derive(Debug, Subcommand)]
enum Command {
    /// Sequencing data: SAM, BAM, VCF and BCF view, indexing, flagstat, sort and export.
    Seq {
        #[command(subcommand)]
        cmd: SeqCmd,
    },
    /// Biological image analysis: segmentation and per-object features of PNG and TIFF.
    Vision {
        #[command(subcommand)]
        cmd: VisionCmd,
    },
    /// Workflow execution from TOML pipeline definitions.
    Pipeline {
        #[command(subcommand)]
        cmd: PipelineCmd,
    },
    /// Statistical analysis of CSV tables: normalization and PCA.
    Stats {
        #[command(subcommand)]
        cmd: StatsCmd,
    },
    /// Summarize a PNG or TIFF image: size, pixel type and per-channel statistics.
    ///
    /// This is a data summary, not a viewer. The interactive 3D viewer is not in the
    /// CLI yet.
    View {
        /// PNG or TIFF image (.png, .tif or .tiff).
        input: PathBuf,
    },
    /// Print a shell completion script to standard output.
    Completions {
        /// Shell to generate for.
        #[arg(value_enum)]
        shell: Shell,
    },
    /// Write one man page per command into a directory.
    Man {
        /// Directory to write into. Created if missing.
        dir: PathBuf,
    },
}

#[derive(Debug, Subcommand)]
enum SeqCmd {
    /// Print SAM, BAM, VCF or BCF records, optionally restricted to a region.
    View(ViewArgs),
    /// Build a BAI index next to a coordinate-sorted BAM.
    Index {
        /// Coordinate-sorted BAM file.
        input: PathBuf,
        /// Output path. Defaults to the input path with `.bai` appended.
        #[arg(short, long)]
        output: Option<PathBuf>,
    },
    /// Count records by flag category, in the style of `samtools flagstat`.
    Flagstat {
        /// SAM or BAM file.
        input: PathBuf,
    },
    /// Sort records by coordinate, in memory, and write BAM.
    Sort {
        /// SAM or BAM file.
        input: PathBuf,
        /// Output BAM file.
        #[arg(short, long)]
        output: PathBuf,
    },
    /// Export SAM or BAM records as CSV or JSON Lines. No Parquet: it needs Arrow.
    Export {
        /// SAM or BAM file.
        input: PathBuf,
        /// Output format. `table` is not supported here and falls back to CSV.
        /// Defaults to the global `--output-format`, then CSV.
        #[arg(short, long, value_enum)]
        format: Option<OutputFormat>,
        /// Output file. Defaults to standard output.
        #[arg(short, long)]
        output: Option<PathBuf>,
    },
}

#[derive(Debug, Args)]
struct ViewArgs {
    /// SAM, BAM, VCF, bgzipped VCF (.vcf.gz) or BCF file. The format is detected.
    input: PathBuf,
    /// Region such as `chr1:1000-2000` (1-based, inclusive). Requires an index.
    #[arg(short, long)]
    region: Option<String>,
    /// Index file to use with `--region`. Defaults to `<input>.bai` for BAM and
    /// `<input>.tbi` for VCF and BCF.
    #[arg(short, long)]
    index: Option<PathBuf>,
    /// Print only the header.
    #[arg(short = 'H', long)]
    header_only: bool,
    /// Write BAM instead of SAM text. SAM and BAM input only.
    #[arg(short = 'b', long)]
    bam: bool,
    /// Output file. Defaults to standard output.
    #[arg(short, long)]
    output: Option<PathBuf>,
}

#[derive(Debug, Subcommand)]
enum VisionCmd {
    /// Threshold one channel, label connected components and write a label image.
    Segment(SegmentArgs),
    /// Measure shape and, optionally, intensity features of each labelled object.
    Features(FeaturesArgs),
}

#[derive(Debug, Clone, Copy, ValueEnum)]
enum ThresholdArg {
    /// Global Otsu threshold over the channel's histogram.
    Otsu,
}

#[derive(Debug, Clone, Copy, ValueEnum)]
enum ConnectivityArg {
    /// Edge-sharing neighbours.
    #[value(name = "4")]
    Four,
    /// Edge- and corner-sharing neighbours.
    #[value(name = "8")]
    Eight,
}

#[derive(Debug, Args)]
struct SegmentArgs {
    /// PNG or TIFF image (.png, .tif or .tiff).
    input: PathBuf,
    /// Output label image, written as 16-bit TIFF. Background is 0.
    #[arg(short, long)]
    output: PathBuf,
    /// Channel to segment, counting from 0.
    #[arg(short, long, default_value_t = 0)]
    channel: usize,
    /// Thresholding method.
    #[arg(short, long, value_enum, default_value_t = ThresholdArg::Otsu)]
    method: ThresholdArg,
    /// Smallest object kept, in pixels.
    #[arg(long, default_value_t = 0, value_name = "PIXELS")]
    min_size: usize,
    /// Drop objects that touch the image border.
    #[arg(long)]
    remove_border: bool,
    /// Pixel neighbourhood used to join foreground pixels into objects.
    #[arg(long, value_enum, default_value_t = ConnectivityArg::Four)]
    connectivity: ConnectivityArg,
}

#[derive(Debug, Args)]
struct FeaturesArgs {
    /// Image to measure intensities from. Alternative to --intensity. Omit both for
    /// shape features only.
    #[arg(conflicts_with = "intensity")]
    input: Option<PathBuf>,
    /// Label image from `vision segment` (TIFF or PNG).
    #[arg(long, value_name = "PATH")]
    labels: PathBuf,
    /// Intensity image with the same size as the labels.
    #[arg(long, value_name = "PATH")]
    intensity: Option<PathBuf>,
    /// Channel of the intensity image, counting from 0.
    #[arg(short, long, default_value_t = 0)]
    channel: usize,
    /// Output CSV, or the format chosen by `--output-format`. Required.
    #[arg(short, long)]
    output: PathBuf,
}

#[derive(Debug, Subcommand)]
enum PipelineCmd {
    /// Run a pipeline definition and print its table outputs. Labels are summarized.
    Run {
        /// Pipeline definition file (TOML).
        definition: PathBuf,
    },
}

#[derive(Debug, Subcommand)]
enum StatsCmd {
    /// Normalize each numeric column, optionally within groups such as plates.
    Normalize(NormalizeArgs),
    /// Project rows onto the top principal components.
    Pca(PcaArgs),
}

#[derive(Debug, Clone, Copy, ValueEnum)]
enum MethodArg {
    /// (x - mean) / standard deviation.
    #[value(name = "zscore")]
    ZScore,
    /// (x - median) / (1.4826 * MAD).
    #[value(name = "robust")]
    Robust,
}

impl From<MethodArg> for Method {
    fn from(m: MethodArg) -> Self {
        match m {
            MethodArg::ZScore => Self::ZScore,
            MethodArg::Robust => Self::Robust,
        }
    }
}

#[derive(Debug, Args)]
struct NormalizeArgs {
    /// CSV file with a header row. Every column except --group must be numeric.
    input: PathBuf,
    /// Centering and scaling rule.
    #[arg(short, long, value_enum, default_value_t = MethodArg::ZScore)]
    method: MethodArg,
    /// Column of group labels, such as a plate name. Statistics are computed per group.
    #[arg(short, long, value_name = "COLUMN")]
    group: Option<String>,
    /// Output CSV. Defaults to standard output.
    #[arg(short, long)]
    output: Option<PathBuf>,
}

#[derive(Debug, Args)]
struct PcaArgs {
    /// CSV file with a header row. Every column must be numeric.
    input: PathBuf,
    /// Number of principal components to keep.
    #[arg(short = 'k', long, default_value_t = 2)]
    components: usize,
    /// Output CSV of scores. Defaults to standard output.
    #[arg(short, long)]
    output: Option<PathBuf>,
}

/// Settings resolved from the command line and the configuration file.
#[derive(Debug)]
struct Settings {
    /// Layout for tabular output. `None` lets each command use its own default.
    output_format: Option<OutputFormat>,
}

fn main() -> ExitCode {
    let cli = Cli::parse();
    let (source, config) = match config::load(cli.config.as_deref()) {
        Ok(loaded) => loaded,
        Err(msg) => {
            eprintln!("tptbio: {msg}");
            return ExitCode::FAILURE;
        }
    };
    // Flags beat the file, and the file beats the default.
    let level = if cli.quiet {
        LogLevel::Error
    } else if cli.verbose > 0 {
        LogLevel::Debug
    } else {
        config.log_level.unwrap_or(LogLevel::Info)
    };
    LOG_LEVEL.store(level as u8, Ordering::Relaxed);
    if let Some(path) = &source {
        log_debug(format_args!("using config {}", path.display()));
    }
    if let Some(threads) = config.threads {
        log_debug(format_args!("threads = {threads} (recorded; the CLI has no thread pool)"));
    }
    let settings = Settings { output_format: cli.output_format.or(config.output_format) };

    let result = match cli.command {
        Command::Seq { cmd: SeqCmd::View(args) } => seq_view(&args, &settings),
        Command::Seq { cmd: SeqCmd::Index { input, output } } => seq_index(&input, output),
        Command::Seq { cmd: SeqCmd::Flagstat { input } } => seq_flagstat(&input),
        Command::Seq { cmd: SeqCmd::Sort { input, output } } => seq_sort(&input, &output),
        Command::Seq { cmd: SeqCmd::Export { input, format, output } } => {
            seq_export(&input, format.or(settings.output_format), output.as_deref())
        }
        Command::Vision { cmd: VisionCmd::Segment(args) } => vision_segment(&args),
        Command::Vision { cmd: VisionCmd::Features(args) } => vision_features(&args, &settings),
        Command::Pipeline { cmd: PipelineCmd::Run { definition } } => {
            pipeline_run(&definition, &settings)
        }
        Command::Stats { cmd: StatsCmd::Normalize(args) } => stats_normalize(&args, &settings),
        Command::Stats { cmd: StatsCmd::Pca(args) } => stats_pca(&args, &settings),
        Command::View { input } => image_view(&input),
        Command::Completions { shell } => write_completions(shell),
        Command::Man { dir } => write_man_pages(&dir),
    };
    match result {
        Ok(()) => ExitCode::SUCCESS,
        Err(CliError::Usage(msg)) => {
            eprintln!("tptbio: {msg}");
            ExitCode::from(2)
        }
        Err(CliError::Seq(e)) => {
            eprintln!("tptbio: {e}");
            ExitCode::FAILURE
        }
        Err(CliError::Runtime(msg)) => {
            eprintln!("tptbio: {msg}");
            ExitCode::FAILURE
        }
    }
}

enum CliError {
    /// Bad arguments or input: exit status 2.
    Usage(String),
    Seq(SeqError),
    /// Runtime failure: exit status 1.
    Runtime(String),
}

impl From<ImageError> for CliError {
    fn from(e: ImageError) -> Self {
        Self::Runtime(e.to_string())
    }
}

impl From<SeqError> for CliError {
    fn from(e: SeqError) -> Self {
        Self::Seq(e)
    }
}

impl From<io::Error> for CliError {
    fn from(e: io::Error) -> Self {
        Self::Seq(SeqError::Io(e))
    }
}

impl From<StatsError> for CliError {
    fn from(e: StatsError) -> Self {
        Self::Runtime(e.to_string())
    }
}

type CliResult<T = ()> = Result<T, CliError>;

/// Whether the file starts with the gzip magic, which BAM always does.
fn is_bam(path: &Path) -> io::Result<bool> {
    let mut magic = [0u8; 2];
    let mut f = File::open(path)?;
    let n = f.read(&mut magic)?;
    Ok(n == 2 && magic == [0x1f, 0x8b])
}

fn default_index_path(bam: &Path) -> PathBuf {
    let mut p = bam.as_os_str().to_owned();
    p.push(".bai");
    PathBuf::from(p)
}

fn open_output(path: Option<&Path>) -> io::Result<Box<dyn Write>> {
    Ok(match path {
        Some(p) => Box::new(BufWriter::new(File::create(p)?)),
        None => Box::new(BufWriter::new(io::stdout().lock())),
    })
}

fn write_completions(shell: Shell) -> CliResult {
    let mut cmd = Cli::command();
    let mut out = io::stdout().lock();
    clap_complete::generate(shell, &mut cmd, "tptbio", &mut out);
    Ok(out.flush()?)
}

fn write_man_pages(dir: &Path) -> CliResult {
    fs::create_dir_all(dir)?;
    let mut written = 0;
    man_page(&Cli::command(), "tptbio", "tptbio", dir, &mut written)?;
    log_info(format_args!("wrote {written} man pages to {}", dir.display()));
    Ok(())
}

/// Writes the page for `cmd`, then recurses into its subcommands. `stem` names the
/// file (`tptbio-seq-view.1`) and `display` is the invocation (`tptbio seq view`).
fn man_page(
    cmd: &clap::Command,
    stem: &str,
    display: &str,
    dir: &Path,
    written: &mut usize,
) -> CliResult {
    // clap stores names as `&'static str`. The process exits after writing the pages,
    // so the few leaked strings do not accumulate.
    let page = cmd
        .clone()
        .name(&*Box::leak(stem.to_owned().into_boxed_str()))
        .bin_name(&*Box::leak(display.to_owned().into_boxed_str()));
    let mut buf = Vec::new();
    clap_mangen::Man::new(page).render(&mut buf)?;
    fs::write(dir.join(format!("{stem}.1")), buf)?;
    *written += 1;
    for sub in cmd.get_subcommands() {
        let name = sub.get_name();
        if name == "help" || sub.is_hide_set() {
            continue;
        }
        man_page(sub, &format!("{stem}-{name}"), &format!("{display} {name}"), dir, written)?;
    }
    Ok(())
}

fn seq_view(args: &ViewArgs, settings: &Settings) -> CliResult {
    match variants::detect(&args.input)? {
        kind @ (Kind::Vcf | Kind::VcfGz | Kind::Bcf) => {
            variants::view(kind, args, settings.output_format)
        }
        Kind::Bam => view_bam(args),
        Kind::Sam => {
            if args.region.is_some() {
                return Err(SeqError::Unsupported("region queries need a BAM input".into()).into());
            }
            view_sam(args)
        }
    }
}

fn view_sam(args: &ViewArgs) -> CliResult {
    let file = BufReader::with_capacity(1 << 16, File::open(&args.input)?);
    let mut reader = SamReader::new(file)?;
    let header = reader.header().clone();
    let mut out = open_output(args.output.as_deref())?;
    if args.header_only {
        out.write_all(header.text().as_bytes())?;
        return Ok(out.flush()?);
    }
    if args.bam {
        let mut w = BamWriter::new(out, header)?;
        while let Some(rec) = reader.next_record()? {
            w.write_record(&rec)?;
        }
        w.finish()?;
    } else {
        let mut w = SamWriter::new(out, header)?;
        while let Some(rec) = reader.next_record()? {
            w.write_record(&rec)?;
        }
        w.finish()?;
    }
    Ok(())
}

fn view_bam(args: &ViewArgs) -> CliResult {
    let file = BufReader::with_capacity(1 << 16, File::open(&args.input)?);
    let mut reader = BamReader::new(file)?;
    let header = reader.header().clone();
    let mut out = open_output(args.output.as_deref())?;
    if args.header_only {
        out.write_all(header.text().as_bytes())?;
        return Ok(out.flush()?);
    }

    let records = match &args.region {
        Some(text) => {
            let region = Region::parse(text)
                .map_err(|e| SeqError::Unsupported(format!("bad region {text}: {e}")))?;
            let index_path = args.index.clone().unwrap_or_else(|| default_index_path(&args.input));
            let index =
                Index::read_from(BufReader::new(File::open(&index_path).map_err(|e| {
                    SeqError::NotFound(format!("index {}: {e}", index_path.display()))
                })?))?;
            let (rid, beg, end) = resolve_region(&region, &header)?;
            Some(bai::fetch(&mut reader, &index, rid, beg, end)?)
        }
        None => None,
    };

    if args.bam {
        let mut w = BamWriter::new(out, header)?;
        write_selected(&mut reader, records, |rec| w.write_record(rec))?;
        w.finish()?;
    } else {
        let mut w = SamWriter::new(out, header)?;
        write_selected(&mut reader, records, |rec| w.write_record(rec))?;
        w.finish()?;
    }
    Ok(())
}

/// Writes either the pre-fetched region records or every record in the file.
fn write_selected<R: Read>(
    reader: &mut BamReader<R>,
    fetched: Option<Vec<tpt_bio_seq::Record>>,
    mut sink: impl FnMut(&tpt_bio_seq::Record) -> tpt_bio_seq::Result<()>,
) -> CliResult {
    match fetched {
        Some(records) => {
            for rec in &records {
                sink(rec)?;
            }
        }
        None => {
            while let Some(rec) = reader.next_record()? {
                sink(&rec)?;
            }
        }
    }
    Ok(())
}

/// Converts a 1-based inclusive region into `(ref id, 0-based beg, 0-based end)`.
fn resolve_region(region: &Region, header: &Header) -> CliResult<(usize, u64, u64)> {
    let rid = header
        .ref_id(region.contig())
        .ok_or_else(|| SeqError::NotFound(format!("reference {}", region.contig())))?;
    let len = header.references()[rid as usize].len;
    let (beg, end) = match region.interval() {
        Some(i) => (i.start(), i.end().min(len)),
        None => (0, len),
    };
    Ok((rid as usize, beg, end))
}

fn seq_index(input: &Path, output: Option<PathBuf>) -> CliResult {
    let mut reader = BamReader::new(BufReader::with_capacity(1 << 16, File::open(input)?))?;
    let index = bai::build_index(&mut reader)?;
    let out_path = output.unwrap_or_else(|| default_index_path(input));
    let mut w = BufWriter::new(File::create(&out_path)?);
    index.write_to(&mut w)?;
    w.flush()?;
    log_info(format_args!("wrote {}", out_path.display()));
    Ok(())
}

/// Reads every record from a SAM or BAM file.
fn read_all(path: &Path) -> CliResult<(Header, Vec<tpt_bio_seq::Record>)> {
    if is_bam(path)? {
        let mut reader = BamReader::new(BufReader::with_capacity(1 << 16, File::open(path)?))?;
        let header = reader.header().clone();
        let mut records = Vec::new();
        while let Some(rec) = reader.next_record()? {
            records.push(rec);
        }
        Ok((header, records))
    } else {
        let mut reader = SamReader::new(BufReader::with_capacity(1 << 16, File::open(path)?))?;
        let header = reader.header().clone();
        let mut records = Vec::new();
        while let Some(rec) = reader.next_record()? {
            records.push(rec);
        }
        Ok((header, records))
    }
}

fn seq_flagstat(input: &Path) -> CliResult {
    let (_, records) = read_all(input)?;
    let s = tpt_bio_seq::Flagstat::from_records(&records);
    let mut out = io::stdout().lock();
    let line = |out: &mut dyn Write, n: u64, label: &str| -> io::Result<()> {
        writeln!(out, "{n} + 0 {label}")
    };
    line(&mut out, s.total, "in total (QC-passed reads + QC-failed reads)")?;
    line(&mut out, s.secondary, "secondary")?;
    line(&mut out, s.supplementary, "supplementary")?;
    line(&mut out, s.duplicates, "duplicates")?;
    line(&mut out, s.mapped, "mapped")?;
    line(&mut out, s.paired, "paired in sequencing")?;
    line(&mut out, s.read1, "read1")?;
    line(&mut out, s.read2, "read2")?;
    line(&mut out, s.properly_paired, "properly paired")?;
    line(&mut out, s.both_mapped, "with itself and mate mapped")?;
    line(&mut out, s.singletons, "singletons")?;
    Ok(())
}

fn seq_sort(input: &Path, output: &Path) -> CliResult {
    let (header, mut records) = read_all(input)?;
    tpt_bio_seq::sort_coordinate(&mut records);
    // The output header must say the file is coordinate-sorted. Any existing SO
    // field on the @HD line is replaced, so the tag never appears twice.
    let text = header.text().to_string();
    let header = if header.hd().and_then(|hd| hd.get("SO")) == Some("coordinate") {
        header
    } else {
        let mut rewritten = String::new();
        let mut done = false;
        for line in text.lines() {
            if !done && line.starts_with("@HD") {
                let mut fields: Vec<&str> =
                    line.split('\t').filter(|f| !f.starts_with("SO:")).collect();
                fields.push("SO:coordinate");
                let _ = writeln!(rewritten, "{}", fields.join("\t"));
                done = true;
            } else {
                rewritten.push_str(line);
                rewritten.push('\n');
            }
        }
        if !done {
            rewritten = format!("@HD\tVN:1.6\tSO:coordinate\n{rewritten}");
        }
        Header::from_text_and_references(&rewritten, header.references().to_vec())?
    };
    let mut w = BamWriter::new(BufWriter::new(File::create(output)?), header)?;
    for rec in &records {
        w.write_record(rec)?;
    }
    w.finish()?;
    Ok(())
}

fn seq_export(input: &Path, format: Option<OutputFormat>, output: Option<&Path>) -> CliResult {
    let format = match format.unwrap_or(OutputFormat::Csv) {
        OutputFormat::Table => {
            log_warn("seq export has no table layout; writing CSV instead");
            OutputFormat::Csv
        }
        other => other,
    };
    let (header, records) = read_all(input)?;
    let text = match format {
        OutputFormat::Jsonl => tpt_bio_seq::records_to_jsonl(&records, &header),
        _ => tpt_bio_seq::records_to_csv(&records, &header)?,
    };
    let mut out = open_output(output)?;
    out.write_all(text.as_bytes())?;
    out.flush()?;
    Ok(())
}

/// A CSV table read into memory. Fields are kept as text until a command needs numbers.
struct Table {
    header: Vec<String>,
    rows: Vec<Vec<String>>,
}

impl Table {
    fn column_index(&self, name: &str) -> CliResult<usize> {
        self.header
            .iter()
            .position(|h| h == name)
            .ok_or_else(|| CliError::Runtime(format!("no column named {name:?}")))
    }

    /// Parses the listed columns of every row as finite numbers.
    fn numeric(&self, columns: &[usize]) -> CliResult<Matrix> {
        let mut data = Vec::with_capacity(self.rows.len() * columns.len());
        for (r, row) in self.rows.iter().enumerate() {
            for &c in columns {
                let field = &row[c];
                let value: f64 = field.parse().map_err(|_| {
                    CliError::Runtime(format!(
                        "row {}, column {:?}: {field:?} is not a number",
                        r + 1,
                        self.header[c]
                    ))
                })?;
                data.push(value);
            }
        }
        Ok(Matrix::new(self.rows.len(), columns.len(), data)?)
    }
}

/// Reads a comma-separated file with a header row. Blank lines are skipped.
/// Fields cannot contain commas or quotes.
fn read_table(path: &Path) -> CliResult<Table> {
    let reader = BufReader::with_capacity(1 << 16, File::open(path)?);
    let mut header: Option<Vec<String>> = None;
    let mut rows = Vec::new();
    for (lineno, line) in reader.lines().enumerate() {
        let line = line?;
        let line = line.trim_end_matches('\r');
        if line.trim().is_empty() {
            continue;
        }
        let fields: Vec<String> = line.split(',').map(|f| f.trim().to_owned()).collect();
        match &header {
            None => header = Some(fields),
            Some(h) => {
                if fields.len() != h.len() {
                    return Err(CliError::Runtime(format!(
                        "{}: line {}: expected {} fields, found {}",
                        path.display(),
                        lineno + 1,
                        h.len(),
                        fields.len()
                    )));
                }
                rows.push(fields);
            }
        }
    }
    let header = header
        .ok_or_else(|| CliError::Runtime(format!("{}: file has no header row", path.display())))?;
    Ok(Table { header, rows })
}

fn write_table(
    output: Option<&Path>,
    format: OutputFormat,
    header: &[String],
    rows: &[Vec<String>],
) -> CliResult {
    let mut out = open_output(output)?;
    output::write_rows(&mut *out, format, header, rows)?;
    out.flush()?;
    Ok(())
}

fn stats_normalize(args: &NormalizeArgs, settings: &Settings) -> CliResult {
    let table = read_table(&args.input)?;
    let group_col = args.group.as_deref().map(|g| table.column_index(g)).transpose()?;
    let feature_cols: Vec<usize> =
        (0..table.header.len()).filter(|&c| Some(c) != group_col).collect();
    let values = table.numeric(&feature_cols)?;
    let method = Method::from(args.method);
    let normalized = match group_col {
        Some(g) => {
            let labels: Vec<&str> = table.rows.iter().map(|r| r[g].as_str()).collect();
            normalize::normalize_by_group(&values, &labels, method)?
        }
        None => normalize::normalize_columns(&values, method)?,
    };

    let mut rows = table.rows.clone();
    for (r, row) in rows.iter_mut().enumerate() {
        for (k, &c) in feature_cols.iter().enumerate() {
            row[c] = normalized.get(r, k).to_string();
        }
    }
    log_info(format_args!(
        "normalized {} rows x {} features ({})",
        table.rows.len(),
        feature_cols.len(),
        if group_col.is_some() { "per group" } else { "all rows" }
    ));
    let format = settings.output_format.unwrap_or(OutputFormat::Csv);
    write_table(args.output.as_deref(), format, &table.header, &rows)
}

fn stats_pca(args: &PcaArgs, settings: &Settings) -> CliResult {
    let table = read_table(&args.input)?;
    let all: Vec<usize> = (0..table.header.len()).collect();
    let values = table.numeric(&all)?;
    let model = Pca::fit(&values, args.components)?;
    let scores = model.transform(&values)?;
    for (k, ratio) in model.explained_variance_ratio().iter().enumerate() {
        log_info(format_args!("PC{}: explained variance ratio {ratio:.4}", k + 1));
    }

    let header: Vec<String> = (1..=model.n_components()).map(|k| format!("PC{k}")).collect();
    let rows: Vec<Vec<String>> = (0..scores.rows())
        .map(|r| scores.row(r).iter().map(ToString::to_string).collect())
        .collect();
    let format = settings.output_format.unwrap_or(OutputFormat::Csv);
    write_table(args.output.as_deref(), format, &header, &rows)
}

/// Number of histogram bins for the Otsu threshold.
const OTSU_BINS: usize = 256;

/// Loads a PNG or TIFF image, chosen by file extension.
fn load_image(path: &Path) -> CliResult<Planes> {
    let ext = path.extension().and_then(OsStr::to_str).map(str::to_ascii_lowercase);
    let is_png = match ext.as_deref() {
        Some("png") => true,
        Some("tif" | "tiff") => false,
        _ => {
            return Err(CliError::Usage(format!(
                "{}: expected a .png, .tif or .tiff file",
                path.display()
            )));
        }
    };
    let source = FileSource::open(path)
        .map_err(|e| CliError::Runtime(format!("{}: {e}", path.display())))?;
    let decoded = if is_png { read_png(&source) } else { read_tiff(&source) };
    decoded.map_err(|e| CliError::Runtime(format!("{}: {e}", path.display())))
}

/// Takes one channel out of a decoded list of channels.
fn take_channel<T>(
    mut channels: Vec<Image<T>>,
    channel: usize,
    path: &Path,
) -> CliResult<Image<T>> {
    if channel >= channels.len() {
        return Err(CliError::Runtime(format!(
            "{}: channel {channel} is out of range, the image has {} channel(s)",
            path.display(),
            channels.len()
        )));
    }
    Ok(channels.swap_remove(channel))
}

/// Otsu-threshold one 2D channel, then label its foreground.
fn segment_channel<T: Copy + Into<f64>>(
    img: &Image<T>,
    connectivity: Connectivity2,
) -> CliResult<(Labeling, Option<f64>)> {
    if img.depth() != 1 {
        return Err(CliError::Runtime(format!(
            "segment needs a single 2D image, found {} pages",
            img.depth()
        )));
    }
    let threshold = otsu_threshold(img, OTSU_BINS)?;
    let mask = if let Some(t) = threshold {
        apply_threshold(img, t)
    } else {
        log_warn("channel is empty or constant: no foreground");
        Image::new_2d(img.width(), img.height(), vec![false; img.len()])?
    };
    Ok((label_2d(&mask, connectivity)?, threshold))
}

fn vision_segment(args: &SegmentArgs) -> CliResult {
    let ThresholdArg::Otsu = args.method;
    let connectivity = match args.connectivity {
        ConnectivityArg::Four => Connectivity2::Four,
        ConnectivityArg::Eight => Connectivity2::Eight,
    };
    let (labeling, threshold) = match load_image(&args.input)? {
        Planes::U8(v) => {
            segment_channel(&take_channel(v, args.channel, &args.input)?, connectivity)?
        }
        Planes::U16(v) => {
            segment_channel(&take_channel(v, args.channel, &args.input)?, connectivity)?
        }
        Planes::F32(v) => {
            segment_channel(&take_channel(v, args.channel, &args.input)?, connectivity)?
        }
    };
    let found = labeling.count;
    let filter =
        ObjectFilter { min_area: args.min_size, max_area: None, remove_border: args.remove_border };
    let kept = filter_objects(&labeling.labels, &filter);
    let kept_count = kept.count;
    let limit = u32::from(u16::MAX);
    if kept_count > limit {
        return Err(CliError::Runtime(format!(
            "{kept_count} objects do not fit a 16-bit label image (limit {limit})"
        )));
    }
    let out: Image<u16> = kept.labels.map(|&v| u16::try_from(v).unwrap_or(u16::MAX));
    save_tiff(&args.output, &out, None, false)
        .map_err(|e| CliError::Runtime(format!("{}: {e}", args.output.display())))?;
    let threshold = threshold.map_or_else(|| "none".to_owned(), |t| format!("{t:.3}"));
    println!(
        "{}: {found} objects found, {kept_count} kept (threshold {threshold})",
        args.input.display()
    );
    Ok(())
}

/// Converts channel `channel` of a decoded image to `u32` labels.
fn label_image(planes: Planes, channel: usize, path: &Path) -> CliResult<Image<u32>> {
    Ok(match planes {
        Planes::U8(v) => take_channel(v, channel, path)?.map(|&x| u32::from(x)),
        Planes::U16(v) => take_channel(v, channel, path)?.map(|&x| u32::from(x)),
        Planes::F32(v) => take_channel(v, channel, path)?
            .map(|&x| if x.is_finite() && x > 0.0 { x.round() as u32 } else { 0 }),
    })
}

/// Converts channel `channel` of a decoded image to `f64` samples.
fn intensity_image(planes: Planes, channel: usize, path: &Path) -> CliResult<Image<f64>> {
    Ok(match planes {
        Planes::U8(v) => take_channel(v, channel, path)?.map(|&x| f64::from(x)),
        Planes::U16(v) => take_channel(v, channel, path)?.map(|&x| f64::from(x)),
        Planes::F32(v) => take_channel(v, channel, path)?.map(|&x| f64::from(x)),
    })
}

fn vision_features(args: &FeaturesArgs, settings: &Settings) -> CliResult {
    let labels = label_image(load_image(&args.labels)?, 0, &args.labels)?;
    let shape = shape_features_2d(&labels)?;

    let intensity_path = args.intensity.as_ref().or(args.input.as_ref());
    let intensity = match intensity_path {
        Some(path) => {
            let img = intensity_image(load_image(path)?, args.channel, path)?;
            Some(intensity_features(&labels, &img)?)
        }
        None => None,
    };

    let mut header: Vec<String> = ["label", "area", "perimeter", "eccentricity", "solidity"]
        .iter()
        .map(|&s| s.to_owned())
        .collect();
    if intensity.is_some() {
        header.extend(
            [
                "intensity_mean",
                "intensity_std",
                "intensity_min",
                "intensity_max",
                "intensity_median",
                "intensity_integrated",
            ]
            .iter()
            .map(|&s| s.to_owned()),
        );
    }

    let rows: Vec<Vec<String>> = shape
        .iter()
        .enumerate()
        .map(|(i, s)| {
            let mut row = vec![
                s.label.to_string(),
                s.area.to_string(),
                s.perimeter.to_string(),
                s.eccentricity.to_string(),
                s.solidity.to_string(),
            ];
            if let Some(f) = intensity.as_ref().and_then(|list| list.get(i)) {
                row.extend([
                    f.mean.to_string(),
                    f.std.to_string(),
                    f.min.to_string(),
                    f.max.to_string(),
                    f.median.to_string(),
                    f.integrated.to_string(),
                ]);
            }
            row
        })
        .collect();
    let format = settings.output_format.unwrap_or(OutputFormat::Csv);
    write_table(Some(&args.output), format, &header, &rows)?;
    log_info(format_args!("{} objects written to {}", rows.len(), args.output.display()));
    Ok(())
}

fn pipeline_run(definition: &Path, settings: &Settings) -> CliResult {
    let text = fs::read_to_string(definition)
        .map_err(|e| CliError::Runtime(format!("{}: {e}", definition.display())))?;
    let pipeline = Pipeline::from_toml(&text)
        .map_err(|e| CliError::Usage(format!("{}: {e}", definition.display())))?;
    let outputs = pipeline.run().map_err(|e| CliError::Runtime(e.to_string()))?;
    let format = settings.output_format.unwrap_or(OutputFormat::Csv);
    let mut tables = 0usize;
    for value in &outputs {
        match value {
            Value::Table(m) => {
                tables += 1;
                let header: Vec<String> = (1..=m.cols()).map(|c| format!("c{c}")).collect();
                let rows: Vec<Vec<String>> = (0..m.rows())
                    .map(|r| m.row(r).iter().map(ToString::to_string).collect())
                    .collect();
                write_table(None, format, &header, &rows)?;
            }
            Value::Frame(frame) => {
                tables += 1;
                let rows: Vec<Vec<String>> = frame
                    .rows
                    .iter()
                    .map(|r| r.iter().map(ToString::to_string).collect())
                    .collect();
                write_table(None, format, &frame.columns, &rows)?;
            }
            Value::Labels(labels) => {
                let clusters = labels.iter().copied().max().map_or(0, |m| m + 1);
                println!("labels: {} rows, {clusters} clusters", labels.len());
            }
            other => println!("output: {}", value_kind(other)),
        }
    }
    log_info(format_args!(
        "{}: {} output(s), {tables} table(s) printed",
        definition.display(),
        outputs.len()
    ));
    Ok(())
}

/// A short name for the kind of a pipeline output, for summaries.
const fn value_kind(value: &Value) -> &'static str {
    match value {
        Value::Table(_) => "table",
        Value::Labels(_) => "labels",
        Value::Fastq(_) => "fastq reads",
        Value::Records(_) => "alignment records",
        Value::File(_) => "file",
        Value::Image(_) => "image",
        Value::Mask(_) => "mask",
        Value::Objects(_) => "objects",
        Value::Frame(_) => "frame",
    }
}

/// Per-channel `(min, max, mean)`, or `None` for an empty channel.
fn channel_stats<T: Copy + Into<f64>>(channels: &[Image<T>]) -> Vec<Option<(f64, f64, f64)>> {
    channels
        .iter()
        .map(|img| {
            if img.is_empty() {
                return None;
            }
            let (mut lo, mut hi, mut sum) = (f64::INFINITY, f64::NEG_INFINITY, 0.0);
            for &v in img.data() {
                let x: f64 = v.into();
                lo = lo.min(x);
                hi = hi.max(x);
                sum += x;
            }
            Some((lo, hi, sum / img.len() as f64))
        })
        .collect()
}

fn image_view(input: &Path) -> CliResult {
    let planes = load_image(input)?;
    let Some((width, height, depth)) = planes.shape() else {
        return Err(CliError::Runtime(format!("{}: image has no channels", input.display())));
    };
    let (pixel, stats) = match &planes {
        Planes::U8(v) => ("u8", channel_stats(v)),
        Planes::U16(v) => ("u16", channel_stats(v)),
        Planes::F32(v) => ("f32", channel_stats(v)),
    };
    println!("file: {}", input.display());
    println!("size: {width} x {height} x {depth} (width x height x depth)");
    println!("pixel type: {pixel}");
    println!("channels: {}", planes.channel_count());
    for (c, s) in stats.iter().enumerate() {
        match s {
            Some((lo, hi, mean)) => println!("channel {c}: min {lo}, max {hi}, mean {mean:.3}"),
            None => println!("channel {c}: empty"),
        }
    }
    Ok(())
}
