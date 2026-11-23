//! `tptbio`: the unified command-line tool.
//!
//! Implemented so far: `seq view` and `seq index`. The other command groups are
//! declared so the tree is stable, but they exit with status 2 until their phase lands.

use std::fmt::Write as _;
use std::fs::File;
use std::io::{self, BufReader, BufWriter, Read, Write};
use std::path::{Path, PathBuf};
use std::process::ExitCode;

use clap::{Parser, Subcommand};
use tpt_bio_core::Region;
use tpt_bio_seq::bai::{self, Index};
use tpt_bio_seq::{BamReader, BamWriter, Header, SamReader, SamWriter, SeqError};

#[derive(Debug, Parser)]
#[command(name = "tptbio", version, about = "High-performance bioinformatics and imaging toolkit")]
struct Cli {
    #[command(subcommand)]
    command: Command,
}

#[derive(Debug, Subcommand)]
enum Command {
    /// Sequencing data: SAM/BAM view and BAI indexing.
    Seq {
        #[command(subcommand)]
        cmd: SeqCmd,
    },
    /// Biological image analysis (not yet implemented).
    Vision {
        #[command(subcommand)]
        cmd: VisionCmd,
    },
    /// Workflow execution (not yet implemented).
    Pipeline {
        #[command(subcommand)]
        cmd: PipelineCmd,
    },
    /// Statistical analysis (not yet implemented).
    Stats,
    /// 3D viewer (not yet implemented).
    View,
}

#[derive(Debug, Subcommand)]
enum SeqCmd {
    /// Print SAM, or write BAM, optionally restricted to a region.
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
    /// Export records as CSV or JSON Lines. No Parquet: it needs Arrow.
    Export {
        /// SAM or BAM file.
        input: PathBuf,
        /// Output format.
        #[arg(short, long, value_enum, default_value_t = ExportFormat::Csv)]
        format: ExportFormat,
        /// Output file. Defaults to standard output.
        #[arg(short, long)]
        output: Option<PathBuf>,
    },
}

#[derive(Debug, Clone, Copy, clap::ValueEnum)]
enum ExportFormat {
    /// Comma-separated values, with a header row (via tpt-strata).
    Csv,
    /// One JSON object per line.
    Jsonl,
}

#[derive(Debug, clap::Args)]
struct ViewArgs {
    /// SAM or BAM file.
    input: PathBuf,
    /// Region such as `chr1:1000-2000` (1-based, inclusive). Requires an index.
    #[arg(short, long)]
    region: Option<String>,
    /// Index file to use with `--region`. Defaults to `<input>.bai`.
    #[arg(short, long)]
    index: Option<PathBuf>,
    /// Print only the header.
    #[arg(short = 'H', long)]
    header_only: bool,
    /// Write BAM instead of SAM text.
    #[arg(short = 'b', long)]
    bam: bool,
    /// Output file. Defaults to standard output.
    #[arg(short, long)]
    output: Option<PathBuf>,
}

#[derive(Debug, Subcommand)]
enum VisionCmd {
    /// Segment images (not yet implemented).
    Segment,
    /// Extract features (not yet implemented).
    Features,
}

#[derive(Debug, Subcommand)]
enum PipelineCmd {
    /// Run a pipeline definition (not yet implemented).
    Run,
}

fn main() -> ExitCode {
    let cli = Cli::parse();
    let result = match cli.command {
        Command::Seq { cmd: SeqCmd::View(args) } => seq_view(&args),
        Command::Seq { cmd: SeqCmd::Index { input, output } } => seq_index(&input, output),
        Command::Seq { cmd: SeqCmd::Flagstat { input } } => seq_flagstat(&input),
        Command::Seq { cmd: SeqCmd::Sort { input, output } } => seq_sort(&input, &output),
        Command::Seq { cmd: SeqCmd::Export { input, format, output } } => {
            seq_export(&input, format, output.as_deref())
        }
        Command::Vision { .. } => not_implemented("vision"),
        Command::Pipeline { .. } => not_implemented("pipeline"),
        Command::Stats => not_implemented("stats"),
        Command::View => not_implemented("view"),
    };
    match result {
        Ok(()) => ExitCode::SUCCESS,
        Err(CliError::NotImplemented(what)) => {
            eprintln!("tptbio: `{what}` is not implemented yet");
            ExitCode::from(2)
        }
        Err(CliError::Seq(e)) => {
            eprintln!("tptbio: {e}");
            ExitCode::FAILURE
        }
    }
}

enum CliError {
    NotImplemented(&'static str),
    Seq(SeqError),
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

type CliResult<T = ()> = Result<T, CliError>;

fn not_implemented(what: &'static str) -> CliResult {
    Err(CliError::NotImplemented(what))
}

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

fn seq_view(args: &ViewArgs) -> CliResult {
    if is_bam(&args.input)? {
        view_bam(args)
    } else {
        if args.region.is_some() {
            return Err(SeqError::Unsupported("region queries need a BAM input".into()).into());
        }
        view_sam(args)
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
    eprintln!("wrote {}", out_path.display());
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
    // The output header must say the file is coordinate-sorted.
    let text = header.text().to_string();
    let header = if text.contains("SO:coordinate") {
        header
    } else {
        let mut rewritten = String::new();
        let mut done = false;
        for line in text.lines() {
            if !done && line.starts_with("@HD") {
                let _ = writeln!(rewritten, "{line}\tSO:coordinate");
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

fn seq_export(input: &Path, format: ExportFormat, output: Option<&Path>) -> CliResult {
    let (header, records) = read_all(input)?;
    let text = match format {
        ExportFormat::Csv => tpt_bio_seq::records_to_csv(&records, &header)?,
        ExportFormat::Jsonl => tpt_bio_seq::records_to_jsonl(&records, &header),
    };
    let mut out = open_output(output)?;
    out.write_all(text.as_bytes())?;
    out.flush()?;
    Ok(())
}
