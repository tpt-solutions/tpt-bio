//! `seq view` for VCF and BCF inputs, and the format detection shared with SAM and BAM.
//!
//! Detection checks the file contents first. The extension must agree with them when
//! it names a known format, so a mislabelled file is an error rather than a mis-parse.

use std::fs::File;
use std::io::{self, BufRead, BufReader, Read, Write};
use std::path::{Path, PathBuf};

use tpt_bio_core::Region;
use tpt_bio_seq::bcf::{self, BcfReader};
use tpt_bio_seq::tabix::{self, TabixIndex};
use tpt_bio_seq::vcf::VcfReader;
use tpt_bio_seq::{BgzfReader, SeqError};

use crate::config::OutputFormat;
use crate::output::write_rows;
use crate::{CliError, CliResult, ViewArgs, open_output};

/// A sequence or variant file format that `seq view` understands.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum Kind {
    Sam,
    Bam,
    Vcf,
    VcfGz,
    Bcf,
}

impl Kind {
    fn name(self) -> &'static str {
        match self {
            Self::Sam => "SAM",
            Self::Bam => "BAM",
            Self::Vcf => "VCF",
            Self::VcfGz => "bgzipped VCF",
            Self::Bcf => "BCF",
        }
    }
}

/// Decides the format from the contents, checked against the extension.
///
/// # Errors
/// Returns a runtime error if the extension and the contents disagree, or an I/O error.
pub(crate) fn detect(path: &Path) -> CliResult<Kind> {
    match (by_extension(path), sniff(path)?) {
        (Some(named), Some(found)) if named != found => Err(CliError::Runtime(format!(
            "{}: extension says {}, but the contents are {}",
            path.display(),
            named.name(),
            found.name()
        ))),
        (_, Some(found)) => Ok(found),
        (Some(named), None) => Ok(named),
        (None, None) => Ok(Kind::Sam),
    }
}

fn by_extension(path: &Path) -> Option<Kind> {
    let name = path.file_name()?.to_str()?.to_ascii_lowercase();
    let (rest, ext) = name.rsplit_once('.')?;
    match ext {
        "gz" => {
            let (_, inner) = rest.rsplit_once('.')?;
            (inner == "vcf").then_some(Kind::VcfGz)
        }
        "vcf" => Some(Kind::Vcf),
        "bcf" => Some(Kind::Bcf),
        "bam" => Some(Kind::Bam),
        "sam" => Some(Kind::Sam),
        _ => None,
    }
}

/// Reads the first bytes of the file. Gzip input is decompressed to look at the
/// payload, since BAM, BCF and bgzipped VCF all start with gzip magic.
fn sniff(path: &Path) -> CliResult<Option<Kind>> {
    let mut head = [0u8; 5];
    let n = read_head(&mut File::open(path)?, &mut head)?;
    if n == 0 {
        return Ok(None);
    }
    if n >= 2 && head[..2] == [0x1f, 0x8b] {
        let mut inner = BgzfReader::new(File::open(path)?);
        let mut payload = [0u8; 5];
        let m = read_head(&mut inner, &mut payload)?;
        let payload = &payload[..m];
        let kind = if payload.starts_with(b"BCF") {
            Kind::Bcf
        } else if payload.starts_with(b"#") {
            Kind::VcfGz
        } else {
            Kind::Bam
        };
        return Ok(Some(kind));
    }
    let kind = if head[..n].starts_with(b"BCF") {
        Kind::Bcf
    } else if head[0] == b'#' {
        Kind::Vcf
    } else {
        Kind::Sam
    };
    Ok(Some(kind))
}

fn read_head(r: &mut impl Read, buf: &mut [u8]) -> io::Result<usize> {
    let mut n = 0;
    while n < buf.len() {
        match r.read(&mut buf[n..])? {
            0 => break,
            k => n += k,
        }
    }
    Ok(n)
}

/// Prints the header or the records of a VCF or BCF file.
///
/// Without `format`, the output is VCF text. With a format, each record becomes a row
/// of the eight fixed VCF columns.
///
/// # Errors
/// Returns a runtime error for unreadable input, a bad region, or a missing index.
pub(crate) fn view(kind: Kind, args: &ViewArgs, format: Option<OutputFormat>) -> CliResult {
    if args.bam {
        return Err(SeqError::Unsupported(
            "-b writes BAM and does not apply to VCF or BCF input".into(),
        )
        .into());
    }
    let region = args.region.as_deref().map(parse_region).transpose()?;
    if region.is_none() && format.is_none() && !args.header_only {
        // Whole-file text output streams each record, so memory does not grow with the file.
        return stream_whole_file(kind, args);
    }
    let (header_text, lines) = read_variants(kind, args, region.as_ref())?;

    let mut out = open_output(args.output.as_deref())?;
    // Tabular output has its own header row, so the VCF header is for text output only.
    if args.header_only || format.is_none() {
        out.write_all(header_text.as_bytes())?;
    }
    if !args.header_only {
        match format {
            None => {
                for line in &lines {
                    writeln!(out, "{line}")?;
                }
            }
            Some(format) => {
                let (columns, rows) = variant_rows(&lines);
                write_rows(&mut *out, format, &columns, &rows)?;
            }
        }
    }
    out.flush()?;
    Ok(())
}

/// Returns the VCF header text and the record lines, which are empty when only the
/// header is wanted.
fn read_variants(
    kind: Kind,
    args: &ViewArgs,
    region: Option<&Region>,
) -> CliResult<(String, Vec<String>)> {
    let path = args.input.as_path();
    let header_only = args.header_only;
    match kind {
        Kind::Vcf => {
            if region.is_some() && !header_only {
                return Err(SeqError::Unsupported(
                    "region queries need a bgzipped VCF (.vcf.gz) with a tabix index".into(),
                )
                .into());
            }
            let mut reader = VcfReader::new(BufReader::with_capacity(1 << 16, File::open(path)?))?;
            let text = reader.header().text().to_owned();
            let lines = if header_only { Vec::new() } else { collect_vcf(&mut reader)? };
            Ok((text, lines))
        }
        Kind::VcfGz => {
            let mut reader = VcfReader::new_bgzf(File::open(path)?)?;
            let text = reader.header().text().to_owned();
            if header_only {
                return Ok((text, Vec::new()));
            }
            match region {
                None => Ok((text, collect_vcf(&mut reader)?)),
                Some(region) => {
                    let index = load_index(path, args.index.as_deref())?;
                    let mut bgzf = BgzfReader::new(File::open(path)?);
                    let (chrom, start, end) = region_bounds(region);
                    let lines = tabix::fetch_text(&mut bgzf, &index, chrom, start, end)?;
                    Ok((text, lines))
                }
            }
        }
        Kind::Bcf => {
            let mut reader = BcfReader::new_bgzf(File::open(path)?)?;
            let text = reader.header().vcf().text().to_owned();
            if header_only {
                return Ok((text, Vec::new()));
            }
            match region {
                None => {
                    let mut lines = Vec::new();
                    while let Some(rec) = reader.next_record()? {
                        lines.push(rec.line().to_owned());
                    }
                    Ok((text, lines))
                }
                Some(region) => {
                    let index = load_index(path, args.index.as_deref())?;
                    let mut bgzf = BgzfReader::new(File::open(path)?);
                    let (chrom, start, end) = region_bounds(region);
                    let lines =
                        bcf::fetch_lines(&mut bgzf, reader.header(), &index, chrom, start, end)?;
                    Ok((text, lines))
                }
            }
        }
        Kind::Sam | Kind::Bam => unreachable!("SAM and BAM are handled by the seq view path"),
    }
}

/// Writes a whole VCF or BCF file as VCF text, one record at a time.
fn stream_whole_file(kind: Kind, args: &ViewArgs) -> CliResult {
    let path = args.input.as_path();
    let mut out = open_output(args.output.as_deref())?;
    match kind {
        Kind::Vcf => {
            let mut reader = VcfReader::new(BufReader::with_capacity(1 << 16, File::open(path)?))?;
            out.write_all(reader.header().text().as_bytes())?;
            while let Some(rec) = reader.next_record()? {
                writeln!(out, "{}", rec.line())?;
            }
        }
        Kind::VcfGz => {
            let mut reader = VcfReader::new_bgzf(File::open(path)?)?;
            out.write_all(reader.header().text().as_bytes())?;
            while let Some(rec) = reader.next_record()? {
                writeln!(out, "{}", rec.line())?;
            }
        }
        Kind::Bcf => {
            let mut reader = BcfReader::new_bgzf(File::open(path)?)?;
            out.write_all(reader.header().vcf().text().as_bytes())?;
            while let Some(rec) = reader.next_record()? {
                writeln!(out, "{}", rec.line())?;
            }
        }
        Kind::Sam | Kind::Bam => unreachable!("SAM and BAM are handled by the seq view path"),
    }
    out.flush()?;
    Ok(())
}

fn collect_vcf<R: BufRead>(reader: &mut VcfReader<R>) -> CliResult<Vec<String>> {
    let mut lines = Vec::new();
    while let Some(rec) = reader.next_record()? {
        lines.push(rec.line().to_owned());
    }
    Ok(lines)
}

fn parse_region(text: &str) -> CliResult<Region> {
    Region::parse(text).map_err(|e| SeqError::Unsupported(format!("bad region {text}: {e}")).into())
}

/// The contig and the 1-based inclusive bounds that tabix queries use.
fn region_bounds(region: &Region) -> (&str, u64, u64) {
    match region.interval() {
        Some(i) => (region.contig(), i.start() + 1, i.end()),
        None => (region.contig(), 1, u64::from(u32::MAX)),
    }
}

/// Reads the index named by `--index`. By default it is `<input>.tbi` (the usual
/// index for VCF) or, if that is absent, `<input>.csi` (the usual index for BCF).
fn load_index(path: &Path, explicit: Option<&Path>) -> CliResult<TabixIndex> {
    let with_suffix = |suffix: &str| {
        let mut s = path.as_os_str().to_owned();
        s.push(suffix);
        PathBuf::from(s)
    };
    let index_path: PathBuf = if let Some(p) = explicit {
        p.to_path_buf()
    } else {
        let tbi = with_suffix(".tbi");
        if tbi.is_file() { tbi } else { with_suffix(".csi") }
    };
    let file = File::open(&index_path)
        .map_err(|e| SeqError::NotFound(format!("index {}: {e}", index_path.display())))?;
    Ok(TabixIndex::read_from(BufReader::new(file))?)
}

/// Splits VCF record lines into the eight fixed columns.
fn variant_rows(lines: &[String]) -> (Vec<String>, Vec<Vec<String>>) {
    let columns: Vec<String> = ["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO"]
        .iter()
        .map(|c| (*c).to_owned())
        .collect();
    let rows = lines
        .iter()
        .map(|line| {
            let mut fields: Vec<String> = line.split('\t').take(8).map(str::to_owned).collect();
            fields.resize(8, ".".to_owned());
            fields
        })
        .collect();
    (columns, rows)
}
