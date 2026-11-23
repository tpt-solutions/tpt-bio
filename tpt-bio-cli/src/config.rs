//! The `tptbio.toml` configuration file and the global settings it can provide.
//!
//! Lookup order: the file named by `--config`, otherwise `tptbio.toml` in the current
//! directory if it exists. Command-line flags always override the file, and the file
//! overrides the built-in defaults.

use std::fs;
use std::path::{Path, PathBuf};

use clap::ValueEnum;
use serde::Deserialize;

/// File name looked up in the current directory when `--config` is not given.
pub(crate) const DEFAULT_FILE: &str = "tptbio.toml";

/// Severity threshold for messages written to standard error. Errors are always shown.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Deserialize)]
#[serde(rename_all = "lowercase")]
pub(crate) enum LogLevel {
    Error,
    Warn,
    Info,
    Debug,
}

/// Layout of tabular results.
#[derive(Debug, Clone, Copy, PartialEq, Eq, ValueEnum, Deserialize)]
#[serde(rename_all = "lowercase")]
pub(crate) enum OutputFormat {
    /// Comma-separated values with a header row.
    Csv,
    /// One JSON object per line.
    Jsonl,
    /// Column-aligned plain text for reading in a terminal.
    Table,
}

/// Keys accepted in `tptbio.toml`. Unknown keys are an error.
#[derive(Debug, Default, Deserialize)]
#[serde(deny_unknown_fields)]
pub(crate) struct Config {
    /// Worker thread count. Recorded only: the CLI has no thread pool to configure.
    pub(crate) threads: Option<u32>,
    /// Default log level, used when no `-q` or `-v` flag is given.
    pub(crate) log_level: Option<LogLevel>,
    /// Default output format, used when `--output-format` is not given.
    pub(crate) output_format: Option<OutputFormat>,
}

/// Loads the configuration. Returns the file used, if any, and its contents.
///
/// An explicit path must exist. The default `tptbio.toml` is optional.
///
/// # Errors
/// Returns a message for an unreadable or invalid file.
pub(crate) fn load(explicit: Option<&Path>) -> Result<(Option<PathBuf>, Config), String> {
    let path = if let Some(p) = explicit {
        p.to_path_buf()
    } else {
        let default = PathBuf::from(DEFAULT_FILE);
        if !default.is_file() {
            return Ok((None, Config::default()));
        }
        default
    };
    let text = fs::read_to_string(&path)
        .map_err(|e| format!("cannot read config {}: {e}", path.display()))?;
    let config: Config =
        toml::from_str(&text).map_err(|e| format!("invalid config {}: {e}", path.display()))?;
    if config.threads == Some(0) {
        return Err(format!("config {}: threads must be at least 1", path.display()));
    }
    Ok((Some(path), config))
}
