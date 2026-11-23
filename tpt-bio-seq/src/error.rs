//! Error type for sequencing I/O.

use std::fmt;
use std::io;

/// Errors from reading or writing sequencing data.
#[derive(Debug)]
pub enum SeqError {
    /// An underlying I/O failure.
    Io(io::Error),
    /// The input is not well-formed. `offset` is a byte offset when known.
    Malformed {
        /// The format being parsed, for example `"BGZF"` or `"FASTQ"`.
        format: &'static str,
        /// Byte offset into the input, when known.
        offset: Option<u64>,
        /// Human-readable description of the problem.
        reason: String,
    },
    /// The data is well-formed but uses a feature this crate does not support.
    Unsupported(String),
    /// A requested index or reference is not present.
    NotFound(String),
}

impl SeqError {
    /// Builds a [`SeqError::Malformed`] without an offset.
    pub fn malformed(format: &'static str, reason: impl Into<String>) -> Self {
        Self::Malformed { format, offset: None, reason: reason.into() }
    }

    /// Builds a [`SeqError::Malformed`] at a byte offset.
    pub fn malformed_at(format: &'static str, offset: u64, reason: impl Into<String>) -> Self {
        Self::Malformed { format, offset: Some(offset), reason: reason.into() }
    }
}

impl fmt::Display for SeqError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Io(e) => write!(f, "I/O error: {e}"),
            Self::Malformed { format, offset: Some(o), reason } => {
                write!(f, "malformed {format} at byte {o}: {reason}")
            }
            Self::Malformed { format, offset: None, reason } => {
                write!(f, "malformed {format}: {reason}")
            }
            Self::Unsupported(what) => write!(f, "unsupported: {what}"),
            Self::NotFound(what) => write!(f, "not found: {what}"),
        }
    }
}

impl std::error::Error for SeqError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Io(e) => Some(e),
            _ => None,
        }
    }
}

impl From<io::Error> for SeqError {
    fn from(e: io::Error) -> Self {
        Self::Io(e)
    }
}

/// Result alias for this crate.
pub type Result<T, E = SeqError> = std::result::Result<T, E>;
