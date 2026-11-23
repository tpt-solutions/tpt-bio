//! Reference sequences for reference-based CRAM decoding and encoding.
//!
//! A [`CramReference`] holds a FASTA file in memory together with its `.fai`
//! index (built on load). Bases are returned upper-cased, which is how CRAM
//! compares and checksums reference regions.

use std::path::Path;

use crate::cram::md5;
use crate::error::Result;
use crate::fasta::Fai;

/// A FASTA reference held in memory with its index.
#[derive(Debug, Clone)]
pub struct CramReference {
    data: Vec<u8>,
    fai: Fai,
}

impl CramReference {
    /// Loads a reference from in-memory FASTA bytes.
    ///
    /// # Errors
    /// Returns an error if the bytes are not FASTA that the index builder accepts.
    pub fn from_fasta_bytes(data: Vec<u8>) -> Result<Self> {
        let fai = Fai::build(&data)?;
        Ok(Self { data, fai })
    }

    /// Loads a reference from a FASTA file on disk.
    ///
    /// # Errors
    /// Returns [`crate::SeqError::Io`] if the file cannot be read, or the errors of
    /// [`CramReference::from_fasta_bytes`].
    pub fn from_path(path: impl AsRef<Path>) -> Result<Self> {
        Self::from_fasta_bytes(std::fs::read(path)?)
    }

    /// Whether the reference has a sequence with this name.
    #[must_use]
    pub fn contains(&self, name: &str) -> bool {
        self.fai.get(name).is_some()
    }

    /// Upper-cased bases `[beg, end)` (0-based) of the named sequence. The end is
    /// clamped to the sequence length.
    ///
    /// # Errors
    /// Returns [`crate::SeqError::NotFound`] for an unknown name.
    pub fn fetch(&self, name: &str, beg: u64, end: u64) -> Result<Vec<u8>> {
        let mut bases = self.fai.fetch(&self.data, name, beg, end)?;
        bases.make_ascii_uppercase();
        Ok(bases)
    }

    /// MD5 of the upper-cased bases `[beg, end)` of the named sequence.
    ///
    /// # Errors
    /// As [`CramReference::fetch`].
    pub fn region_md5(&self, name: &str, beg: u64, end: u64) -> Result<[u8; 16]> {
        Ok(md5::digest(&self.fetch(name, beg, end)?))
    }
}
