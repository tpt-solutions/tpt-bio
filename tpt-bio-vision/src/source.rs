//! Random-access byte sources for image and container formats.
//!
//! A [`ByteSource`] is a read-only, sized byte store that supports reads at an
//! arbitrary offset. Format readers (TIFF, Zarr chunks, OME-XML) take a source
//! instead of a file path, so the same code runs on local files, in-memory
//! buffers and, later, HTTP range reads.
//!
//! Memory-mapped access is provided by `tpt-bio-seq`, which owns the crate's only
//! `unsafe` wrapper; this crate forbids `unsafe`.

use std::fs::File;
use std::io::{self, Read, Seek, SeekFrom};
use std::path::Path;
use std::sync::Mutex;

/// A sized, read-only byte store with reads at arbitrary offsets.
pub trait ByteSource {
    /// Total length in bytes.
    fn len(&self) -> u64;

    /// Whether the source is empty.
    fn is_empty(&self) -> bool {
        self.len() == 0
    }

    /// Fills `buf` from `offset`. Returns [`io::ErrorKind::UnexpectedEof`] if the
    /// range extends past the end of the source.
    ///
    /// # Errors
    /// Returns an I/O error if the range is out of bounds or the read fails.
    fn read_exact_at(&self, offset: u64, buf: &mut [u8]) -> io::Result<()>;

    /// Reads `len` bytes starting at `offset` into a new vector.
    ///
    /// # Errors
    /// Same as [`ByteSource::read_exact_at`].
    fn read_vec_at(&self, offset: u64, len: usize) -> io::Result<Vec<u8>> {
        let mut buf = vec![0u8; len];
        self.read_exact_at(offset, &mut buf)?;
        Ok(buf)
    }
}

/// Checks that `offset + len` lies within a source of length `total`.
fn check_range(offset: u64, len: usize, total: u64) -> io::Result<()> {
    let end = offset
        .checked_add(len as u64)
        .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidInput, "range overflows u64"))?;
    if end > total {
        return Err(io::Error::new(
            io::ErrorKind::UnexpectedEof,
            format!("read of {len} bytes at {offset} exceeds source length {total}"),
        ));
    }
    Ok(())
}

impl ByteSource for [u8] {
    fn len(&self) -> u64 {
        <[u8]>::len(self) as u64
    }

    fn read_exact_at(&self, offset: u64, buf: &mut [u8]) -> io::Result<()> {
        check_range(offset, buf.len(), ByteSource::len(self))?;
        // Bounds were checked above, so the conversion and slice are in range.
        let start = offset as usize;
        buf.copy_from_slice(&self[start..start + buf.len()]);
        Ok(())
    }
}

impl ByteSource for Vec<u8> {
    fn len(&self) -> u64 {
        self.as_slice().len() as u64
    }

    fn read_exact_at(&self, offset: u64, buf: &mut [u8]) -> io::Result<()> {
        self.as_slice().read_exact_at(offset, buf)
    }
}

/// A local file read through a seek and a lock.
///
/// Reads take a mutex because seeking needs `&mut File`. The file length is
/// captured at open time; the file must not be truncated while the source is alive.
#[derive(Debug)]
pub struct FileSource {
    file: Mutex<File>,
    len: u64,
}

impl FileSource {
    /// Opens a file for reading.
    ///
    /// # Errors
    /// Returns an I/O error if the file cannot be opened or its metadata read.
    pub fn open(path: &Path) -> io::Result<Self> {
        let file = File::open(path)?;
        let len = file.metadata()?.len();
        Ok(Self { file: Mutex::new(file), len })
    }
}

impl ByteSource for FileSource {
    fn len(&self) -> u64 {
        self.len
    }

    fn read_exact_at(&self, offset: u64, buf: &mut [u8]) -> io::Result<()> {
        check_range(offset, buf.len(), self.len)?;
        let mut file = self
            .file
            .lock()
            .map_err(|_| io::Error::other("file source lock poisoned"))?;
        file.seek(SeekFrom::Start(offset))?;
        file.read_exact(buf)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;

    #[test]
    fn slice_reads_within_bounds() {
        let data: &[u8] = b"0123456789";
        assert_eq!(ByteSource::len(data), 10);
        assert_eq!(data.read_vec_at(2, 3).unwrap(), b"234");
        assert!(data.read_vec_at(8, 3).is_err());
        assert!(data.read_vec_at(u64::MAX, 2).is_err());
    }

    #[test]
    fn vec_source_matches_slice() {
        let v = vec![1u8, 2, 3, 4];
        assert_eq!(v.read_vec_at(1, 2).unwrap(), [2, 3]);
        assert!(!v.is_empty());
        assert!(Vec::<u8>::new().is_empty());
    }

    #[test]
    fn file_source_reads_at_offsets() {
        let dir = std::env::temp_dir().join(format!("tpt-bio-source-{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let path = dir.join("bytes.bin");
        {
            let mut f = File::create(&path).unwrap();
            f.write_all(b"abcdefghij").unwrap();
        }
        let src = FileSource::open(&path).unwrap();
        assert_eq!(src.len(), 10);
        assert_eq!(src.read_vec_at(7, 3).unwrap(), b"hij");
        assert_eq!(src.read_vec_at(0, 2).unwrap(), b"ab");
        assert!(src.read_vec_at(9, 2).is_err());
        let _ = std::fs::remove_dir_all(&dir);
    }
}
