//! Memory-mapped file access.
//!
//! This is the only module in `tpt-bio-seq` that contains `unsafe` code. Mapping a
//! file is `unsafe` because another process may truncate or modify it while it is
//! mapped, which would be undefined behavior. The API is therefore documented as
//! requiring that the file is not modified while the map is alive.

use std::fs::File;
use std::io;
use std::ops::Deref;
use std::path::Path;

/// A read-only view of a file's bytes, memory-mapped when possible.
#[derive(Debug)]
pub struct MappedFile {
    inner: Inner,
}

#[derive(Debug)]
enum Inner {
    Mapped(memmap2::Mmap),
    Owned(Vec<u8>),
    Empty,
}

impl MappedFile {
    /// Maps a file read-only.
    ///
    /// # Errors
    /// Returns an I/O error if the file cannot be opened or mapped.
    ///
    /// # Safety
    /// The caller must not truncate or modify the file while the returned value
    /// is alive. Modification by another process while mapped is undefined behavior.
    #[allow(unsafe_code)]
    pub unsafe fn open(path: &Path) -> io::Result<Self> {
        let file = File::open(path)?;
        if file.metadata()?.len() == 0 {
            // Zero-length mappings are rejected on some platforms.
            return Ok(Self { inner: Inner::Empty });
        }
        // SAFETY: the file is opened read-only and the documented contract above
        // requires that it is not modified while mapped.
        let map = unsafe { memmap2::Mmap::map(&file)? };
        Ok(Self { inner: Inner::Mapped(map) })
    }

    /// Wraps bytes already in memory, without mapping anything. This is safe, and it lets
    /// code that takes a `MappedFile` be tested without a file on disk.
    #[must_use]
    pub fn from_vec(bytes: Vec<u8>) -> Self {
        if bytes.is_empty() {
            Self { inner: Inner::Empty }
        } else {
            Self { inner: Inner::Owned(bytes) }
        }
    }

    /// The file contents.
    #[must_use]
    pub fn as_slice(&self) -> &[u8] {
        match &self.inner {
            Inner::Mapped(m) => m,
            Inner::Owned(v) => v,
            Inner::Empty => &[],
        }
    }
}

impl Deref for MappedFile {
    type Target = [u8];

    fn deref(&self) -> &[u8] {
        self.as_slice()
    }
}

#[cfg(test)]
// The tests map files, which needs `unsafe`; the files are not modified while mapped.
#[allow(unsafe_code)]
mod tests {
    use super::*;
    use std::io::Write;

    #[test]
    fn maps_file_contents() {
        let dir = std::env::temp_dir().join(format!("tpt-bio-mmap-{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let path = dir.join("data.bin");
        {
            let mut f = File::create(&path).unwrap();
            f.write_all(b"hello mmap").unwrap();
        }
        // SAFETY: the test file is not modified while it is mapped.
        let m = unsafe { MappedFile::open(&path) }.unwrap();
        assert_eq!(&m[..], b"hello mmap");

        let empty = dir.join("empty.bin");
        File::create(&empty).unwrap();
        // SAFETY: the test file is not modified while it is mapped.
        assert!(unsafe { MappedFile::open(&empty) }.unwrap().is_empty());
        std::fs::remove_dir_all(&dir).unwrap();
    }
}
