//! Image file input and output.
//!
//! Readers take a [`ByteSource`] so the same code runs on
//! local files, in-memory buffers and remote range reads. Decoded samples are
//! returned as [`Planes`]: one [`Image`] per channel, typed by sample type.
//!
//! Supported formats:
//! - [`png`]: 8-bit and 16-bit PNG, one plane per channel.
//! - [`tiff`]: striped and tiled TIFF and `BigTIFF`, 8-bit, 16-bit and float samples.
//!   Multi-page files are stacked into a volume with one page per `z` slice.
//! - [`ome`]: OME-TIFF with multiple series, read plane by `(z, c, t)`.
//! - [`mapped`]: memory-mapped access to uncompressed TIFF pages.
//! - [`zarr`]: Zarr v2 and v3 arrays, with the `bytes`, `gzip` and no-compression codecs.
//! - [`hcs`]: OME-NGFF high-content screening plate and well metadata.
//! - [`ndarray`]: conversion to and from `ndarray` arrays.
//! - [`mod@write`]: PNG and TIFF encoders, and Zarr v2 writing.

use std::fmt;
use std::io::{self, Read, Seek, SeekFrom};

use crate::image::{Image, ImageError};
use crate::source::ByteSource;

pub mod hcs;
pub mod mapped;
pub mod ndarray;
pub mod ome;
pub mod png;
pub mod sample;
pub mod tiff;
pub mod write;
pub mod zarr;

pub use sample::Sample;

/// Errors from image and container formats.
#[derive(Debug)]
pub enum FormatError {
    /// An I/O error from the underlying source or file.
    Io(io::Error),
    /// An error from the image container.
    Image(ImageError),
    /// The data could not be decoded.
    Decode(String),
    /// The data is well formed but uses a feature this crate does not implement.
    Unsupported(String),
    /// The data is inconsistent with what the format requires.
    Invalid(String),
}

/// Result type for image and container formats.
pub type FormatResult<T> = Result<T, FormatError>;

impl fmt::Display for FormatError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Io(e) => write!(f, "I/O error: {e}"),
            Self::Image(e) => write!(f, "image error: {e}"),
            Self::Decode(m) => write!(f, "decode error: {m}"),
            Self::Unsupported(m) => write!(f, "unsupported: {m}"),
            Self::Invalid(m) => write!(f, "invalid data: {m}"),
        }
    }
}

impl std::error::Error for FormatError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Io(e) => Some(e),
            Self::Image(e) => Some(e),
            _ => None,
        }
    }
}

impl From<io::Error> for FormatError {
    fn from(e: io::Error) -> Self {
        Self::Io(e)
    }
}

impl From<ImageError> for FormatError {
    fn from(e: ImageError) -> Self {
        Self::Image(e)
    }
}

/// Decoded samples, one [`Image`] per channel. Each channel is 2D (`depth == 1`)
/// for a single PNG or TIFF page, or a volume with one page per `z` slice.
#[derive(Debug, Clone, PartialEq)]
pub enum Planes {
    /// 8-bit unsigned samples.
    U8(Vec<Image<u8>>),
    /// 16-bit unsigned samples.
    U16(Vec<Image<u16>>),
    /// 32-bit float samples.
    F32(Vec<Image<f32>>),
}

impl Planes {
    /// Number of channels.
    #[must_use]
    pub fn channel_count(&self) -> usize {
        match self {
            Self::U8(v) => v.len(),
            Self::U16(v) => v.len(),
            Self::F32(v) => v.len(),
        }
    }

    /// Shape `(width, height, depth)` shared by all channels, or `None` if there
    /// are no channels.
    #[must_use]
    pub fn shape(&self) -> Option<(usize, usize, usize)> {
        match self {
            Self::U8(v) => v.first().map(Image::shape),
            Self::U16(v) => v.first().map(Image::shape),
            Self::F32(v) => v.first().map(Image::shape),
        }
    }
}

/// Adapts a [`ByteSource`] to [`Read`] and [`Seek`], so container parsers that
/// need a stream (such as the `tiff` crate's decoder) can use any source.
#[derive(Debug)]
pub struct SourceReader<'a, S: ByteSource + ?Sized> {
    source: &'a S,
    pos: u64,
}

impl<'a, S: ByteSource + ?Sized> SourceReader<'a, S> {
    /// Creates a reader positioned at offset 0.
    #[must_use]
    pub fn new(source: &'a S) -> Self {
        Self { source, pos: 0 }
    }
}

impl<S: ByteSource + ?Sized> Read for SourceReader<'_, S> {
    fn read(&mut self, buf: &mut [u8]) -> io::Result<usize> {
        let len = ByteSource::len(self.source);
        if buf.is_empty() || self.pos >= len {
            return Ok(0);
        }
        let remaining = len - self.pos;
        let n = usize::try_from(remaining.min(buf.len() as u64)).unwrap_or(buf.len());
        self.source.read_exact_at(self.pos, &mut buf[..n])?;
        self.pos += n as u64;
        Ok(n)
    }
}

impl<S: ByteSource + ?Sized> Seek for SourceReader<'_, S> {
    fn seek(&mut self, from: SeekFrom) -> io::Result<u64> {
        let len = i128::from(ByteSource::len(self.source));
        let target: i128 = match from {
            SeekFrom::Start(n) => i128::from(n),
            SeekFrom::End(d) => len + i128::from(d),
            SeekFrom::Current(d) => i128::from(self.pos) + i128::from(d),
        };
        let target = u64::try_from(target)
            .map_err(|_| io::Error::new(io::ErrorKind::InvalidInput, "seek before start"))?;
        self.pos = target;
        Ok(target)
    }
}

/// Reads an entire source into memory.
///
/// # Errors
/// Returns an I/O error if the source cannot be read.
pub fn read_all<S: ByteSource + ?Sized>(source: &S) -> FormatResult<Vec<u8>> {
    let len = usize::try_from(ByteSource::len(source))
        .map_err(|_| FormatError::Unsupported("source too large for memory".into()))?;
    Ok(source.read_vec_at(0, len)?)
}

/// Builds a [`FormatError::Decode`] from any displayable error.
pub(crate) fn decode_err(e: impl fmt::Display) -> FormatError {
    FormatError::Decode(e.to_string())
}

/// Splits interleaved samples (`samples` per pixel, x fastest) into one image per channel.
pub(crate) fn split_channels<T: Copy>(
    data: &[T],
    width: usize,
    height: usize,
    samples: usize,
) -> FormatResult<Vec<Image<T>>> {
    let expected = width * height * samples;
    if data.len() != expected {
        return Err(FormatError::Invalid(format!(
            "decoded {} samples, expected {expected} for {width}x{height}x{samples}",
            data.len()
        )));
    }
    (0..samples)
        .map(|c| {
            let plane: Vec<T> = data.iter().skip(c).step_by(samples).copied().collect();
            Ok(Image::new_2d(width, height, plane)?)
        })
        .collect()
}

impl From<::tiff::TiffError> for FormatError {
    fn from(e: ::tiff::TiffError) -> Self {
        Self::Decode(format!("TIFF: {e}"))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn source_reader_reads_and_seeks() {
        let data: Vec<u8> = (0u8..20).collect();
        let mut r = SourceReader::new(&data);
        let mut buf = [0u8; 4];
        r.read_exact(&mut buf).unwrap();
        assert_eq!(buf, [0, 1, 2, 3]);
        r.seek(SeekFrom::End(-2)).unwrap();
        let mut tail = Vec::new();
        r.read_to_end(&mut tail).unwrap();
        assert_eq!(tail, [18, 19]);
        assert!(r.seek(SeekFrom::Current(-100)).is_err());
    }

    #[test]
    fn planes_report_shape() {
        let p = Planes::U8(vec![Image::new_2d(2, 3, vec![0u8; 6]).unwrap()]);
        assert_eq!(p.channel_count(), 1);
        assert_eq!(p.shape(), Some((2, 3, 1)));
    }
}
