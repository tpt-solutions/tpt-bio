//! Memory-mapped access to TIFF pages.
//!
//! A [`MappedTiff`] maps a TIFF file with [`tpt_bio_seq::mmap::MappedFile`] and
//! records, for each page, whether its pixel data is one contiguous, uncompressed
//! byte range. For such pages [`MappedTiff::plane_bytes`] returns a slice of the
//! mapping with no copy. Other pages (compressed, tiled, or split into
//! non-adjacent strips) are decoded on request through the `tiff` crate.
//!
//! The mapping is read-only, but the operating system does not stop another
//! process from changing the file. **The caller must not truncate or modify the
//! file while a `MappedTiff` is alive.** Modifying it during that time is
//! undefined behaviour, so only map files that are not being written.

use std::ops::Deref;

use ::tiff::ColorType;
use ::tiff::decoder::Decoder;
use ::tiff::tags::Tag;
use tpt_bio_seq::mmap::MappedFile;

use super::tiff::read_page;
use super::{FormatError, FormatResult, Planes, SourceReader};
use crate::image::Image;

/// Sample type of a TIFF page.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SampleKind {
    /// 8-bit unsigned.
    U8,
    /// 16-bit unsigned.
    U16,
    /// 32-bit float.
    F32,
}

/// Layout of one TIFF page within the mapped file.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct MappedPage {
    /// Width in pixels.
    pub width: usize,
    /// Height in pixels.
    pub height: usize,
    /// Samples per pixel. Interleaved when greater than 1.
    pub samples: usize,
    /// Sample type.
    pub kind: SampleKind,
    /// Whether the file stores multi-byte samples little-endian.
    pub little_endian: bool,
    /// Byte range `(start, end)` of the pixel data in the file, when it is one
    /// contiguous uncompressed run. `None` otherwise.
    pub byte_range: Option<(usize, usize)>,
}

/// A memory-mapped TIFF file with per-page layout information.
///
/// See the [module documentation](self) for the mapping contract.
#[derive(Debug)]
pub struct MappedTiff<B: Deref<Target = [u8]> = MappedFile> {
    data: B,
    pages: Vec<MappedPage>,
}

impl<B: Deref<Target = [u8]>> MappedTiff<B> {
    /// Scans the pages of a TIFF held in `data`.
    ///
    /// Pass a [`MappedFile`] to use a memory map, or a `Vec<u8>` for an in-memory
    /// file. Mapping a file requires `unsafe` in `tpt-bio-seq`, so this crate takes
    /// the mapped bytes rather than a path.
    ///
    /// # Errors
    /// Returns [`FormatError::Unsupported`] or [`FormatError::Decode`] for pages the
    /// reader cannot describe.
    pub fn from_backing(data: B) -> FormatResult<Self> {
        let pages = scan_pages(&data)?;
        Ok(Self { data, pages })
    }

    /// Layout of every page, in IFD order.
    #[must_use]
    pub fn pages(&self) -> &[MappedPage] {
        &self.pages
    }

    /// Number of pages.
    #[must_use]
    pub fn page_count(&self) -> usize {
        self.pages.len()
    }

    /// The raw pixel bytes of a page as a slice of the mapping, without copying.
    ///
    /// The bytes are in the file's byte order and interleaved if the page has
    /// several samples per pixel.
    ///
    /// # Errors
    /// Returns [`FormatError::Invalid`] for an unknown page, and
    /// [`FormatError::Unsupported`] if the page is not one contiguous uncompressed
    /// run (use [`MappedTiff::plane_u8`], [`MappedTiff::plane_u16`] or
    /// [`MappedTiff::plane_f32`] to decode those).
    pub fn plane_bytes(&self, page: usize) -> FormatResult<&[u8]> {
        let p = self.page(page)?;
        let (start, end) = p.byte_range.ok_or_else(|| {
            FormatError::Unsupported(format!("page {page} is not contiguous and uncompressed"))
        })?;
        Ok(&self.data[start..end])
    }

    /// Decodes a page as 8-bit samples. Copies only the page's own bytes.
    ///
    /// # Errors
    /// Returns [`FormatError::Invalid`] for a page of another sample type or with
    /// several samples per pixel.
    pub fn plane_u8(&self, page: usize) -> FormatResult<Image<u8>> {
        let p = self.single_sample_page(page, SampleKind::U8)?;
        match p.byte_range {
            Some(_) => Ok(Image::new_2d(p.width, p.height, self.plane_bytes(page)?.to_vec())?),
            None => take_first(read_page(&self.data[..], page)?, |planes| match planes {
                Planes::U8(v) => Some(v),
                _ => None,
            }),
        }
    }

    /// Decodes a page as 16-bit samples, converting from the file's byte order.
    ///
    /// # Errors
    /// Returns [`FormatError::Invalid`] for a page of another sample type or with
    /// several samples per pixel.
    pub fn plane_u16(&self, page: usize) -> FormatResult<Image<u16>> {
        let p = self.single_sample_page(page, SampleKind::U16)?;
        match p.byte_range {
            Some(_) => {
                let little = p.little_endian;
                let data: Vec<u16> = self
                    .plane_bytes(page)?
                    .chunks_exact(2)
                    .map(|c| {
                        let pair = [c[0], c[1]];
                        if little { u16::from_le_bytes(pair) } else { u16::from_be_bytes(pair) }
                    })
                    .collect();
                Ok(Image::new_2d(p.width, p.height, data)?)
            }
            None => take_first(read_page(&self.data[..], page)?, |planes| match planes {
                Planes::U16(v) => Some(v),
                _ => None,
            }),
        }
    }

    /// Decodes a page as 32-bit float samples, converting from the file's byte order.
    ///
    /// # Errors
    /// Returns [`FormatError::Invalid`] for a page of another sample type or with
    /// several samples per pixel.
    pub fn plane_f32(&self, page: usize) -> FormatResult<Image<f32>> {
        let p = self.single_sample_page(page, SampleKind::F32)?;
        match p.byte_range {
            Some(_) => {
                let little = p.little_endian;
                let data: Vec<f32> = self
                    .plane_bytes(page)?
                    .chunks_exact(4)
                    .map(|c| {
                        let quad = [c[0], c[1], c[2], c[3]];
                        if little { f32::from_le_bytes(quad) } else { f32::from_be_bytes(quad) }
                    })
                    .collect();
                Ok(Image::new_2d(p.width, p.height, data)?)
            }
            None => take_first(read_page(&self.data[..], page)?, |planes| match planes {
                Planes::F32(v) => Some(v),
                _ => None,
            }),
        }
    }

    fn page(&self, index: usize) -> FormatResult<&MappedPage> {
        self.pages
            .get(index)
            .ok_or_else(|| FormatError::Invalid(format!("page {index} out of range")))
    }

    fn single_sample_page(&self, index: usize, kind: SampleKind) -> FormatResult<&MappedPage> {
        let p = self.page(index)?;
        if p.kind != kind {
            return Err(FormatError::Invalid(format!(
                "page {index} is {:?}, not {kind:?}",
                p.kind
            )));
        }
        if p.samples != 1 {
            return Err(FormatError::Invalid(format!("page {index} has {} samples", p.samples)));
        }
        Ok(p)
    }
}

fn take_first<T: Copy>(
    planes: Planes,
    pick: fn(Planes) -> Option<Vec<Image<T>>>,
) -> FormatResult<Image<T>> {
    pick(planes)
        .and_then(|mut v| if v.is_empty() { None } else { Some(v.remove(0)) })
        .ok_or_else(|| FormatError::Invalid("page sample type changed during decode".into()))
}

fn bits_of(ct: ColorType) -> u8 {
    match ct {
        ColorType::Gray(b)
        | ColorType::RGB(b)
        | ColorType::Palette(b)
        | ColorType::GrayA(b)
        | ColorType::RGBA(b)
        | ColorType::CMYK(b)
        | ColorType::YCbCr(b) => b,
    }
}

fn scan_pages(bytes: &[u8]) -> FormatResult<Vec<MappedPage>> {
    let little_endian = if bytes.starts_with(b"II") {
        true
    } else if bytes.starts_with(b"MM") {
        false
    } else {
        return Err(FormatError::Invalid("not a TIFF file".into()));
    };
    let mut decoder = Decoder::new(SourceReader::new(bytes))?;
    let mut pages = Vec::new();
    loop {
        pages.push(page_layout(&mut decoder, bytes, little_endian)?);
        if !decoder.more_images() {
            break;
        }
        decoder.next_image()?;
    }
    Ok(pages)
}

fn page_layout<R: std::io::Read + std::io::Seek>(
    decoder: &mut Decoder<R>,
    bytes: &[u8],
    little_endian: bool,
) -> FormatResult<MappedPage> {
    let (w, h) = decoder.dimensions()?;
    let (width, height) = (w as usize, h as usize);
    let samples = decoder.find_tag_unsigned::<u32>(Tag::SamplesPerPixel)?.unwrap_or(1) as usize;
    let format = decoder.find_tag_unsigned::<u16>(Tag::SampleFormat)?.unwrap_or(1);
    let bits = bits_of(decoder.colortype()?);
    let kind = match (bits, format) {
        (8, 1) => SampleKind::U8,
        (16, 1) => SampleKind::U16,
        (32, 3) => SampleKind::F32,
        _ => {
            return Err(FormatError::Unsupported(format!(
                "TIFF page with {bits}-bit samples of format {format}"
            )));
        }
    };
    let bytes_per_sample = usize::from(bits / 8);
    let expected = width * height * samples * bytes_per_sample;

    let compression = decoder.find_tag_unsigned::<u16>(Tag::Compression)?.unwrap_or(1);
    let tiled = decoder.find_tag(Tag::TileOffsets)?.is_some();
    let planar = decoder.find_tag_unsigned::<u16>(Tag::PlanarConfiguration)?.unwrap_or(1);
    let mut byte_range = None;
    if compression == 1 && !tiled && (planar == 1 || samples == 1) {
        let offsets = decoder.find_tag_unsigned_vec::<u64>(Tag::StripOffsets)?;
        let counts = decoder.find_tag_unsigned_vec::<u64>(Tag::StripByteCounts)?;
        if let (Some(offsets), Some(counts)) = (offsets, counts) {
            byte_range = contiguous_range(&offsets, &counts, expected, bytes.len());
        }
    }
    Ok(MappedPage { width, height, samples, kind, little_endian, byte_range })
}

/// The byte range covered by strips that abut one another, if their total size is `expected`.
fn contiguous_range(
    offsets: &[u64],
    counts: &[u64],
    expected: usize,
    file_len: usize,
) -> Option<(usize, usize)> {
    if offsets.is_empty() || offsets.len() != counts.len() {
        return None;
    }
    let total: u64 = counts.iter().sum();
    if total != expected as u64 {
        return None;
    }
    for i in 1..offsets.len() {
        if offsets[i] != offsets[i - 1].checked_add(counts[i - 1])? {
            return None;
        }
    }
    let start = usize::try_from(offsets[0]).ok()?;
    let end = start.checked_add(expected)?;
    (end <= file_len).then_some((start, end))
}

#[cfg(test)]
mod tests {
    use super::*;
    use ::tiff::encoder::{TiffEncoder, colortype};
    use std::fs;
    use std::io::{Cursor, Write};
    use std::path::Path;

    fn temp_dir(tag: &str) -> std::path::PathBuf {
        let dir = std::env::temp_dir().join(format!("tpt-bio-mapped-{tag}-{}", std::process::id()));
        fs::create_dir_all(&dir).unwrap();
        dir
    }

    fn write_u16_stack(path: &Path, pages: &[Vec<u16>]) {
        let mut out = Cursor::new(Vec::new());
        {
            let mut enc = TiffEncoder::new(&mut out).unwrap();
            for p in pages {
                enc.write_image::<colortype::Gray16>(2, 2, p).unwrap();
            }
        }
        fs::File::create(path).unwrap().write_all(&out.into_inner()).unwrap();
    }

    #[test]
    fn uncompressed_pages_are_zero_copy_slices() {
        let dir = temp_dir("zero-copy");
        let path = dir.join("stack.tif");
        write_u16_stack(&path, &[vec![1, 2, 3, 4], vec![500, 600, 700, 800]]);

        let mapped = MappedTiff::from_backing(fs::read(&path).unwrap()).unwrap();
        assert_eq!(mapped.page_count(), 2);
        assert!(mapped.pages()[0].byte_range.is_some());
        let bytes = mapped.plane_bytes(1).unwrap();
        assert_eq!(bytes.len(), 8);
        let decoded: Vec<u16> =
            bytes.chunks_exact(2).map(|c| u16::from_le_bytes([c[0], c[1]])).collect();
        assert_eq!(decoded, [500, 600, 700, 800]);

        let img = mapped.plane_u16(0).unwrap();
        assert_eq!(img.data(), &[1, 2, 3, 4]);
        assert!(mapped.plane_u8(0).is_err());
        assert!(mapped.plane_bytes(2).is_err());
        drop(mapped);
        fs::remove_dir_all(&dir).unwrap();
    }

    #[test]
    fn compressed_pages_fall_back_to_decoding() {
        use ::tiff::encoder::compression::Deflate;
        let dir = temp_dir("compressed");
        let path = dir.join("deflate.tif");
        let mut out = Cursor::new(Vec::new());
        {
            let mut enc = TiffEncoder::new(&mut out).unwrap();
            enc.write_image_with_compression::<colortype::Gray16, _>(
                2,
                2,
                Deflate::default(),
                &[9u16, 8, 7, 6],
            )
            .unwrap();
        }
        fs::File::create(&path).unwrap().write_all(&out.into_inner()).unwrap();

        let mapped = MappedTiff::from_backing(fs::read(&path).unwrap()).unwrap();
        assert!(mapped.pages()[0].byte_range.is_none());
        assert!(mapped.plane_bytes(0).is_err());
        assert_eq!(mapped.plane_u16(0).unwrap().data(), &[9, 8, 7, 6]);
        drop(mapped);
        fs::remove_dir_all(&dir).unwrap();
    }

    #[test]
    fn contiguity_check_rejects_gaps() {
        assert_eq!(contiguous_range(&[8, 12], &[4, 4], 8, 100), Some((8, 16)));
        assert_eq!(contiguous_range(&[8, 20], &[4, 4], 8, 100), None);
        assert_eq!(contiguous_range(&[8], &[4], 8, 100), None);
        assert_eq!(contiguous_range(&[90], &[10], 10, 95), None);
    }
}
