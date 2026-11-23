//! TIFF and `BigTIFF` reading, through the `tiff` crate.
//!
//! Striped and tiled layouts, and classic and `BigTIFF` containers, are all handled
//! by the decoder. Each sample of a pixel becomes one channel. Planar-configuration
//! files with more than one sample per pixel are rejected as unsupported.

use std::io::{Read, Seek};

use ::tiff::decoder::{Decoder, DecodingResult};
use ::tiff::tags::Tag;

use super::{FormatError, FormatResult, Planes, SourceReader, split_channels};
use crate::image::Image;
use crate::source::ByteSource;

/// Reads every page of a TIFF or `BigTIFF` file.
///
/// Pages are stacked into a volume: page `i` becomes slice `z = i`. Each channel
/// of the result is one [`Image`] of depth equal to the page count.
///
/// # Errors
/// Returns a decode error for malformed files, [`FormatError::Unsupported`] for
/// sample types and layouts that are not implemented, and [`FormatError::Invalid`]
/// when pages differ in shape or channel count.
pub fn read_tiff<S: ByteSource + ?Sized>(source: &S) -> FormatResult<Planes> {
    let mut decoder = Decoder::new(SourceReader::new(source))?;
    let mut pages = Vec::new();
    loop {
        pages.push(decode_page(&mut decoder)?);
        if !decoder.more_images() {
            break;
        }
        decoder.next_image()?;
    }
    stack_pages(pages)
}

/// Reads a single page (IFD) by its index in the file's IFD chain.
///
/// Each channel of the result has depth 1.
///
/// # Errors
/// Same as [`read_tiff`], and [`FormatError::Invalid`] if the index is past the
/// last IFD.
pub fn read_page<S: ByteSource + ?Sized>(source: &S, index: usize) -> FormatResult<Planes> {
    let mut decoder = Decoder::new(SourceReader::new(source))?;
    decoder.seek_to_image(index).map_err(|e| FormatError::Invalid(format!("IFD {index}: {e}")))?;
    decode_page(&mut decoder)
}

fn decode_page<R: Read + Seek>(decoder: &mut Decoder<R>) -> FormatResult<Planes> {
    let (w, h) = decoder.dimensions()?;
    let samples = decoder.find_tag_unsigned::<u32>(Tag::SamplesPerPixel)?.unwrap_or(1) as usize;
    let planar = decoder.find_tag_unsigned::<u16>(Tag::PlanarConfiguration)?.unwrap_or(1);
    if samples > 1 && planar != 1 {
        return Err(FormatError::Unsupported(
            "planar-configuration TIFF with several samples per pixel".into(),
        ));
    }
    let (w, h) = (w as usize, h as usize);
    match decoder.read_image()? {
        DecodingResult::U8(v) => Ok(Planes::U8(split_channels(&v, w, h, samples)?)),
        DecodingResult::U16(v) => Ok(Planes::U16(split_channels(&v, w, h, samples)?)),
        DecodingResult::F32(v) => Ok(Planes::F32(split_channels(&v, w, h, samples)?)),
        other => Err(FormatError::Unsupported(format!("TIFF sample type {}", sample_name(&other)))),
    }
}

fn sample_name(result: &DecodingResult) -> &'static str {
    match result {
        DecodingResult::U8(_) => "8-bit unsigned",
        DecodingResult::U16(_) => "16-bit unsigned",
        DecodingResult::U32(_) => "32-bit unsigned",
        DecodingResult::U64(_) => "64-bit unsigned",
        DecodingResult::F32(_) => "32-bit float",
        DecodingResult::F64(_) => "64-bit float",
        DecodingResult::I8(_) => "8-bit signed",
        DecodingResult::I16(_) => "16-bit signed",
        DecodingResult::I32(_) => "32-bit signed",
        DecodingResult::I64(_) => "64-bit signed",
    }
}

/// Stacks per-page planes into per-channel volumes.
pub(crate) fn stack_pages(pages: Vec<Planes>) -> FormatResult<Planes> {
    let mut iter = pages.into_iter();
    let first = iter.next().ok_or_else(|| FormatError::Invalid("TIFF has no pages".into()))?;
    let rest: Vec<Planes> = iter.collect();
    Ok(match first {
        Planes::U8(v) => Planes::U8(stack_typed(v, rest, |p| match p {
            Planes::U8(v) => Some(v),
            _ => None,
        })?),
        Planes::U16(v) => Planes::U16(stack_typed(v, rest, |p| match p {
            Planes::U16(v) => Some(v),
            _ => None,
        })?),
        Planes::F32(v) => Planes::F32(stack_typed(v, rest, |p| match p {
            Planes::F32(v) => Some(v),
            _ => None,
        })?),
    })
}

fn stack_typed<T: Copy>(
    first: Vec<Image<T>>,
    rest: Vec<Planes>,
    pick: fn(Planes) -> Option<Vec<Image<T>>>,
) -> FormatResult<Vec<Image<T>>> {
    let mut all = vec![first];
    for p in rest {
        all.push(
            pick(p).ok_or_else(|| FormatError::Invalid("pages differ in sample type".into()))?,
        );
    }
    let channels = all[0].len();
    let (w, h, _) = all[0].first().map_or((0, 0, 0), Image::shape);
    for page in &all {
        if page.len() != channels || page.iter().any(|c| c.shape() != (w, h, 1)) {
            return Err(FormatError::Invalid("pages differ in shape or channel count".into()));
        }
    }
    let depth = all.len();
    (0..channels)
        .map(|c| {
            let mut data = Vec::with_capacity(w * h * depth);
            for page in &all {
                data.extend_from_slice(page[c].data());
            }
            Ok(Image::new_3d(w, h, depth, data)?)
        })
        .collect()
}

#[cfg(test)]
pub(crate) mod tests {
    use super::*;
    use ::tiff::encoder::{TiffEncoder, colortype};
    use std::io::Cursor;

    /// Encodes gray `u16` pages into TIFF bytes, classic or `BigTIFF`.
    pub(crate) fn encode_gray16_pages(w: u32, h: u32, pages: &[Vec<u16>], big: bool) -> Vec<u8> {
        let mut out = Cursor::new(Vec::new());
        if big {
            let mut enc = TiffEncoder::new_big(&mut out).unwrap();
            for p in pages {
                enc.write_image::<colortype::Gray16>(w, h, p).unwrap();
            }
        } else {
            let mut enc = TiffEncoder::new(&mut out).unwrap();
            for p in pages {
                enc.write_image::<colortype::Gray16>(w, h, p).unwrap();
            }
        }
        out.into_inner()
    }

    #[test]
    fn reads_multi_page_u16_as_volume() {
        let bytes = encode_gray16_pages(2, 2, &[vec![1, 2, 3, 4], vec![5, 6, 7, 8]], false);
        let Planes::U16(ch) = read_tiff(&bytes[..]).unwrap() else { panic!("expected u16") };
        assert_eq!(ch.len(), 1);
        assert_eq!(ch[0].shape(), (2, 2, 2));
        assert_eq!(ch[0].data(), &[1, 2, 3, 4, 5, 6, 7, 8]);
    }

    #[test]
    fn reads_bigtiff() {
        let bytes = encode_gray16_pages(3, 1, &[vec![9, 8, 7]], true);
        assert_eq!(&bytes[0..2], b"II");
        let Planes::U16(ch) = read_tiff(&bytes[..]).unwrap() else { panic!("expected u16") };
        assert_eq!(ch[0].data(), &[9, 8, 7]);
    }

    #[test]
    fn reads_single_page_by_index() {
        let bytes = encode_gray16_pages(2, 1, &[vec![10, 11], vec![20, 21]], false);
        let Planes::U16(ch) = read_page(&bytes[..], 1).unwrap() else { panic!("expected u16") };
        assert_eq!(ch[0].shape(), (2, 1, 1));
        assert_eq!(ch[0].data(), &[20, 21]);
        assert!(read_page(&bytes[..], 5).is_err());
    }

    #[test]
    fn reads_rgb_as_three_channels() {
        let mut out = Cursor::new(Vec::new());
        {
            let mut enc = TiffEncoder::new(&mut out).unwrap();
            enc.write_image::<colortype::RGB8>(2, 1, &[1u8, 2, 3, 4, 5, 6]).unwrap();
        }
        let bytes = out.into_inner();
        let Planes::U8(ch) = read_tiff(&bytes[..]).unwrap() else { panic!("expected u8") };
        assert_eq!(ch.len(), 3);
        assert_eq!(ch[0].data(), &[1, 4]);
        assert_eq!(ch[2].data(), &[3, 6]);
    }

    #[test]
    fn reads_float_pages() {
        let mut out = Cursor::new(Vec::new());
        {
            let mut enc = TiffEncoder::new(&mut out).unwrap();
            enc.write_image::<colortype::Gray32Float>(2, 1, &[0.5f32, -1.25]).unwrap();
        }
        let bytes = out.into_inner();
        let Planes::F32(ch) = read_tiff(&bytes[..]).unwrap() else { panic!("expected f32") };
        assert_eq!(ch[0].data(), &[0.5, -1.25]);
    }

    /// Builds a little-endian, uncompressed, 8-bit TIFF with 16x16 tiles. The image
    /// is `w` x 16 and `pixels` holds it row-major. Edge tiles are zero-padded.
    pub(crate) fn tiled_u8_tiff(w: u32, pixels: &[u8]) -> Vec<u8> {
        const TILE: usize = 16;
        let tiles_across = (w as usize).div_ceil(TILE);
        let mut tiles = Vec::new();
        for tx in 0..tiles_across {
            for r in 0..TILE {
                for c in 0..TILE {
                    let x = tx * TILE + c;
                    tiles.push(if x < w as usize { pixels[r * w as usize + x] } else { 0 });
                }
            }
        }
        let tile_bytes = TILE * TILE;
        let data_start = 8usize;
        let ifd_start = data_start + tiles.len();
        let offsets_start = ifd_start + 2 + 10 * 12 + 4;
        let counts_start = offsets_start + 4 * tiles_across;

        let mut out = b"II*\0".to_vec();
        out.extend_from_slice(&(ifd_start as u32).to_le_bytes());
        out.extend_from_slice(&tiles);
        assert_eq!(out.len(), ifd_start);

        let entry = |tag: u16, ty: u16, count: u32, value: u32| {
            let mut e = Vec::new();
            e.extend_from_slice(&tag.to_le_bytes());
            e.extend_from_slice(&ty.to_le_bytes());
            e.extend_from_slice(&count.to_le_bytes());
            e.extend_from_slice(&value.to_le_bytes());
            e
        };
        // SHORT values with count 1 sit in the low two bytes of the value field.
        let short = |v: u16| u32::from(v);
        out.extend_from_slice(&10u16.to_le_bytes());
        out.extend(entry(256, 4, 1, w));
        out.extend(entry(257, 4, 1, 16));
        out.extend(entry(258, 3, 1, short(8)));
        out.extend(entry(259, 3, 1, short(1)));
        out.extend(entry(262, 3, 1, short(1)));
        out.extend(entry(277, 3, 1, short(1)));
        out.extend(entry(322, 4, 1, TILE as u32));
        out.extend(entry(323, 4, 1, TILE as u32));
        out.extend(entry(324, 4, tiles_across as u32, offsets_start as u32));
        out.extend(entry(325, 4, tiles_across as u32, counts_start as u32));
        out.extend_from_slice(&0u32.to_le_bytes());
        for tx in 0..tiles_across {
            out.extend_from_slice(&((data_start + tx * tile_bytes) as u32).to_le_bytes());
        }
        for _ in 0..tiles_across {
            out.extend_from_slice(&(tile_bytes as u32).to_le_bytes());
        }
        out
    }

    #[test]
    fn reads_tiled_tiff() {
        // 20 x 16 image: two 16 x 16 tiles across, the second one partly padding.
        let w = 20u32;
        let pixels: Vec<u8> = (0..w * 16).map(|i| (i % 200) as u8).collect();
        let bytes = tiled_u8_tiff(w, &pixels);
        let Planes::U8(ch) = read_tiff(&bytes[..]).unwrap() else { panic!("expected u8") };
        assert_eq!(ch[0].shape(), (20, 16, 1));
        assert_eq!(ch[0].data(), &pixels[..]);
    }

    #[test]
    fn rejects_truncated_file() {
        let bytes = encode_gray16_pages(2, 2, &[vec![1, 2, 3, 4]], false);
        assert!(read_tiff(&bytes[..10]).is_err());
    }
}
