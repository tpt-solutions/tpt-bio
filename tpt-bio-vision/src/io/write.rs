//! PNG and TIFF encoders for 8-bit, 16-bit and float images.
//!
//! PNG stores one 2D image. A 3D [`Image`] is written as a multi-page TIFF with one
//! page per `z` slice, which is the usual layout for microscopy stacks. PNG cannot
//! store float samples, so [`encode_png`] rejects `f32` images.

use std::fs;
use std::io::{Cursor, Seek, Write};
use std::path::Path;

use ::image::ImageFormat;
use ::tiff::encoder::{TiffEncoder, TiffKind, TiffKindBig, TiffKindStandard};

use super::{FormatError, FormatResult, Sample};
use crate::image::Image;

/// Encodes a 2D image as PNG bytes.
///
/// # Errors
/// Returns [`FormatError::Invalid`] if the image has depth greater than 1, and
/// [`FormatError::Unsupported`] for `f32` samples.
pub fn encode_png<T: Sample>(img: &Image<T>) -> FormatResult<Vec<u8>> {
    if img.depth() != 1 {
        return Err(FormatError::Invalid(format!(
            "PNG holds one 2D image, got depth {}",
            img.depth()
        )));
    }
    let dynamic = T::to_png(img)?;
    let mut out = Cursor::new(Vec::new());
    dynamic
        .write_to(&mut out, ImageFormat::Png)
        .map_err(|e| FormatError::Decode(format!("PNG encode: {e}")))?;
    Ok(out.into_inner())
}

/// Writes a 2D image to a PNG file.
///
/// # Errors
/// Same as [`encode_png`], and an I/O error if the file cannot be written.
pub fn save_png<T: Sample>(path: &Path, img: &Image<T>) -> FormatResult<()> {
    fs::write(path, encode_png(img)?)?;
    Ok(())
}

/// Encodes an image as TIFF bytes, one page per `z` slice.
///
/// `description`, when given, is stored as the `ImageDescription` tag of the first
/// page (the place OME-TIFF readers look for its XML). `big` selects `BigTIFF`.
///
/// # Errors
/// Returns a decode error if the encoder fails.
pub fn encode_tiff<T: Sample>(
    img: &Image<T>,
    description: Option<&str>,
    big: bool,
) -> FormatResult<Vec<u8>> {
    let mut out = Cursor::new(Vec::new());
    if big {
        write_volume::<T, _, TiffKindBig>(&mut out, img, description)?;
    } else {
        write_volume::<T, _, TiffKindStandard>(&mut out, img, description)?;
    }
    Ok(out.into_inner())
}

/// Writes an image to a TIFF file, one page per `z` slice.
///
/// # Errors
/// Same as [`encode_tiff`], and an I/O error if the file cannot be written.
pub fn save_tiff<T: Sample>(
    path: &Path,
    img: &Image<T>,
    description: Option<&str>,
    big: bool,
) -> FormatResult<()> {
    fs::write(path, encode_tiff(img, description, big)?)?;
    Ok(())
}

fn write_volume<T: Sample, W: Write + Seek, K: TiffKind>(
    writer: W,
    img: &Image<T>,
    description: Option<&str>,
) -> FormatResult<()> {
    let mut enc: TiffEncoder<W, K> = TiffEncoder::new_generic(writer)
        .map_err(|e| FormatError::Decode(format!("TIFF encode: {e}")))?;
    let (w, h, depth) = img.shape();
    let plane_len = w * h;
    for z in 0..depth {
        let page = Image::new_2d(w, h, img.data()[z * plane_len..(z + 1) * plane_len].to_vec())?;
        let desc = if z == 0 { description } else { None };
        T::write_tiff_page(&mut enc, &page, desc)?;
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::io::{Planes, png::decode_png, tiff::read_tiff};

    fn ramp_u8(w: usize, h: usize, d: usize) -> Image<u8> {
        let data: Vec<u8> = (0..w * h * d).map(|i| (i % 251) as u8).collect();
        Image::new_3d(w, h, d, data).unwrap()
    }

    #[test]
    fn png_round_trip_u8() {
        let img = ramp_u8(5, 4, 1);
        let bytes = encode_png(&img).unwrap();
        let Planes::U8(ch) = decode_png(&bytes).unwrap() else { panic!("expected u8") };
        assert_eq!(ch[0], img);
    }

    #[test]
    fn png_round_trip_u16() {
        let data: Vec<u16> = (0..12).map(|i| i * 5000).collect();
        let img = Image::new_2d(4, 3, data).unwrap();
        let bytes = encode_png(&img).unwrap();
        let Planes::U16(ch) = decode_png(&bytes).unwrap() else { panic!("expected u16") };
        assert_eq!(ch[0], img);
    }

    #[test]
    fn png_rejects_volume_and_float() {
        assert!(matches!(encode_png(&ramp_u8(2, 2, 2)), Err(FormatError::Invalid(_))));
        let f = Image::new_2d(1, 1, vec![0.5f32]).unwrap();
        assert!(matches!(encode_png(&f), Err(FormatError::Unsupported(_))));
    }

    #[test]
    fn tiff_round_trip_volume_u8() {
        let img = ramp_u8(3, 2, 4);
        let bytes = encode_tiff(&img, None, false).unwrap();
        let Planes::U8(ch) = read_tiff(&bytes[..]).unwrap() else { panic!("expected u8") };
        assert_eq!(ch[0], img);
    }

    #[test]
    fn tiff_round_trip_u16_bigtiff() {
        let data: Vec<u16> = (0..24).map(|i| u16::MAX - i * 7).collect();
        let img = Image::new_3d(4, 3, 2, data).unwrap();
        let bytes = encode_tiff(&img, None, true).unwrap();
        let Planes::U16(ch) = read_tiff(&bytes[..]).unwrap() else { panic!("expected u16") };
        assert_eq!(ch[0], img);
    }

    #[test]
    fn tiff_round_trip_f32() {
        let data = vec![0.25f32, -3.5, 1e6, f32::MIN_POSITIVE];
        let img = Image::new_2d(2, 2, data).unwrap();
        let bytes = encode_tiff(&img, Some("note"), false).unwrap();
        let Planes::F32(ch) = read_tiff(&bytes[..]).unwrap() else { panic!("expected f32") };
        assert_eq!(ch[0], img);
    }

    #[test]
    fn save_and_reload_tiff_file() {
        let dir = std::env::temp_dir().join(format!("tpt-bio-write-{}", std::process::id()));
        fs::create_dir_all(&dir).unwrap();
        let path = dir.join("stack.tif");
        let img = ramp_u8(2, 2, 2);
        save_tiff(&path, &img, None, false).unwrap();
        let bytes = fs::read(&path).unwrap();
        let Planes::U8(ch) = read_tiff(&bytes[..]).unwrap() else { panic!("expected u8") };
        assert_eq!(ch[0], img);
        fs::remove_dir_all(&dir).unwrap();
    }
}
