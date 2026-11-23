//! PNG reading. Each channel of the PNG becomes one plane.
//!
//! Greyscale, grey-alpha, RGB and RGBA PNGs at 8 or 16 bits per sample are supported.
//! A multi-channel PNG returns one plane per channel in file order (for RGBA:
//! red, green, blue, alpha). Palette images are expanded by the decoder to RGB or RGBA.

use std::io::Cursor;

use ::image::{DynamicImage, ImageFormat, ImageReader};

use super::{FormatError, FormatResult, Planes, decode_err, read_all, split_channels};
use crate::source::ByteSource;

/// Reads a PNG from a byte source.
///
/// # Errors
/// Returns an I/O error, a decode error, or [`FormatError::Unsupported`] for
/// sample layouts this reader does not map.
pub fn read_png<S: ByteSource + ?Sized>(source: &S) -> FormatResult<Planes> {
    decode_png(&read_all(source)?)
}

/// Decodes PNG bytes into planes.
///
/// # Errors
/// Returns a decode error if the bytes are not a valid PNG.
pub fn decode_png(bytes: &[u8]) -> FormatResult<Planes> {
    let img = ImageReader::with_format(Cursor::new(bytes), ImageFormat::Png)
        .decode()
        .map_err(decode_err)?;
    dynamic_to_planes(img)
}

fn dynamic_to_planes(img: DynamicImage) -> FormatResult<Planes> {
    let w = img.width() as usize;
    let h = img.height() as usize;
    let planes = match img {
        DynamicImage::ImageLuma8(b) => Planes::U8(split_channels(b.as_raw().as_slice(), w, h, 1)?),
        DynamicImage::ImageLumaA8(b) => Planes::U8(split_channels(b.as_raw().as_slice(), w, h, 2)?),
        DynamicImage::ImageRgb8(b) => Planes::U8(split_channels(b.as_raw().as_slice(), w, h, 3)?),
        DynamicImage::ImageRgba8(b) => Planes::U8(split_channels(b.as_raw().as_slice(), w, h, 4)?),
        DynamicImage::ImageLuma16(b) => {
            Planes::U16(split_channels(b.as_raw().as_slice(), w, h, 1)?)
        }
        DynamicImage::ImageLumaA16(b) => {
            Planes::U16(split_channels(b.as_raw().as_slice(), w, h, 2)?)
        }
        DynamicImage::ImageRgb16(b) => Planes::U16(split_channels(b.as_raw().as_slice(), w, h, 3)?),
        DynamicImage::ImageRgba16(b) => {
            Planes::U16(split_channels(b.as_raw().as_slice(), w, h, 4)?)
        }
        other => {
            return Err(FormatError::Unsupported(format!("PNG colour type {:?}", other.color())));
        }
    };
    Ok(planes)
}

#[cfg(test)]
mod tests {
    use super::*;
    use ::image::{ImageBuffer, Luma, Rgb};

    fn png_bytes(img: &DynamicImage) -> Vec<u8> {
        let mut out = Cursor::new(Vec::new());
        img.write_to(&mut out, ImageFormat::Png).unwrap();
        out.into_inner()
    }

    #[test]
    fn reads_8_bit_greyscale() {
        let buf = ImageBuffer::<Luma<u8>, _>::from_raw(3, 2, vec![1u8, 2, 3, 4, 5, 6]).unwrap();
        let planes = read_png(&png_bytes(&DynamicImage::ImageLuma8(buf))).unwrap();
        let Planes::U8(ch) = planes else { panic!("expected u8") };
        assert_eq!(ch.len(), 1);
        assert_eq!(ch[0].shape(), (3, 2, 1));
        assert_eq!(ch[0].data(), &[1, 2, 3, 4, 5, 6]);
    }

    #[test]
    fn reads_16_bit_greyscale() {
        let buf = ImageBuffer::<Luma<u16>, _>::from_raw(2, 1, vec![300u16, 65535]).unwrap();
        let Planes::U16(ch) = read_png(&png_bytes(&DynamicImage::ImageLuma16(buf))).unwrap() else {
            panic!("expected u16")
        };
        assert_eq!(ch[0].data(), &[300, 65535]);
    }

    #[test]
    fn splits_rgb_into_three_planes() {
        let buf =
            ImageBuffer::<Rgb<u8>, _>::from_raw(2, 1, vec![10u8, 20, 30, 40, 50, 60]).unwrap();
        let Planes::U8(ch) = read_png(&png_bytes(&DynamicImage::ImageRgb8(buf))).unwrap() else {
            panic!("expected u8")
        };
        assert_eq!(ch.len(), 3);
        assert_eq!(ch[0].data(), &[10, 40]);
        assert_eq!(ch[1].data(), &[20, 50]);
        assert_eq!(ch[2].data(), &[30, 60]);
    }

    #[test]
    fn rejects_non_png() {
        assert!(decode_png(b"not a png at all").is_err());
    }
}
