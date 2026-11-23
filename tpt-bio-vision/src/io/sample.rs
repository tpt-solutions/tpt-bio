//! The sample types that image formats can read and write.

use std::io::{Seek, Write};

use ::image::{DynamicImage, ImageBuffer, Luma};
use ::tiff::encoder::{TiffEncoder, TiffKind, colortype};
use ::tiff::tags::Tag;

use super::ome::PixelType;
use super::{FormatError, FormatResult, Planes};
use crate::image::Image;

mod sealed {
    /// Prevents implementation of [`super::Sample`] outside this crate.
    pub trait Sealed {}
}

/// A sample type that can be read from and written to image files.
///
/// Implemented for `u8`, `u16` and `f32`. The methods are plumbing for the format
/// modules; callers use the typed readers and writers instead.
pub trait Sample: Copy + Default + Send + Sync + 'static + sealed::Sealed {
    /// OME-XML pixel type for this sample type.
    const PIXEL_TYPE: PixelType;

    /// Extracts channels of this type, or fails if `planes` holds another type.
    ///
    /// # Errors
    /// Returns [`FormatError::Invalid`] when the sample type does not match.
    fn from_planes(planes: Planes) -> FormatResult<Vec<Image<Self>>>;

    /// Wraps channels of this type.
    fn into_planes(channels: Vec<Image<Self>>) -> Planes;

    /// Appends one 2D page to a TIFF encoder, with an optional `ImageDescription`.
    ///
    /// # Errors
    /// Returns a decode error if the encoder fails.
    fn write_tiff_page<W: Write + Seek, K: TiffKind>(
        enc: &mut TiffEncoder<W, K>,
        page: &Image<Self>,
        description: Option<&str>,
    ) -> FormatResult<()>;

    /// Converts a 2D image to a PNG-compatible dynamic image.
    ///
    /// # Errors
    /// Returns [`FormatError::Unsupported`] for types PNG cannot store.
    fn to_png(page: &Image<Self>) -> FormatResult<DynamicImage>;
}

fn planes_name(planes: &Planes) -> &'static str {
    match planes {
        Planes::U8(_) => "8-bit",
        Planes::U16(_) => "16-bit",
        Planes::F32(_) => "float",
    }
}

fn mismatch(expected: &str, planes: &Planes) -> FormatError {
    FormatError::Invalid(format!("expected {expected} samples, found {}", planes_name(planes)))
}

fn tiff_err(e: impl std::fmt::Display) -> FormatError {
    FormatError::Decode(format!("TIFF encode: {e}"))
}

impl sealed::Sealed for u8 {}
impl sealed::Sealed for u16 {}
impl sealed::Sealed for f32 {}

impl Sample for u8 {
    const PIXEL_TYPE: PixelType = PixelType::Uint8;

    fn from_planes(planes: Planes) -> FormatResult<Vec<Image<u8>>> {
        match planes {
            Planes::U8(v) => Ok(v),
            other => Err(mismatch("8-bit", &other)),
        }
    }

    fn into_planes(channels: Vec<Image<u8>>) -> Planes {
        Planes::U8(channels)
    }

    fn write_tiff_page<W: Write + Seek, K: TiffKind>(
        enc: &mut TiffEncoder<W, K>,
        page: &Image<u8>,
        description: Option<&str>,
    ) -> FormatResult<()> {
        let (w, h, _) = page.shape();
        let mut image = enc.new_image::<colortype::Gray8>(w as u32, h as u32).map_err(tiff_err)?;
        if let Some(d) = description {
            image.encoder().write_tag(Tag::ImageDescription, d).map_err(tiff_err)?;
        }
        image.write_data(page.data()).map_err(tiff_err)
    }

    fn to_png(page: &Image<u8>) -> FormatResult<DynamicImage> {
        let (w, h, _) = page.shape();
        let buf =
            ImageBuffer::<Luma<u8>, Vec<u8>>::from_raw(w as u32, h as u32, page.data().to_vec())
                .ok_or_else(|| FormatError::Invalid("PNG buffer size mismatch".into()))?;
        Ok(DynamicImage::ImageLuma8(buf))
    }
}

impl Sample for u16 {
    const PIXEL_TYPE: PixelType = PixelType::Uint16;

    fn from_planes(planes: Planes) -> FormatResult<Vec<Image<u16>>> {
        match planes {
            Planes::U16(v) => Ok(v),
            other => Err(mismatch("16-bit", &other)),
        }
    }

    fn into_planes(channels: Vec<Image<u16>>) -> Planes {
        Planes::U16(channels)
    }

    fn write_tiff_page<W: Write + Seek, K: TiffKind>(
        enc: &mut TiffEncoder<W, K>,
        page: &Image<u16>,
        description: Option<&str>,
    ) -> FormatResult<()> {
        let (w, h, _) = page.shape();
        let mut image = enc.new_image::<colortype::Gray16>(w as u32, h as u32).map_err(tiff_err)?;
        if let Some(d) = description {
            image.encoder().write_tag(Tag::ImageDescription, d).map_err(tiff_err)?;
        }
        image.write_data(page.data()).map_err(tiff_err)
    }

    fn to_png(page: &Image<u16>) -> FormatResult<DynamicImage> {
        let (w, h, _) = page.shape();
        let buf =
            ImageBuffer::<Luma<u16>, Vec<u16>>::from_raw(w as u32, h as u32, page.data().to_vec())
                .ok_or_else(|| FormatError::Invalid("PNG buffer size mismatch".into()))?;
        Ok(DynamicImage::ImageLuma16(buf))
    }
}

impl Sample for f32 {
    const PIXEL_TYPE: PixelType = PixelType::Float;

    fn from_planes(planes: Planes) -> FormatResult<Vec<Image<f32>>> {
        match planes {
            Planes::F32(v) => Ok(v),
            other => Err(mismatch("float", &other)),
        }
    }

    fn into_planes(channels: Vec<Image<f32>>) -> Planes {
        Planes::F32(channels)
    }

    fn write_tiff_page<W: Write + Seek, K: TiffKind>(
        enc: &mut TiffEncoder<W, K>,
        page: &Image<f32>,
        description: Option<&str>,
    ) -> FormatResult<()> {
        let (w, h, _) = page.shape();
        let mut image =
            enc.new_image::<colortype::Gray32Float>(w as u32, h as u32).map_err(tiff_err)?;
        if let Some(d) = description {
            image.encoder().write_tag(Tag::ImageDescription, d).map_err(tiff_err)?;
        }
        image.write_data(page.data()).map_err(tiff_err)
    }

    fn to_png(_page: &Image<f32>) -> FormatResult<DynamicImage> {
        Err(FormatError::Unsupported("PNG cannot store float samples".into()))
    }
}
