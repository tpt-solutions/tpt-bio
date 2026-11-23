//! OME-TIFF reading.
//!
//! An OME-TIFF file stores OME-XML in the `ImageDescription` tag of its first IFD.
//! The XML describes one or more `Image` elements (series). Each `Pixels` element
//! gives the plane counts and the `DimensionOrder`, which says how `(z, c, t)`
//! maps to a linear plane index. `TiffData` elements map those plane indices to
//! IFD indices in the file. Without `TiffData`, the planes of series `i` follow
//! those of series `i - 1` in IFD order, which is how Bio-Formats lays them out.
//!
//! Limitations: only single-file OME-TIFF is handled. A `UUID` element is
//! ignored, so a multi-file set is read as if all its planes were in this file.
//! Pixel types without a matching [`Sample`] type (for example `int16`) can be
//! parsed but not read as pixels.

use ::tiff::tags::Tag;
use quick_xml::Reader;
use quick_xml::events::{BytesStart, Event};

use super::tiff::read_page;
use super::{FormatError, FormatResult, Sample};
use crate::image::Image;
use crate::source::ByteSource;

/// OME pixel types.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PixelType {
    /// Signed 8-bit integer.
    Int8,
    /// Unsigned 8-bit integer.
    Uint8,
    /// Signed 16-bit integer.
    Int16,
    /// Unsigned 16-bit integer.
    Uint16,
    /// Signed 32-bit integer.
    Int32,
    /// Unsigned 32-bit integer.
    Uint32,
    /// 32-bit float.
    Float,
    /// 64-bit float.
    Double,
    /// One bit per pixel.
    Bit,
    /// Complex pair of 32-bit floats.
    Complex,
    /// Complex pair of 64-bit floats.
    DoubleComplex,
}

impl PixelType {
    /// Parses the OME `Type` attribute.
    ///
    /// # Errors
    /// Returns [`FormatError::Invalid`] for unknown names.
    pub fn parse(s: &str) -> FormatResult<Self> {
        Ok(match s {
            "int8" => Self::Int8,
            "uint8" => Self::Uint8,
            "int16" => Self::Int16,
            "uint16" => Self::Uint16,
            "int32" => Self::Int32,
            "uint32" => Self::Uint32,
            "float" => Self::Float,
            "double" => Self::Double,
            "bit" => Self::Bit,
            "complex" => Self::Complex,
            "double-complex" => Self::DoubleComplex,
            other => return Err(FormatError::Invalid(format!("unknown OME pixel type {other:?}"))),
        })
    }
}

/// One of the three non-spatial dimensions of an OME image.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Dim {
    /// Focal plane.
    Z,
    /// Channel.
    C,
    /// Time point.
    T,
}

/// OME `DimensionOrder`: the order of Z, C and T from fastest to slowest, after X and Y.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum DimensionOrder {
    /// X, Y, Z, C, T.
    XYZCT,
    /// X, Y, Z, T, C.
    XYZTC,
    /// X, Y, C, T, Z.
    XYCTZ,
    /// X, Y, C, Z, T.
    XYCZT,
    /// X, Y, T, C, Z.
    XYTCZ,
    /// X, Y, T, Z, C.
    XYTZC,
}

impl DimensionOrder {
    /// Parses the OME `DimensionOrder` attribute.
    ///
    /// # Errors
    /// Returns [`FormatError::Invalid`] for anything other than the six orders.
    pub fn parse(s: &str) -> FormatResult<Self> {
        Ok(match s {
            "XYZCT" => Self::XYZCT,
            "XYZTC" => Self::XYZTC,
            "XYCTZ" => Self::XYCTZ,
            "XYCZT" => Self::XYCZT,
            "XYTCZ" => Self::XYTCZ,
            "XYTZC" => Self::XYTZC,
            other => {
                return Err(FormatError::Invalid(format!("unknown DimensionOrder {other:?}")));
            }
        })
    }

    /// The order of Z, C and T from fastest to slowest.
    #[must_use]
    pub fn fast_to_slow(self) -> [Dim; 3] {
        use Dim::{C, T, Z};
        match self {
            Self::XYZCT => [Z, C, T],
            Self::XYZTC => [Z, T, C],
            Self::XYCTZ => [C, T, Z],
            Self::XYCZT => [C, Z, T],
            Self::XYTCZ => [T, C, Z],
            Self::XYTZC => [T, Z, C],
        }
    }
}

/// One `TiffData` element: a run of planes stored in consecutive IFDs.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct TiffData {
    /// Index of the first IFD in the file's IFD chain.
    pub ifd: usize,
    /// Z index of the first plane.
    pub first_z: usize,
    /// C index of the first plane.
    pub first_c: usize,
    /// T index of the first plane.
    pub first_t: usize,
    /// Number of consecutive planes (default 1).
    pub plane_count: usize,
}

/// The `Pixels` element of an OME image.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct OmePixels {
    /// Element ID.
    pub id: String,
    /// Width in pixels (`SizeX`).
    pub size_x: usize,
    /// Height in pixels (`SizeY`).
    pub size_y: usize,
    /// Number of focal planes (`SizeZ`).
    pub size_z: usize,
    /// Number of channels (`SizeC`).
    pub size_c: usize,
    /// Number of time points (`SizeT`).
    pub size_t: usize,
    /// Plane order (`DimensionOrder`).
    pub dimension_order: DimensionOrder,
    /// Sample type (`Type`).
    pub pixel_type: PixelType,
    /// Whether stored samples are big-endian (`BigEndian`). Informational: the TIFF
    /// decoder reads byte order from the file itself.
    pub big_endian: bool,
    /// The `TiffData` elements, in document order.
    pub tiff_data: Vec<TiffData>,
}

impl OmePixels {
    /// Total number of planes, `SizeZ * SizeC * SizeT`.
    #[must_use]
    pub fn plane_count(&self) -> usize {
        self.size_z * self.size_c * self.size_t
    }

    /// Linear plane index of `(z, c, t)` in this image's dimension order, or
    /// `None` if a coordinate is out of range.
    #[must_use]
    pub fn plane_index(&self, z: usize, c: usize, t: usize) -> Option<usize> {
        if z >= self.size_z || c >= self.size_c || t >= self.size_t {
            return None;
        }
        let mut index = 0;
        let mut stride = 1;
        for dim in self.dimension_order.fast_to_slow() {
            let (coord, size) = match dim {
                Dim::Z => (z, self.size_z),
                Dim::C => (c, self.size_c),
                Dim::T => (t, self.size_t),
            };
            index += coord * stride;
            stride *= size;
        }
        Some(index)
    }
}

/// One `Image` element (a series) of an OME-XML document.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct OmeImage {
    /// Element ID.
    pub id: String,
    /// Optional `Name` attribute.
    pub name: Option<String>,
    /// The pixel description.
    pub pixels: OmePixels,
}

fn xml_err(e: impl std::fmt::Display) -> FormatError {
    FormatError::Decode(format!("OME-XML: {e}"))
}

type Attrs = Vec<(String, String)>;

fn read_attrs(e: &BytesStart<'_>) -> FormatResult<Attrs> {
    let mut out = Vec::new();
    for attr in e.attributes() {
        let attr = attr.map_err(xml_err)?;
        let key = String::from_utf8_lossy(attr.key.local_name().as_ref()).into_owned();
        let value = attr.unescape_value().map_err(xml_err)?.into_owned();
        out.push((key, value));
    }
    Ok(out)
}

fn attr<'a>(attrs: &'a Attrs, key: &str) -> Option<&'a str> {
    attrs.iter().find(|(k, _)| k == key).map(|(_, v)| v.as_str())
}

fn required<'a>(attrs: &'a Attrs, element: &str, key: &str) -> FormatResult<&'a str> {
    attr(attrs, key).ok_or_else(|| FormatError::Invalid(format!("{element} is missing {key}")))
}

fn usize_attr(attrs: &Attrs, element: &str, key: &str) -> FormatResult<usize> {
    required(attrs, element, key)?
        .parse()
        .map_err(|_| FormatError::Invalid(format!("{element} {key} is not a count")))
}

fn usize_attr_or(attrs: &Attrs, element: &str, key: &str, default: usize) -> FormatResult<usize> {
    match attr(attrs, key) {
        None => Ok(default),
        Some(_) => usize_attr(attrs, element, key),
    }
}

fn parse_pixels(attrs: &Attrs) -> FormatResult<OmePixels> {
    let e = "Pixels";
    let big_endian = attr(attrs, "BigEndian") == Some("true");
    Ok(OmePixels {
        id: attr(attrs, "ID").unwrap_or_default().to_owned(),
        size_x: usize_attr(attrs, e, "SizeX")?,
        size_y: usize_attr(attrs, e, "SizeY")?,
        size_z: usize_attr(attrs, e, "SizeZ")?,
        size_c: usize_attr(attrs, e, "SizeC")?,
        size_t: usize_attr(attrs, e, "SizeT")?,
        dimension_order: DimensionOrder::parse(required(attrs, e, "DimensionOrder")?)?,
        pixel_type: PixelType::parse(required(attrs, e, "Type")?)?,
        big_endian,
        tiff_data: Vec::new(),
    })
}

fn parse_tiff_data(attrs: &Attrs) -> FormatResult<TiffData> {
    let e = "TiffData";
    Ok(TiffData {
        ifd: usize_attr_or(attrs, e, "IFD", 0)?,
        first_z: usize_attr_or(attrs, e, "FirstZ", 0)?,
        first_c: usize_attr_or(attrs, e, "FirstC", 0)?,
        first_t: usize_attr_or(attrs, e, "FirstT", 0)?,
        plane_count: usize_attr_or(attrs, e, "PlaneCount", 1)?,
    })
}

/// Parses the `Image` elements of an OME-XML document.
///
/// # Errors
/// Returns a decode error for malformed XML, and [`FormatError::Invalid`] when a
/// required attribute is missing or has a bad value.
pub fn parse_ome_xml(xml: &str) -> FormatResult<Vec<OmeImage>> {
    let mut reader = Reader::from_str(xml);
    reader.config_mut().trim_text(true);

    let mut images = Vec::new();
    // The Image element being parsed, with its Pixels once seen.
    let mut image: Option<(String, Option<String>, Option<OmePixels>)> = None;

    loop {
        match reader.read_event().map_err(xml_err)? {
            Event::Eof => break,
            Event::Start(e) | Event::Empty(e) => {
                let name = e.local_name().as_ref().to_vec();
                if name == b"Image" {
                    let attrs = read_attrs(&e)?;
                    image = Some((
                        attr(&attrs, "ID").unwrap_or_default().to_owned(),
                        attr(&attrs, "Name").map(str::to_owned),
                        None,
                    ));
                } else if name == b"Pixels" {
                    if let Some((_, _, pixels)) = image.as_mut() {
                        *pixels = Some(parse_pixels(&read_attrs(&e)?)?);
                    }
                } else if name == b"TiffData"
                    && let Some(pixels) = image.as_mut().and_then(|(_, _, p)| p.as_mut())
                {
                    pixels.tiff_data.push(parse_tiff_data(&read_attrs(&e)?)?);
                }
            }
            Event::End(e) if e.local_name().as_ref() == b"Image" => {
                let (id, name, pixels) = image
                    .take()
                    .ok_or_else(|| FormatError::Invalid("unbalanced Image element".into()))?;
                let pixels = pixels
                    .ok_or_else(|| FormatError::Invalid(format!("Image {id} has no Pixels")))?;
                images.push(OmeImage { id, name, pixels });
            }
            _ => {}
        }
    }
    Ok(images)
}

/// An OME-TIFF file read through a [`ByteSource`].
///
/// The source is owned. Each plane read opens a fresh TIFF decoder and seeks to
/// the plane's IFD, so the struct holds no decoder state.
#[derive(Debug)]
pub struct OmeTiff<S: ByteSource> {
    source: S,
    images: Vec<OmeImage>,
    /// Per series: IFD index for each linear plane index.
    tables: Vec<Vec<Option<usize>>>,
}

impl<S: ByteSource> OmeTiff<S> {
    /// Opens an OME-TIFF from a byte source, parsing its OME-XML.
    ///
    /// # Errors
    /// Returns [`FormatError::Invalid`] if the first IFD has no OME-XML description
    /// or the plane mapping is inconsistent, and a decode error for malformed files.
    pub fn open(source: S) -> FormatResult<Self> {
        let description = {
            let mut decoder = ::tiff::decoder::Decoder::new(super::SourceReader::new(&source))?;
            decoder
                .get_tag_ascii_string(Tag::ImageDescription)
                .map_err(|e| FormatError::Invalid(format!("no ImageDescription: {e}")))?
        };
        if !description.trim_start().starts_with('<') {
            return Err(FormatError::Invalid("ImageDescription is not OME-XML".into()));
        }
        let images = parse_ome_xml(&description)?;
        if images.is_empty() {
            return Err(FormatError::Invalid("OME-XML has no Image elements".into()));
        }

        let mut tables = Vec::with_capacity(images.len());
        let mut base = 0usize;
        for image in &images {
            let p = &image.pixels;
            let total = p.plane_count();
            let mut table = vec![None; total];
            if p.tiff_data.is_empty() {
                for (i, slot) in table.iter_mut().enumerate() {
                    *slot = Some(base + i);
                }
            } else {
                for td in &p.tiff_data {
                    let start =
                        p.plane_index(td.first_z, td.first_c, td.first_t).ok_or_else(|| {
                            FormatError::Invalid(format!("TiffData start outside {}", p.id))
                        })?;
                    for k in 0..td.plane_count {
                        let slot = table.get_mut(start + k).ok_or_else(|| {
                            FormatError::Invalid(format!("TiffData runs past {}", p.id))
                        })?;
                        *slot = Some(td.ifd + k);
                    }
                }
            }
            base += total;
            tables.push(table);
        }
        Ok(Self { source, images, tables })
    }

    /// The series (`Image` elements) in document order.
    #[must_use]
    pub fn images(&self) -> &[OmeImage] {
        &self.images
    }

    /// Number of series.
    #[must_use]
    pub fn series_count(&self) -> usize {
        self.images.len()
    }

    /// The IFD index holding plane `(z, c, t)` of series `series`.
    ///
    /// # Errors
    /// Returns [`FormatError::Invalid`] if the series or coordinates are out of
    /// range, or if no IFD is mapped to that plane.
    pub fn plane_ifd(&self, series: usize, z: usize, c: usize, t: usize) -> FormatResult<usize> {
        let image = self
            .images
            .get(series)
            .ok_or_else(|| FormatError::Invalid(format!("series {series} out of range")))?;
        let index = image.pixels.plane_index(z, c, t).ok_or_else(|| {
            FormatError::Invalid(format!("plane (z={z}, c={c}, t={t}) out of range"))
        })?;
        self.tables[series][index]
            .ok_or_else(|| FormatError::Invalid(format!("plane (z={z}, c={c}, t={t}) has no IFD")))
    }

    /// Reads plane `(z, c, t)` of series `series` as a 2D image.
    ///
    /// # Errors
    /// Returns [`FormatError::Invalid`] if the pixel type of the series is not
    /// `T::PIXEL_TYPE`, if the coordinates are out of range, or if the plane's
    /// shape disagrees with the XML. Returns [`FormatError::Unsupported`] for
    /// planes with several samples per pixel.
    pub fn read_plane<T: Sample>(
        &self,
        series: usize,
        z: usize,
        c: usize,
        t: usize,
    ) -> FormatResult<Image<T>> {
        let image = self
            .images
            .get(series)
            .ok_or_else(|| FormatError::Invalid(format!("series {series} out of range")))?;
        let pixels = &image.pixels;
        if pixels.pixel_type != T::PIXEL_TYPE {
            return Err(FormatError::Invalid(format!(
                "series {series} is {:?}, not {:?}",
                pixels.pixel_type,
                T::PIXEL_TYPE
            )));
        }
        let ifd = self.plane_ifd(series, z, c, t)?;
        let mut channels = T::from_planes(read_page(&self.source, ifd)?)?;
        if channels.len() != 1 {
            return Err(FormatError::Unsupported(
                "OME-TIFF plane with several samples per pixel".into(),
            ));
        }
        let plane = channels.remove(0);
        if plane.shape() != (pixels.size_x, pixels.size_y, 1) {
            return Err(FormatError::Invalid(format!(
                "plane is {:?}, XML says {}x{}",
                plane.shape(),
                pixels.size_x,
                pixels.size_y
            )));
        }
        Ok(plane)
    }
}

impl PixelType {
    /// The OME `Type` attribute value for this pixel type.
    #[must_use]
    pub fn as_str(self) -> &'static str {
        match self {
            Self::Int8 => "int8",
            Self::Uint8 => "uint8",
            Self::Int16 => "int16",
            Self::Uint16 => "uint16",
            Self::Int32 => "int32",
            Self::Uint32 => "uint32",
            Self::Float => "float",
            Self::Double => "double",
            Self::Bit => "bit",
            Self::Complex => "complex",
            Self::DoubleComplex => "double-complex",
        }
    }
}

/// Builds the OME-XML for one single-series image of `size_x` by `size_y` pixels with
/// `size_z * size_c * size_t` planes, in `DimensionOrder="XYZCT"`.
///
/// Plane `i` of the file is `(z, c, t)` with `i = z + size_z * (c + size_c * t)`: z is
/// fastest, then c, then t. One `TiffData` element per plane maps plane `i` to IFD `i`,
/// so the XML matches the page order that [`encode_ome_tiff`] writes. `name`, when given,
/// becomes the `Image` element's `Name`, with XML escaping.
#[must_use]
pub fn ome_xml_xyzct(
    name: Option<&str>,
    (size_x, size_y): (usize, usize),
    (size_z, size_c, size_t): (usize, usize, usize),
    pixel_type: PixelType,
) -> String {
    use std::fmt::Write as _;
    let mut xml = String::new();
    xml.push_str("<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n");
    xml.push_str("<OME xmlns=\"http://www.openmicroscopy.org/Schemas/OME/2016-06\">\n");
    xml.push_str("  <Image ID=\"Image:0\"");
    if let Some(name) = name {
        // Writing to a String cannot fail, so the fmt result is ignored.
        let _ = write!(xml, " Name=\"{}\"", xml_escape(name));
    }
    xml.push_str(">\n");
    let _ = writeln!(
        xml,
        "    <Pixels ID=\"Pixels:0\" DimensionOrder=\"XYZCT\" Type=\"{}\" \
         SizeX=\"{size_x}\" SizeY=\"{size_y}\" SizeZ=\"{size_z}\" SizeC=\"{size_c}\" \
         SizeT=\"{size_t}\" BigEndian=\"false\">",
        pixel_type.as_str()
    );
    let planes = size_z * size_c * size_t;
    for i in 0..planes {
        let z = i % size_z;
        let c = (i / size_z) % size_c;
        let t = i / (size_z * size_c);
        let _ = writeln!(
            xml,
            "      <TiffData IFD=\"{i}\" FirstZ=\"{z}\" FirstC=\"{c}\" FirstT=\"{t}\" PlaneCount=\"1\"/>"
        );
    }
    xml.push_str("    </Pixels>\n  </Image>\n</OME>\n");
    xml
}

fn xml_escape(s: &str) -> String {
    let mut out = String::with_capacity(s.len());
    for ch in s.chars() {
        match ch {
            '&' => out.push_str("&amp;"),
            '<' => out.push_str("&lt;"),
            '>' => out.push_str("&gt;"),
            '"' => out.push_str("&quot;"),
            '\'' => out.push_str("&apos;"),
            other => out.push(other),
        }
    }
    out
}

/// Encodes a volume as an OME-TIFF: one TIFF page per plane, with OME-XML in the first
/// page's `ImageDescription`.
///
/// `volume` has depth `size_z * size_c * size_t`, and its `z` slices are the planes in
/// the order [`ome_xml_xyzct`] describes (z fastest, then c, then t). `big` selects
/// `BigTIFF`.
///
/// # Errors
/// Returns [`FormatError::Invalid`] if a size is zero or the depth does not equal the
/// plane count, and a decode error if the TIFF encoder fails.
pub fn encode_ome_tiff<T: Sample>(
    volume: &Image<T>,
    (size_z, size_c, size_t): (usize, usize, usize),
    name: Option<&str>,
    big: bool,
) -> FormatResult<Vec<u8>> {
    if size_z == 0 || size_c == 0 || size_t == 0 {
        return Err(FormatError::Invalid("OME sizes must be at least 1".into()));
    }
    let planes = size_z * size_c * size_t;
    if volume.depth() != planes {
        return Err(FormatError::Invalid(format!(
            "volume has {} planes, sizes give {planes}",
            volume.depth()
        )));
    }
    let xml = ome_xml_xyzct(
        name,
        (volume.width(), volume.height()),
        (size_z, size_c, size_t),
        T::PIXEL_TYPE,
    );
    super::write::encode_tiff(volume, Some(&xml), big)
}

#[cfg(test)]
mod tests {
    use super::*;
    use ::tiff::encoder::{TiffEncoder, colortype};
    use std::io::Cursor;

    const TWO_SERIES: &str = r#"<?xml version="1.0" encoding="UTF-8"?>
<OME xmlns="http://www.openmicroscopy.org/Schemas/OME/2016-06">
  <Image ID="Image:0" Name="stack">
    <Pixels ID="Pixels:0" DimensionOrder="XYZCT" Type="uint16" SizeX="2" SizeY="2"
            SizeZ="2" SizeC="2" SizeT="2" BigEndian="false"/>
  </Image>
  <Image ID="Image:1" Name="split">
    <Pixels ID="Pixels:1" DimensionOrder="XYCZT" Type="uint16" SizeX="2" SizeY="2"
            SizeZ="1" SizeC="2" SizeT="2">
      <TiffData IFD="8" PlaneCount="2" FirstT="0"/>
      <TiffData IFD="10" PlaneCount="2" FirstT="1"/>
    </Pixels>
  </Image>
</OME>"#;

    #[test]
    fn parses_images_and_tiff_data() {
        let images = parse_ome_xml(TWO_SERIES).unwrap();
        assert_eq!(images.len(), 2);
        assert_eq!(images[0].name.as_deref(), Some("stack"));
        let p0 = &images[0].pixels;
        assert_eq!((p0.size_x, p0.size_z, p0.size_c, p0.size_t), (2, 2, 2, 2));
        assert_eq!(p0.dimension_order, DimensionOrder::XYZCT);
        assert_eq!(p0.pixel_type, PixelType::Uint16);
        assert!(p0.tiff_data.is_empty());
        let td = images[1].pixels.tiff_data[1];
        assert_eq!((td.ifd, td.plane_count, td.first_t), (10, 2, 1));
    }

    #[test]
    fn plane_index_follows_dimension_order() {
        let images = parse_ome_xml(TWO_SERIES).unwrap();
        // XYZCT: Z fastest, then C, then T.
        assert_eq!(images[0].pixels.plane_index(1, 0, 0), Some(1));
        assert_eq!(images[0].pixels.plane_index(0, 1, 0), Some(2));
        assert_eq!(images[0].pixels.plane_index(0, 0, 1), Some(4));
        // XYCZT: C fastest; with SizeZ = 1, plane index = c + 2 * t.
        assert_eq!(images[1].pixels.plane_index(0, 1, 1), Some(3));
        assert_eq!(images[0].pixels.plane_index(2, 0, 0), None);
    }

    #[test]
    fn rejects_bad_dimension_order_and_type() {
        assert!(DimensionOrder::parse("XYZTX").is_err());
        assert!(PixelType::parse("uint12").is_err());
    }

    #[test]
    fn rejects_non_ome_description() {
        let xml = "<Image ID='x'><Pixels/></Image>";
        assert!(parse_ome_xml(xml).is_err());
    }

    /// Builds a two-series OME-TIFF: series 0 is a 2x2x2x2 stack in XYZCT order with
    /// implicit IFD mapping (IFDs 0-7); series 1 is 2x2 with explicit `TiffData`
    /// (IFDs 8-11, planes in XYCZT order). Each page holds its own IFD index as
    /// every sample, so a read identifies the plane it came from.
    fn build_ome_tiff(big: bool) -> Vec<u8> {
        let mut out = Cursor::new(Vec::new());
        let page_values: Vec<Vec<u16>> = (0..12u16).map(|i| vec![i; 4]).collect();
        if big {
            let mut enc = TiffEncoder::new_big(&mut out).unwrap();
            write_pages(&mut enc, &page_values, TWO_SERIES);
        } else {
            let mut enc = TiffEncoder::new(&mut out).unwrap();
            write_pages(&mut enc, &page_values, TWO_SERIES);
        }
        out.into_inner()
    }

    fn write_pages<W, K>(
        enc: &mut ::tiff::encoder::TiffEncoder<W, K>,
        pages: &[Vec<u16>],
        xml: &str,
    ) where
        W: std::io::Write + std::io::Seek,
        K: ::tiff::encoder::TiffKind,
    {
        for (i, values) in pages.iter().enumerate() {
            let mut image = enc.new_image::<colortype::Gray16>(2, 2).unwrap();
            if i == 0 {
                image.encoder().write_tag(Tag::ImageDescription, xml).unwrap();
            }
            image.write_data(values).unwrap();
        }
    }

    #[test]
    fn reads_planes_by_zct_across_series() {
        let bytes = build_ome_tiff(false);
        let ome = OmeTiff::open(bytes).unwrap();
        assert_eq!(ome.series_count(), 2);

        // Series 0, XYZCT: plane index = z + 2c + 4t, and IFD = index.
        for t in 0..2 {
            for c in 0..2 {
                for z in 0..2 {
                    let img: Image<u16> = ome.read_plane(0, z, c, t).unwrap();
                    let expected = (z + 2 * c + 4 * t) as u16;
                    assert_eq!(img.data(), &[expected; 4], "z={z} c={c} t={t}");
                }
            }
        }
        // Series 1: IFD 8 holds t=0 (c=0,1), IFD 10 holds t=1.
        let img: Image<u16> = ome.read_plane(1, 0, 1, 0).unwrap();
        assert_eq!(img.data(), &[9u16; 4]);
        let img: Image<u16> = ome.read_plane(1, 0, 0, 1).unwrap();
        assert_eq!(img.data(), &[10u16; 4]);
        let img: Image<u16> = ome.read_plane(1, 0, 1, 1).unwrap();
        assert_eq!(img.data(), &[11u16; 4]);
    }

    #[test]
    fn reads_bigtiff_ome() {
        let bytes = build_ome_tiff(true);
        let ome = OmeTiff::open(bytes).unwrap();
        let img: Image<u16> = ome.read_plane(0, 1, 1, 1).unwrap();
        assert_eq!(img.data(), &[7u16; 4]);
    }

    #[test]
    fn rejects_wrong_sample_type_and_range() {
        let ome = OmeTiff::open(build_ome_tiff(false)).unwrap();
        assert!(ome.read_plane::<u8>(0, 0, 0, 0).is_err());
        assert!(ome.read_plane::<u16>(0, 2, 0, 0).is_err());
        assert!(ome.read_plane::<u16>(5, 0, 0, 0).is_err());
    }

    #[test]
    fn rejects_file_without_ome_xml() {
        let mut out = Cursor::new(Vec::new());
        {
            let mut enc = TiffEncoder::new(&mut out).unwrap();
            enc.write_image::<colortype::Gray16>(1, 1, &[5u16]).unwrap();
        }
        let bytes = out.into_inner();
        assert!(matches!(OmeTiff::open(bytes), Err(FormatError::Invalid(_))));
    }

    #[test]
    fn written_ome_tiff_reads_back_every_plane() {
        // 3x2 planes, Z=2, C=3, T=2. Plane k holds k * 10 + pixel index.
        let (w, h, sz, sc, st) = (3usize, 2usize, 2usize, 3usize, 2usize);
        let planes = sz * sc * st;
        let data: Vec<u16> =
            (0..planes).flat_map(|k| (0..w * h).map(move |j| (k * 10 + j) as u16)).collect();
        let volume = Image::new_3d(w, h, planes, data.clone()).unwrap();
        for big in [false, true] {
            let bytes = encode_ome_tiff(&volume, (sz, sc, st), Some("a&b <c>"), big).unwrap();
            let ome = OmeTiff::open(bytes).unwrap();
            assert_eq!(ome.images()[0].name.as_deref(), Some("a&b <c>"));
            for t in 0..st {
                for c in 0..sc {
                    for z in 0..sz {
                        let k = z + sz * (c + sc * t);
                        let plane: Image<u16> = ome.read_plane(0, z, c, t).unwrap();
                        assert_eq!(
                            plane.data(),
                            &data[k * w * h..(k + 1) * w * h],
                            "z={z} c={c} t={t}"
                        );
                    }
                }
            }
        }
    }

    #[test]
    fn written_xml_maps_each_plane_to_its_page() {
        let xml = ome_xml_xyzct(None, (4, 4), (2, 2, 1), PixelType::Float);
        let images = parse_ome_xml(&xml).unwrap();
        let p = &images[0].pixels;
        assert_eq!(p.dimension_order, DimensionOrder::XYZCT);
        assert_eq!(p.pixel_type, PixelType::Float);
        assert_eq!((p.size_x, p.size_y, p.size_z, p.size_c, p.size_t), (4, 4, 2, 2, 1));
        assert_eq!(p.tiff_data.len(), 4);
        // Page 2 is (z=0, c=1, t=0): z is fastest, so index 2 = 0 + 2 * 1.
        assert_eq!(
            p.tiff_data[2],
            TiffData { ifd: 2, first_z: 0, first_c: 1, first_t: 0, plane_count: 1 }
        );
        assert_eq!(p.plane_index(1, 1, 0), Some(3));
    }

    #[test]
    fn encode_ome_tiff_rejects_plane_count_mismatch() {
        let volume = Image::new_3d(1, 1, 3, vec![0u8; 3]).unwrap();
        assert!(encode_ome_tiff(&volume, (2, 1, 1), None, false).is_err());
        assert!(encode_ome_tiff(&volume, (3, 0, 1), None, false).is_err());
    }
}
