//! Zarr v2 and v3 arrays over a key-value store.
//!
//! Reading supports:
//! - Zarr v2: `.zarray` metadata, `C` and `F` chunk order, `fill_value`, dot or
//!   slash chunk separators, and the `null` and `gzip` compressors.
//! - Zarr v3: `zarr.json` metadata with a `regular` chunk grid, the default and
//!   `v2` chunk key encodings, and the `bytes` codec with optional `gzip`.
//!
//! Any other codec or compressor (for example `blosc`, `zstd`, `lz4`, `transpose`,
//! `sharding_indexed`) is reported as [`FormatError::Unsupported`] when the array
//! is opened. Zarr v3 writing is not implemented.
//!
//! Writing produces Zarr v2 arrays with `C` order and optional gzip compression.
//!
//! The store is a trait so the same code reads a directory tree, an in-memory map
//! or a remote object store. [`DirectoryStore`] is the local-directory store.

use std::fs;
use std::io::{Read, Write};
use std::path::{Path, PathBuf};

use flate2::Compression as GzLevel;
use flate2::read::GzDecoder;
use flate2::write::GzEncoder;
use serde::Deserialize;
use serde_json::{Value, json};

use super::{FormatError, FormatResult};
use crate::image::Image;

/// A key-value store holding Zarr metadata and chunk objects.
///
/// Keys use `/` as separator and never start with `/`.
pub trait ZarrStore {
    /// Reads the object at `key`, or `None` if it does not exist.
    ///
    /// # Errors
    /// Returns an I/O error if the store cannot be read.
    fn get(&self, key: &str) -> FormatResult<Option<Vec<u8>>>;

    /// Writes the object at `key`, replacing any existing object.
    ///
    /// # Errors
    /// Returns an I/O error if the store cannot be written.
    fn set(&self, key: &str, value: &[u8]) -> FormatResult<()>;
}

impl<T: ZarrStore + ?Sized> ZarrStore for &T {
    fn get(&self, key: &str) -> FormatResult<Option<Vec<u8>>> {
        (**self).get(key)
    }

    fn set(&self, key: &str, value: &[u8]) -> FormatResult<()> {
        (**self).set(key, value)
    }
}

/// A Zarr store backed by a directory: each key is a file below the root.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct DirectoryStore {
    root: PathBuf,
}

impl DirectoryStore {
    /// Creates a store rooted at `root`. The directory is created on first write.
    #[must_use]
    pub fn new(root: impl Into<PathBuf>) -> Self {
        Self { root: root.into() }
    }

    /// The root directory.
    #[must_use]
    pub fn root(&self) -> &Path {
        &self.root
    }

    fn path_of(&self, key: &str) -> FormatResult<PathBuf> {
        check_key(key)?;
        Ok(key.split('/').fold(self.root.clone(), |p, part| p.join(part)))
    }
}

fn check_key(key: &str) -> FormatResult<()> {
    if key.is_empty() || key.contains('\\') {
        return Err(FormatError::Invalid(format!("bad store key {key:?}")));
    }
    for part in key.split('/') {
        if part.is_empty() || part == "." || part == ".." {
            return Err(FormatError::Invalid(format!("bad store key {key:?}")));
        }
    }
    Ok(())
}

impl ZarrStore for DirectoryStore {
    fn get(&self, key: &str) -> FormatResult<Option<Vec<u8>>> {
        match fs::read(self.path_of(key)?) {
            Ok(bytes) => Ok(Some(bytes)),
            Err(e) if e.kind() == std::io::ErrorKind::NotFound => Ok(None),
            Err(e) => Err(e.into()),
        }
    }

    fn set(&self, key: &str, value: &[u8]) -> FormatResult<()> {
        let path = self.path_of(key)?;
        if let Some(parent) = path.parent() {
            fs::create_dir_all(parent)?;
        }
        fs::write(path, value)?;
        Ok(())
    }
}

/// Kind of a stored element type.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum DTypeKind {
    /// Unsigned integer.
    Unsigned,
    /// Signed integer.
    Signed,
    /// IEEE float.
    Float,
}

/// An array element type as stored: kind, size in bytes and byte order.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct DType {
    /// Integer or float.
    pub kind: DTypeKind,
    /// Element size in bytes.
    pub size: usize,
    /// Whether elements are stored big-endian. Irrelevant for one-byte types.
    pub big_endian: bool,
}

/// A Rust element type that a Zarr array can hold.
///
/// Implemented for the signed and unsigned integers of 8 to 64 bits, and `f32`
/// and `f64`.
pub trait ZarrElement: Copy + Default + Send + Sync + 'static {
    /// Stored kind of this type.
    const KIND: DTypeKind;
    /// Size in bytes.
    const SIZE: usize;

    /// Converts a fill value. Integer types saturate and map NaN to zero.
    fn from_f64(v: f64) -> Self;

    /// Decodes one element from at least `SIZE` bytes.
    fn decode(bytes: &[u8], big_endian: bool) -> Self;

    /// Appends the encoded element to `out`.
    fn encode(self, big_endian: bool, out: &mut Vec<u8>);
}

macro_rules! impl_zarr_element {
    ($t:ty, $kind:ident, $size:expr) => {
        impl ZarrElement for $t {
            const KIND: DTypeKind = DTypeKind::$kind;
            const SIZE: usize = $size;

            fn from_f64(v: f64) -> Self {
                v as $t
            }

            fn decode(bytes: &[u8], big_endian: bool) -> Self {
                let mut buf = [0u8; $size];
                buf.copy_from_slice(&bytes[..$size]);
                if big_endian { <$t>::from_be_bytes(buf) } else { <$t>::from_le_bytes(buf) }
            }

            fn encode(self, big_endian: bool, out: &mut Vec<u8>) {
                let bytes = if big_endian { self.to_be_bytes() } else { self.to_le_bytes() };
                out.extend_from_slice(&bytes);
            }
        }
    };
}

impl_zarr_element!(u8, Unsigned, 1);
impl_zarr_element!(i8, Signed, 1);
impl_zarr_element!(u16, Unsigned, 2);
impl_zarr_element!(i16, Signed, 2);
impl_zarr_element!(u32, Unsigned, 4);
impl_zarr_element!(i32, Signed, 4);
impl_zarr_element!(u64, Unsigned, 8);
impl_zarr_element!(i64, Signed, 8);
impl_zarr_element!(f32, Float, 4);
impl_zarr_element!(f64, Float, 8);

/// Order in which elements of a chunk (or array) are laid out.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Order {
    /// Row-major: the last dimension varies fastest.
    C,
    /// Column-major: the first dimension varies fastest.
    F,
}

/// Chunk compression.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Compressor {
    /// Chunks are stored uncompressed.
    None,
    /// Chunks are gzip-compressed (RFC 1952), with the given level in `0..=9`.
    Gzip {
        /// Compression level.
        level: u32,
    },
}

/// How chunk indices become object keys.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct ChunkKeys {
    /// Whether keys start with `c` (Zarr v3 default encoding).
    pub prefixed: bool,
    /// Separator between indices.
    pub separator: char,
}

/// Parsed array metadata, the same for Zarr v2 and v3.
#[derive(Debug, Clone, PartialEq)]
pub struct ArrayMeta {
    /// Zarr format version, 2 or 3.
    pub zarr_format: u8,
    /// Array shape.
    pub shape: Vec<usize>,
    /// Chunk shape. Edge chunks are stored at full size.
    pub chunks: Vec<usize>,
    /// Element type.
    pub dtype: DType,
    /// Chunk element order.
    pub order: Order,
    /// Chunk compression.
    pub compressor: Compressor,
    /// Value of missing chunks, as a float. A `null` fill value is read as zero.
    pub fill_value: f64,
    /// Chunk key encoding.
    pub keys: ChunkKeys,
}

#[derive(Deserialize)]
struct V2Compressor {
    id: String,
    level: Option<u32>,
}

#[derive(Deserialize)]
struct V2Array {
    zarr_format: u8,
    shape: Vec<usize>,
    chunks: Vec<usize>,
    dtype: String,
    compressor: Option<V2Compressor>,
    fill_value: Option<Value>,
    #[serde(default)]
    filters: Option<Value>,
    #[serde(default = "default_order")]
    order: String,
    dimension_separator: Option<String>,
}

fn default_order() -> String {
    "C".into()
}

#[derive(Deserialize)]
struct NamedConfig {
    name: String,
    configuration: Option<Value>,
}

#[derive(Deserialize)]
struct V3Node {
    zarr_format: u8,
    node_type: String,
    shape: Vec<usize>,
    data_type: String,
    chunk_grid: NamedConfig,
    chunk_key_encoding: Option<NamedConfig>,
    fill_value: Option<Value>,
    codecs: Vec<NamedConfig>,
}

fn parse_v2_dtype(s: &str) -> FormatResult<DType> {
    let mut chars = s.chars();
    let endian = chars.next().unwrap_or('?');
    let kind = chars.next().unwrap_or('?');
    let size: usize = chars.as_str().parse().unwrap_or(0);
    let big_endian = match endian {
        '<' | '|' => false,
        '>' => true,
        _ => return Err(FormatError::Unsupported(format!("Zarr dtype {s:?}"))),
    };
    let kind = match kind {
        'u' => DTypeKind::Unsigned,
        'i' => DTypeKind::Signed,
        'f' => DTypeKind::Float,
        _ => return Err(FormatError::Unsupported(format!("Zarr dtype {s:?}"))),
    };
    if size == 0 || !matches!(size, 1 | 2 | 4 | 8) {
        return Err(FormatError::Unsupported(format!("Zarr dtype {s:?}")));
    }
    Ok(DType { kind, size, big_endian })
}

fn parse_v3_dtype(s: &str, big_endian: bool) -> FormatResult<DType> {
    let (kind, size) = match s {
        "uint8" => (DTypeKind::Unsigned, 1),
        "uint16" => (DTypeKind::Unsigned, 2),
        "uint32" => (DTypeKind::Unsigned, 4),
        "uint64" => (DTypeKind::Unsigned, 8),
        "int8" => (DTypeKind::Signed, 1),
        "int16" => (DTypeKind::Signed, 2),
        "int32" => (DTypeKind::Signed, 4),
        "int64" => (DTypeKind::Signed, 8),
        "float32" => (DTypeKind::Float, 4),
        "float64" => (DTypeKind::Float, 8),
        other => return Err(FormatError::Unsupported(format!("Zarr v3 data type {other:?}"))),
    };
    Ok(DType { kind, size, big_endian })
}

/// Reads a fill value. `null` means no fill value, which is read as zero.
fn parse_fill(v: Option<&Value>) -> FormatResult<f64> {
    match v {
        None | Some(Value::Null) => Ok(0.0),
        Some(Value::Bool(b)) => Ok(f64::from(u8::from(*b))),
        Some(Value::Number(n)) => {
            n.as_f64().ok_or_else(|| FormatError::Invalid(format!("fill value {n} is not a float")))
        }
        Some(Value::String(s)) => match s.as_str() {
            "NaN" => Ok(f64::NAN),
            "Infinity" => Ok(f64::INFINITY),
            "-Infinity" => Ok(f64::NEG_INFINITY),
            other => Err(FormatError::Invalid(format!("fill value {other:?}"))),
        },
        Some(other) => Err(FormatError::Invalid(format!("fill value {other}"))),
    }
}

fn check_grid(shape: &[usize], chunks: &[usize]) -> FormatResult<()> {
    if shape.len() != chunks.len() {
        return Err(FormatError::Invalid("shape and chunks differ in length".into()));
    }
    if chunks.contains(&0) {
        return Err(FormatError::Invalid("chunk size of zero".into()));
    }
    Ok(())
}

fn parse_v2_meta(bytes: &[u8]) -> FormatResult<ArrayMeta> {
    let raw: V2Array = serde_json::from_slice(bytes)
        .map_err(|e| FormatError::Decode(format!("Zarr v2 .zarray: {e}")))?;
    if raw.zarr_format != 2 {
        return Err(FormatError::Invalid(format!(
            "expected zarr_format 2, got {}",
            raw.zarr_format
        )));
    }
    if raw.filters.as_ref().is_some_and(|f| !f.is_null()) {
        return Err(FormatError::Unsupported("Zarr v2 filters".into()));
    }
    check_grid(&raw.shape, &raw.chunks)?;
    let compressor = match raw.compressor {
        None => Compressor::None,
        Some(c) if c.id == "gzip" => Compressor::Gzip { level: c.level.unwrap_or(1).min(9) },
        Some(c) => return Err(FormatError::Unsupported(format!("Zarr compressor {:?}", c.id))),
    };
    let order = match raw.order.as_str() {
        "C" => Order::C,
        "F" => Order::F,
        other => return Err(FormatError::Invalid(format!("Zarr order {other:?}"))),
    };
    let separator = match raw.dimension_separator.as_deref() {
        None | Some(".") => '.',
        Some("/") => '/',
        Some(other) => return Err(FormatError::Invalid(format!("dimension_separator {other:?}"))),
    };
    Ok(ArrayMeta {
        zarr_format: 2,
        shape: raw.shape,
        chunks: raw.chunks,
        dtype: parse_v2_dtype(&raw.dtype)?,
        order,
        compressor,
        fill_value: parse_fill(raw.fill_value.as_ref())?,
        keys: ChunkKeys { prefixed: false, separator },
    })
}

fn parse_v3_meta(bytes: &[u8]) -> FormatResult<ArrayMeta> {
    let raw: V3Node = serde_json::from_slice(bytes)
        .map_err(|e| FormatError::Decode(format!("Zarr v3 zarr.json: {e}")))?;
    if raw.zarr_format != 3 {
        return Err(FormatError::Invalid(format!(
            "expected zarr_format 3, got {}",
            raw.zarr_format
        )));
    }
    if raw.node_type != "array" {
        return Err(FormatError::Invalid(format!("node_type is {:?}, not array", raw.node_type)));
    }
    if raw.chunk_grid.name != "regular" {
        return Err(FormatError::Unsupported(format!("chunk grid {:?}", raw.chunk_grid.name)));
    }
    let chunks: Vec<usize> = raw
        .chunk_grid
        .configuration
        .as_ref()
        .and_then(|c| c.get("chunk_shape"))
        .map(|v| serde_json::from_value(v.clone()))
        .transpose()
        .map_err(|e| FormatError::Invalid(format!("chunk_shape: {e}")))?
        .ok_or_else(|| FormatError::Invalid("regular chunk grid has no chunk_shape".into()))?;
    check_grid(&raw.shape, &chunks)?;

    // Codecs: an array-to-bytes codec (`bytes`) and optionally one bytes-to-bytes
    // codec (`gzip`) after it.
    let mut big_endian = false;
    let mut have_bytes = false;
    let mut compressor = Compressor::None;
    for codec in &raw.codecs {
        match codec.name.as_str() {
            "bytes" => {
                have_bytes = true;
                big_endian = codec
                    .configuration
                    .as_ref()
                    .and_then(|c| c.get("endian"))
                    .and_then(Value::as_str)
                    == Some("big");
            }
            "gzip" if have_bytes => {
                let level = codec
                    .configuration
                    .as_ref()
                    .and_then(|c| c.get("level"))
                    .and_then(Value::as_u64)
                    .unwrap_or(1);
                compressor = Compressor::Gzip { level: u32::try_from(level).unwrap_or(9).min(9) };
            }
            other => return Err(FormatError::Unsupported(format!("Zarr v3 codec {other:?}"))),
        }
    }
    if !have_bytes {
        return Err(FormatError::Unsupported("Zarr v3 array without a bytes codec".into()));
    }

    let (separator, prefixed) = match &raw.chunk_key_encoding {
        None => ('/', true),
        Some(NamedConfig { name, configuration }) => {
            let sep =
                configuration.as_ref().and_then(|c| c.get("separator")).and_then(Value::as_str);
            match name.as_str() {
                "default" => (sep.unwrap_or("/").chars().next().unwrap_or('/'), true),
                "v2" => (sep.unwrap_or(".").chars().next().unwrap_or('.'), false),
                other => {
                    return Err(FormatError::Unsupported(format!("chunk key encoding {other:?}")));
                }
            }
        }
    };

    Ok(ArrayMeta {
        zarr_format: 3,
        shape: raw.shape,
        chunks,
        dtype: parse_v3_dtype(&raw.data_type, big_endian)?,
        order: Order::C,
        compressor,
        fill_value: parse_fill(raw.fill_value.as_ref())?,
        keys: ChunkKeys { prefixed, separator },
    })
}

fn join_key(prefix: &str, name: &str) -> String {
    if prefix.is_empty() { name.to_owned() } else { format!("{prefix}/{name}") }
}

fn decompress(compressor: Compressor, raw: Vec<u8>) -> FormatResult<Vec<u8>> {
    match compressor {
        Compressor::None => Ok(raw),
        Compressor::Gzip { .. } => {
            let mut out = Vec::new();
            GzDecoder::new(raw.as_slice())
                .read_to_end(&mut out)
                .map_err(|e| FormatError::Decode(format!("gzip chunk: {e}")))?;
            Ok(out)
        }
    }
}

fn compress(compressor: Compressor, raw: &[u8]) -> FormatResult<Vec<u8>> {
    match compressor {
        Compressor::None => Ok(raw.to_vec()),
        Compressor::Gzip { level } => {
            let mut enc = GzEncoder::new(Vec::new(), GzLevel::new(level.min(9)));
            enc.write_all(raw)?;
            Ok(enc.finish()?)
        }
    }
}

/// Calls `f` for every multi-index of `dims`, in the given order. The fastest-varying
/// dimension is the last one for `C` and the first one for `F`.
fn for_each_index(
    dims: &[usize],
    order: Order,
    mut f: impl FnMut(&[usize]) -> FormatResult<()>,
) -> FormatResult<()> {
    if dims.contains(&0) {
        return Ok(());
    }
    let n = dims.len();
    let mut idx = vec![0usize; n];
    loop {
        f(&idx)?;
        let mut carried = true;
        for step in 0..n {
            let d = match order {
                Order::C => n - 1 - step,
                Order::F => step,
            };
            idx[d] += 1;
            if idx[d] < dims[d] {
                carried = false;
                break;
            }
            idx[d] = 0;
        }
        if carried {
            return Ok(());
        }
    }
}

/// C-order strides for a shape.
fn c_strides(shape: &[usize]) -> Vec<usize> {
    let mut strides = vec![1usize; shape.len()];
    for d in (0..shape.len().saturating_sub(1)).rev() {
        strides[d] = strides[d + 1] * shape[d + 1];
    }
    strides
}

/// An opened Zarr array: metadata plus the store it lives in.
#[derive(Debug)]
pub struct ZarrArray<St: ZarrStore> {
    store: St,
    prefix: String,
    meta: ArrayMeta,
}

impl<St: ZarrStore> ZarrArray<St> {
    /// Opens the array at `path` (empty for the store root). Looks for a Zarr v2
    /// `.zarray` first, then a Zarr v3 `zarr.json`.
    ///
    /// # Errors
    /// Returns [`FormatError::Invalid`] if no array metadata exists at `path`,
    /// [`FormatError::Unsupported`] for codecs or compressors this reader does not
    /// implement, and a decode error for malformed metadata.
    pub fn open(store: St, path: &str) -> FormatResult<Self> {
        let prefix = path.trim_matches('/').to_owned();
        let meta = if let Some(bytes) = store.get(&join_key(&prefix, ".zarray"))? {
            parse_v2_meta(&bytes)?
        } else if let Some(bytes) = store.get(&join_key(&prefix, "zarr.json"))? {
            parse_v3_meta(&bytes)?
        } else {
            return Err(FormatError::Invalid(format!("no Zarr array at {path:?}")));
        };
        Ok(Self { store, prefix, meta })
    }

    /// The array metadata.
    #[must_use]
    pub fn meta(&self) -> &ArrayMeta {
        &self.meta
    }

    /// The array shape.
    #[must_use]
    pub fn shape(&self) -> &[usize] {
        &self.meta.shape
    }

    /// Reads the whole array in C order.
    ///
    /// # Errors
    /// Returns [`FormatError::Invalid`] if `T` does not match the stored element
    /// type, and decode or I/O errors for chunks that cannot be read.
    pub fn read<T: ZarrElement>(&self) -> FormatResult<Vec<T>> {
        let dt = self.meta.dtype;
        if dt.kind != T::KIND || dt.size != T::SIZE {
            return Err(FormatError::Invalid(format!(
                "array holds {:?} x {} bytes, requested {:?} x {}",
                dt.kind,
                dt.size,
                T::KIND,
                T::SIZE
            )));
        }
        let shape = &self.meta.shape;
        let chunks = &self.meta.chunks;
        let total: usize = shape.iter().product();
        let mut out = vec![T::from_f64(self.meta.fill_value); total];
        if total == 0 {
            return Ok(out);
        }
        let grid: Vec<usize> = shape.iter().zip(chunks).map(|(&s, &c)| s.div_ceil(c)).collect();
        let chunk_elems: usize = chunks.iter().product();
        let strides = c_strides(shape);
        let big = dt.big_endian;

        for_each_index(&grid, Order::C, |ci| {
            let key = self.chunk_key(ci);
            let Some(stored) = self.store.get(&key)? else {
                return Ok(());
            };
            let raw = decompress(self.meta.compressor, stored)?;
            if raw.len() != chunk_elems * T::SIZE {
                return Err(FormatError::Invalid(format!(
                    "chunk {key} holds {} bytes, expected {}",
                    raw.len(),
                    chunk_elems * T::SIZE
                )));
            }
            let mut k = 0usize;
            for_each_index(chunks, self.meta.order, |j| {
                let mut lin = 0usize;
                let mut inside = true;
                for d in 0..shape.len() {
                    let g = ci[d] * chunks[d] + j[d];
                    if g >= shape[d] {
                        inside = false;
                        break;
                    }
                    lin += g * strides[d];
                }
                if inside {
                    out[lin] = T::decode(&raw[k * T::SIZE..], big);
                }
                k += 1;
                Ok(())
            })
        })?;
        Ok(out)
    }

    /// Reads a 2D (`[y, x]`) or 3D (`[z, y, x]`) array as an [`Image`].
    ///
    /// # Errors
    /// Returns [`FormatError::Invalid`] for other dimensionalities, and the errors
    /// of [`ZarrArray::read`].
    pub fn read_image<T: ZarrElement>(&self) -> FormatResult<Image<T>> {
        let s = &self.meta.shape;
        let (w, h, d) = match s.len() {
            2 => (s[1], s[0], 1),
            3 => (s[2], s[1], s[0]),
            n => {
                return Err(FormatError::Invalid(format!(
                    "a {n}-D array is not an image; use read"
                )));
            }
        };
        Ok(Image::new_3d(w, h, d, self.read::<T>()?)?)
    }

    fn chunk_key(&self, index: &[usize]) -> String {
        let sep = self.meta.keys.separator.to_string();
        let mut parts: Vec<String> = Vec::with_capacity(index.len() + 1);
        if self.meta.keys.prefixed {
            parts.push("c".into());
        }
        parts.extend(index.iter().map(usize::to_string));
        join_key(&self.prefix, &parts.join(&sep))
    }
}

/// Writes a Zarr v2 group (`.zgroup`) and, if given, its attributes (`.zattrs`).
///
/// # Errors
/// Returns an I/O error if the store cannot be written.
pub fn write_group_v2<St: ZarrStore>(
    store: &St,
    path: &str,
    attrs: Option<&Value>,
) -> FormatResult<()> {
    let prefix = path.trim_matches('/');
    let group = json!({ "zarr_format": 2 });
    store.set(&join_key(prefix, ".zgroup"), group.to_string().as_bytes())?;
    if let Some(attrs) = attrs {
        store.set(&join_key(prefix, ".zattrs"), attrs.to_string().as_bytes())?;
    }
    Ok(())
}

/// Writes an in-memory array as a Zarr v2 array with C order.
///
/// `data` is in C order and must have `product(shape)` elements. Chunks at the
/// edge are stored at full chunk size, padded with the zero fill value.
///
/// # Errors
/// Returns [`FormatError::Invalid`] for inconsistent sizes, and I/O errors.
pub fn write_array_v2<St: ZarrStore, T: ZarrElement>(
    store: &St,
    path: &str,
    shape: &[usize],
    chunks: &[usize],
    data: &[T],
    compressor: Compressor,
) -> FormatResult<()> {
    write_array_v2_order(store, path, shape, chunks, data, compressor, Order::C)
}

fn write_array_v2_order<St: ZarrStore, T: ZarrElement>(
    store: &St,
    path: &str,
    shape: &[usize],
    chunks: &[usize],
    data: &[T],
    compressor: Compressor,
    order: Order,
) -> FormatResult<()> {
    check_grid(shape, chunks)?;
    let total: usize = shape.iter().product();
    if data.len() != total {
        return Err(FormatError::Invalid(format!(
            "array has {} elements, shape needs {total}",
            data.len()
        )));
    }
    let prefix = path.trim_matches('/');
    let dtype = match (T::KIND, T::SIZE) {
        (DTypeKind::Unsigned, 1) => "|u1".to_owned(),
        (DTypeKind::Signed, 1) => "|i1".to_owned(),
        (kind, size) => {
            let k = match kind {
                DTypeKind::Unsigned => 'u',
                DTypeKind::Signed => 'i',
                DTypeKind::Float => 'f',
            };
            format!("<{k}{size}")
        }
    };
    let order_str = match order {
        Order::C => "C",
        Order::F => "F",
    };
    let compressor_json = match compressor {
        Compressor::None => Value::Null,
        Compressor::Gzip { level } => json!({ "id": "gzip", "level": level }),
    };
    let meta = json!({
        "zarr_format": 2,
        "shape": shape,
        "chunks": chunks,
        "dtype": dtype,
        "compressor": compressor_json,
        "fill_value": 0,
        "order": order_str,
        "filters": null,
        "dimension_separator": ".",
    });
    store.set(&join_key(prefix, ".zarray"), meta.to_string().as_bytes())?;

    let strides = c_strides(shape);
    let grid: Vec<usize> = shape.iter().zip(chunks).map(|(&s, &c)| s.div_ceil(c)).collect();
    for_each_index(&grid, Order::C, |ci| {
        let mut raw = Vec::with_capacity(chunks.iter().product::<usize>() * T::SIZE);
        for_each_index(chunks, order, |j| {
            let mut lin = 0usize;
            let mut inside = true;
            for d in 0..shape.len() {
                let g = ci[d] * chunks[d] + j[d];
                if g >= shape[d] {
                    inside = false;
                    break;
                }
                lin += g * strides[d];
            }
            let value = if inside { data[lin] } else { T::from_f64(0.0) };
            value.encode(false, &mut raw);
            Ok(())
        })?;
        let key_parts: Vec<String> = ci.iter().map(usize::to_string).collect();
        let key = join_key(prefix, &key_parts.join("."));
        store.set(&key, &compress(compressor, &raw)?)
    })
}

/// Writes a 2D or 3D image as a Zarr v2 array shaped `[y, x]` or `[z, y, x]`.
///
/// `chunks` gives the chunk shape in the same dimension order and must have
/// length 2 for a 2D image or 3 for a volume.
///
/// # Errors
/// Same as [`write_array_v2`].
pub fn write_image_v2<St: ZarrStore, T: ZarrElement>(
    store: &St,
    path: &str,
    img: &Image<T>,
    chunks: &[usize],
    compressor: Compressor,
) -> FormatResult<()> {
    let (w, h, d) = img.shape();
    let shape: Vec<usize> = if d == 1 { vec![h, w] } else { vec![d, h, w] };
    write_array_v2(store, path, &shape, chunks, img.data(), compressor)
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::sync::atomic::{AtomicUsize, Ordering};

    /// A unique temporary directory, removed on drop.
    struct TempDir(PathBuf);

    impl TempDir {
        fn new(tag: &str) -> Self {
            static COUNTER: AtomicUsize = AtomicUsize::new(0);
            let n = COUNTER.fetch_add(1, Ordering::SeqCst);
            let dir =
                std::env::temp_dir().join(format!("tpt-bio-zarr-{tag}-{}-{n}", std::process::id()));
            let _ = fs::remove_dir_all(&dir);
            fs::create_dir_all(&dir).unwrap();
            Self(dir)
        }

        fn path(&self) -> &Path {
            &self.0
        }
    }

    impl Drop for TempDir {
        fn drop(&mut self) {
            let _ = fs::remove_dir_all(&self.0);
        }
    }

    #[test]
    fn v2_round_trip_with_edge_chunks_and_gzip() {
        let tmp = TempDir::new("v2rt");
        let store = DirectoryStore::new(tmp.path());
        // 3 x 5 x 7 array with chunks that do not divide the shape.
        let shape = [3usize, 5, 7];
        let data: Vec<u16> = (0..105).map(|i| i as u16 * 3).collect();
        write_array_v2(&store, "a", &shape, &[2, 4, 3], &data, Compressor::Gzip { level: 5 })
            .unwrap();
        let arr = ZarrArray::open(&store, "a").unwrap();
        assert_eq!(arr.meta().zarr_format, 2);
        assert_eq!(arr.shape(), &shape);
        assert_eq!(arr.meta().compressor, Compressor::Gzip { level: 5 });
        assert_eq!(arr.read::<u16>().unwrap(), data);
        assert!(arr.read::<u8>().is_err());
    }

    #[test]
    fn v2_uncompressed_image_round_trip() {
        let tmp = TempDir::new("v2img");
        let store = DirectoryStore::new(tmp.path());
        let img = Image::new_3d(4, 3, 2, (0..24).map(|i| i as f32 * 0.5).collect()).unwrap();
        write_image_v2(&store, "group/img", &img, &[1, 2, 3], Compressor::None).unwrap();
        let arr = ZarrArray::open(&store, "group/img").unwrap();
        assert_eq!(arr.shape(), &[2, 3, 4]);
        assert_eq!(arr.read_image::<f32>().unwrap(), img);
    }

    #[test]
    fn v2_missing_chunk_reads_as_fill_value() {
        let tmp = TempDir::new("v2fill");
        let store = DirectoryStore::new(tmp.path());
        write_array_v2(&store, "a", &[2, 2], &[1, 2], &[1u8, 2, 3, 4], Compressor::None).unwrap();
        // Remove the chunk for row 1; the reader must return the fill value (zero).
        fs::remove_file(tmp.path().join("a").join("1.0")).unwrap();
        let arr = ZarrArray::open(&store, "a").unwrap();
        assert_eq!(arr.read::<u8>().unwrap(), [1, 2, 0, 0]);
    }

    #[test]
    fn v2_fortran_order_chunks_are_reordered() {
        let tmp = TempDir::new("v2f");
        let store = DirectoryStore::new(tmp.path());
        // Logical 2 x 3 array [[0,1,2],[3,4,5]], one chunk stored in F order.
        let logical: Vec<i32> = (0..6).collect();
        let meta = json!({
            "zarr_format": 2, "shape": [2, 3], "chunks": [2, 3], "dtype": "<i4",
            "compressor": null, "fill_value": null, "order": "F", "filters": null
        });
        store.set(".zarray", meta.to_string().as_bytes()).unwrap();
        let mut raw = Vec::new();
        for j in 0..3 {
            for i in 0..2 {
                logical[i * 3 + j].encode(false, &mut raw);
            }
        }
        store.set("0.0", &raw).unwrap();
        let arr = ZarrArray::open(&store, "").unwrap();
        assert_eq!(arr.meta().order, Order::F);
        assert_eq!(arr.read::<i32>().unwrap(), logical);
    }

    #[test]
    fn v2_rejects_unknown_compressor() {
        let tmp = TempDir::new("v2blosc");
        let store = DirectoryStore::new(tmp.path());
        let meta = json!({
            "zarr_format": 2, "shape": [4], "chunks": [4], "dtype": "<u1",
            "compressor": {"id": "blosc", "cname": "lz4", "clevel": 5, "shuffle": 1, "blocksize": 0},
            "fill_value": 0, "order": "C", "filters": null
        });
        store.set("a/.zarray", meta.to_string().as_bytes()).unwrap();
        assert!(matches!(ZarrArray::open(&store, "a"), Err(FormatError::Unsupported(_))));
    }

    #[test]
    fn v3_bytes_and_gzip_codecs() {
        let tmp = TempDir::new("v3");
        let store = DirectoryStore::new(tmp.path());
        let meta = json!({
            "zarr_format": 3,
            "node_type": "array",
            "shape": [2, 3],
            "data_type": "uint16",
            "chunk_grid": {"name": "regular", "configuration": {"chunk_shape": [2, 2]}},
            "chunk_key_encoding": {"name": "default", "configuration": {"separator": "/"}},
            "fill_value": 7,
            "codecs": [
                {"name": "bytes", "configuration": {"endian": "big"}},
                {"name": "gzip", "configuration": {"level": 3}}
            ]
        });
        store.set("arr/zarr.json", meta.to_string().as_bytes()).unwrap();
        // Chunk (0, 0) holds [[1, 2], [4, 5]] big-endian, gzip-compressed.
        let mut raw = Vec::new();
        for v in [1u16, 2, 4, 5] {
            v.encode(true, &mut raw);
        }
        store.set("arr/c/0/0", &compress(Compressor::Gzip { level: 3 }, &raw).unwrap()).unwrap();
        // Chunk (0, 1) is missing and reads as the fill value 7. The array is 2 x 3
        // with 2 x 2 chunks, so the second chunk column is only partly inside.
        let arr = ZarrArray::open(&store, "arr").unwrap();
        assert_eq!(arr.meta().zarr_format, 3);
        assert!(arr.meta().dtype.big_endian);
        assert_eq!(arr.read::<u16>().unwrap(), [1, 2, 7, 4, 5, 7]);
    }

    #[test]
    fn v3_unsupported_codecs_are_reported() {
        let tmp = TempDir::new("v3zstd");
        let store = DirectoryStore::new(tmp.path());
        let meta = json!({
            "zarr_format": 3, "node_type": "array", "shape": [4], "data_type": "uint8",
            "chunk_grid": {"name": "regular", "configuration": {"chunk_shape": [4]}},
            "chunk_key_encoding": {"name": "default"}, "fill_value": 0,
            "codecs": [{"name": "bytes"}, {"name": "zstd", "configuration": {"level": 1}}]
        });
        store.set("zarr.json", meta.to_string().as_bytes()).unwrap();
        match ZarrArray::open(&store, "") {
            Err(FormatError::Unsupported(msg)) => assert!(msg.contains("zstd")),
            other => panic!("expected Unsupported, got {other:?}"),
        }
    }

    #[test]
    fn group_writes_zgroup_and_attrs() {
        let tmp = TempDir::new("group");
        let store = DirectoryStore::new(tmp.path());
        write_group_v2(&store, "g", Some(&json!({"name": "plate"}))).unwrap();
        let zgroup = store.get("g/.zgroup").unwrap().unwrap();
        assert_eq!(serde_json::from_slice::<Value>(&zgroup).unwrap()["zarr_format"], 2);
        let zattrs = store.get("g/.zattrs").unwrap().unwrap();
        assert_eq!(serde_json::from_slice::<Value>(&zattrs).unwrap()["name"], "plate");
    }

    #[test]
    fn directory_store_rejects_escaping_keys() {
        let tmp = TempDir::new("keys");
        let store = DirectoryStore::new(tmp.path());
        assert!(store.get("../outside").is_err());
        assert!(store.set("a/../b", b"x").is_err());
        assert!(store.get("/abs").is_err());
    }

    #[test]
    fn for_each_index_visits_in_requested_order() {
        let mut c = Vec::new();
        for_each_index(&[2, 2], Order::C, |i| {
            c.push((i[0], i[1]));
            Ok(())
        })
        .unwrap();
        assert_eq!(c, [(0, 0), (0, 1), (1, 0), (1, 1)]);
        let mut f = Vec::new();
        for_each_index(&[2, 2], Order::F, |i| {
            f.push((i[0], i[1]));
            Ok(())
        })
        .unwrap();
        assert_eq!(f, [(0, 0), (1, 0), (0, 1), (1, 1)]);
    }
}
