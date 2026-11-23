//! `PyO3` bindings for tpt-bio.
//!
//! The module is `tpt_bio`. Heavy work runs under `py.detach`, so the GIL is
//! released while it executes. Sequence readers (FASTQ, indexed BAM) and
//! vision routines (Otsu threshold, thresholding, 2D labelling, watershed) are
//! exposed; NumPy arrays are used for image inputs and outputs.

use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter};
use std::path::Path;
use std::sync::Mutex;

use numpy::ndarray::{Array1, Array2};
use numpy::{Element, PyArray1, PyArray2, PyReadonlyArray2};
use pyo3::exceptions::{PyIOError, PyRuntimeError, PyValueError};
use pyo3::prelude::*;
use pyo3::types::{PyBytes, PyDict, PyList};

use tpt_bio_core::coords::Region;
use tpt_bio_pipeline::{Pipeline as TptPipeline, Value as PipelineValue};
use tpt_bio_seq::bai::{self, Index};
use tpt_bio_seq::fastq::{self, FastqBuf};
use tpt_bio_seq::{BamReader, Record, SeqError};
use tpt_bio_vision::features::{intensity_features, shape_features_2d};
use tpt_bio_vision::image::{Image, ImageError};
use tpt_bio_vision::io::{Planes, png, tiff};
use tpt_bio_vision::label::{self, Connectivity2};
use tpt_bio_vision::source::FileSource;
use tpt_bio_vision::{threshold, watershed};

fn seq_err(e: SeqError) -> PyErr {
    match e {
        SeqError::Io(io) => PyIOError::new_err(io.to_string()),
        other => PyValueError::new_err(other.to_string()),
    }
}

fn image_err(e: ImageError) -> PyErr {
    PyValueError::new_err(e.to_string())
}

/// Reverse complement of an IUPAC DNA sequence. Case is preserved.
#[pyfunction]
fn reverse_complement(py: Python<'_>, seq: &[u8]) -> PyResult<Vec<u8>> {
    py.detach(|| tpt_bio_core::nucleotide::reverse_complement(seq))
        .map_err(|e| PyValueError::new_err(e.to_string()))
}

/// Fraction of unambiguous bases that are G or C, or `None` if there are none.
#[pyfunction]
fn gc_content(py: Python<'_>, seq: &[u8]) -> Option<f64> {
    py.detach(|| tpt_bio_core::stats::gc_content(seq))
}

/// Canonical k-mers of a sequence, as `(offset, packed)` pairs.
#[pyfunction]
fn canonical_kmers(py: Python<'_>, seq: &[u8], k: usize) -> PyResult<Vec<(usize, u64)>> {
    py.detach(|| {
        tpt_bio_core::kmer::Kmers::new(seq, k)
            .map(|it| it.canonical().collect())
            .map_err(|e| PyValueError::new_err(e.to_string()))
    })
}

// ---------------------------------------------------------------------------
// FASTQ
// ---------------------------------------------------------------------------

/// One FASTQ record. `name` excludes the leading `@`; `qual` is Phred+33 text.
#[pyclass(name = "FastqRecord", module = "tpt_bio", frozen)]
pub struct PyFastqRecord {
    #[pyo3(get)]
    name: Py<PyBytes>,
    #[pyo3(get)]
    seq: Py<PyBytes>,
    #[pyo3(get)]
    qual: Py<PyBytes>,
}

#[pymethods]
impl PyFastqRecord {
    fn __repr__(&self, py: Python<'_>) -> String {
        format!(
            "FastqRecord(name={}, len={})",
            String::from_utf8_lossy(self.name.bind(py).as_bytes()),
            self.seq.bind(py).len().unwrap_or(0)
        )
    }
}

fn to_py_record(py: Python<'_>, buf: &FastqBuf) -> PyFastqRecord {
    PyFastqRecord {
        name: PyBytes::new(py, &buf.name).unbind(),
        seq: PyBytes::new(py, &buf.seq).unbind(),
        qual: PyBytes::new(py, &buf.qual).unbind(),
    }
}

/// Streaming FASTQ reader. Gzip-compressed files are detected from their magic
/// bytes. Iterating yields `FastqRecord` objects.
#[pyclass(name = "FastqReader", module = "tpt_bio")]
pub struct PyFastqReader {
    inner: Mutex<fastq::FastqReader<Box<dyn BufRead + Send>>>,
    buf: FastqBuf,
}

fn lock_err<T>(_: std::sync::PoisonError<T>) -> PyErr {
    PyValueError::new_err("FastqReader lock poisoned by an earlier panic")
}

#[pymethods]
impl PyFastqReader {
    #[new]
    fn new(py: Python<'_>, path: &str) -> PyResult<Self> {
        let src = py.detach(|| fastq::open_reader(Path::new(path))).map_err(seq_err)?;
        Ok(Self { inner: Mutex::new(fastq::FastqReader::new(src)), buf: FastqBuf::default() })
    }

    fn __iter__(slf: PyRef<'_, Self>) -> PyRef<'_, Self> {
        slf
    }

    fn __next__(&mut self, py: Python<'_>) -> PyResult<Option<PyFastqRecord>> {
        let inner = &self.inner;
        let buf = &mut self.buf;
        let more = py.detach(|| -> Result<bool, PyErr> {
            let mut guard = inner.lock().map_err(lock_err)?;
            guard.read_record(buf).map_err(seq_err)
        })?;
        Ok(more.then(|| to_py_record(py, &self.buf)))
    }

    /// Reads up to `n` records. Returns fewer (or none) at end of file.
    fn read_batch(&mut self, py: Python<'_>, n: usize) -> PyResult<Vec<PyFastqRecord>> {
        let inner = &self.inner;
        let records = py.detach(|| -> Result<Vec<FastqBuf>, PyErr> {
            let mut guard = inner.lock().map_err(lock_err)?;
            let mut out = Vec::new();
            while out.len() < n {
                let mut buf = FastqBuf::default();
                if !guard.read_record(&mut buf).map_err(seq_err)? {
                    break;
                }
                out.push(buf);
            }
            Ok(out)
        })?;
        Ok(records.iter().map(|b| to_py_record(py, b)).collect())
    }
}

// ---------------------------------------------------------------------------
// Indexed BAM
// ---------------------------------------------------------------------------

/// One alignment record from a BAM file.
///
/// `pos` is 0-based (`-1` when unset). `seq` is ASCII bases. `qual` holds raw
/// Phred scores (not +33); empty means unavailable.
#[pyclass(name = "BamRecord", module = "tpt_bio", frozen)]
pub struct PyBamRecord {
    #[pyo3(get)]
    name: Py<PyBytes>,
    #[pyo3(get)]
    flags: u16,
    #[pyo3(get)]
    ref_id: i32,
    #[pyo3(get)]
    pos: i64,
    #[pyo3(get)]
    mapq: u8,
    #[pyo3(get)]
    cigar: String,
    #[pyo3(get)]
    next_ref_id: i32,
    #[pyo3(get)]
    next_pos: i64,
    #[pyo3(get)]
    tlen: i64,
    #[pyo3(get)]
    seq: Py<PyBytes>,
    #[pyo3(get)]
    qual: Py<PyBytes>,
}

fn to_py_bam_record(py: Python<'_>, r: &Record) -> PyBamRecord {
    PyBamRecord {
        name: PyBytes::new(py, &r.name).unbind(),
        flags: r.flags.0,
        ref_id: r.ref_id,
        pos: r.pos,
        mapq: r.mapq,
        cigar: tpt_bio_seq::record::format_cigar(&r.cigar),
        next_ref_id: r.next_ref_id,
        next_pos: r.next_pos,
        tlen: r.tlen,
        seq: PyBytes::new(py, &r.seq).unbind(),
        qual: PyBytes::new(py, &r.qual).unbind(),
    }
}

/// A coordinate-sorted BAM file with its BAI index, supporting region queries.
#[pyclass(name = "IndexedBam", module = "tpt_bio")]
pub struct PyIndexedBam {
    reader: BamReader<BufReader<File>>,
    index: Index,
}

#[pymethods]
impl PyIndexedBam {
    #[new]
    fn new(py: Python<'_>, bam_path: &str, bai_path: &str) -> PyResult<Self> {
        py.detach(|| -> Result<Self, SeqError> {
            let reader = BamReader::new(BufReader::new(File::open(bam_path)?))?;
            let index = Index::read_from(BufReader::new(File::open(bai_path)?))?;
            Ok(Self { reader, index })
        })
        .map_err(seq_err)
    }

    /// Reference names and lengths, in header order.
    #[getter]
    fn references(&self) -> Vec<(String, u64)> {
        self.reader.header().references().iter().map(|r| (r.name.clone(), r.len)).collect()
    }

    /// Records overlapping `region`, given as `contig`, `contig:start` or
    /// `contig:start-end` (1-based, inclusive).
    fn fetch(&mut self, py: Python<'_>, region: &str) -> PyResult<Vec<PyBamRecord>> {
        let parsed = Region::parse(region).map_err(|e| PyValueError::new_err(e.to_string()))?;
        let ref_id = self.reader.header().ref_id(parsed.contig()).ok_or_else(|| {
            PyValueError::new_err(format!("unknown reference: {}", parsed.contig()))
        })?;
        let (beg, end) = match parsed.interval() {
            Some(iv) => (iv.start(), iv.end()),
            None => {
                let len = self.reader.header().references()[ref_id as usize].len;
                (0, len)
            }
        };
        let reader = &mut self.reader;
        let index = &self.index;
        let records =
            py.detach(|| bai::fetch(reader, index, ref_id as usize, beg, end)).map_err(seq_err)?;
        Ok(records.iter().map(|r| to_py_bam_record(py, r)).collect())
    }
}

/// Builds a BAI index for a coordinate-sorted BAM file and writes it to `bai_path`.
#[pyfunction]
fn build_bai(py: Python<'_>, bam_path: &str, bai_path: &str) -> PyResult<()> {
    py.detach(|| -> Result<(), SeqError> {
        let mut reader = BamReader::new(BufReader::new(File::open(bam_path)?))?;
        let index = bai::build_index(&mut reader)?;
        index.write_to(BufWriter::new(File::create(bai_path)?))?;
        Ok(())
    })
    .map_err(seq_err)
}

// ---------------------------------------------------------------------------
// Vision
// ---------------------------------------------------------------------------

/// Copies a 2D NumPy array into a tpt-bio image. Images own their samples, so
/// this is one copy on input; NumPy arrays of any strides are accepted.
fn image_from<T: Element + Copy>(a: &PyReadonlyArray2<'_, T>) -> PyResult<Image<T>> {
    let view = a.as_array();
    let (h, w) = view.dim();
    let data: Vec<T> = view.iter().copied().collect();
    Image::new_2d(w, h, data).map_err(image_err)
}

/// Converts any supported 2D numeric array to an `f64` image.
fn image_f64(obj: &Bound<'_, PyAny>) -> PyResult<Image<f64>> {
    macro_rules! try_types {
        ($($t:ty),*) => {$(
            if let Ok(a) = obj.extract::<PyReadonlyArray2<'_, $t>>() {
                return Ok(image_from(&a)?.map(|&v| f64::from(v)));
            }
        )*};
    }
    try_types!(u8, u16, u32, i8, i16, i32, f32, f64);
    Err(PyValueError::new_err(
        "expected a 2D array of dtype uint8, uint16, uint32, int8, int16, int32, float32 or float64",
    ))
}

fn label_array<'py>(py: Python<'py>, labels: Image<u32>) -> PyResult<Bound<'py, PyArray2<u32>>> {
    let (w, h, _) = labels.shape();
    let arr = Array2::from_shape_vec((h, w), labels.into_data())
        .map_err(|e| PyValueError::new_err(e.to_string()))?;
    Ok(PyArray2::from_owned_array(py, arr))
}

/// Otsu threshold of a 2D image, or `None` if the image is empty or constant.
#[pyfunction]
#[pyo3(signature = (image, bins=256))]
fn otsu_threshold(py: Python<'_>, image: &Bound<'_, PyAny>, bins: usize) -> PyResult<Option<f64>> {
    let img = image_f64(image)?;
    py.detach(|| threshold::otsu_threshold(&img, bins)).map_err(image_err)
}

/// Foreground mask: `True` where `image > threshold_value`.
#[pyfunction]
fn apply_threshold<'py>(
    py: Python<'py>,
    image: &Bound<'py, PyAny>,
    threshold_value: f64,
) -> PyResult<Bound<'py, PyArray2<bool>>> {
    let img = image_f64(image)?;
    let mask = py.detach(|| threshold::apply_threshold(&img, threshold_value));
    let (w, h, _) = mask.shape();
    let arr = Array2::from_shape_vec((h, w), mask.into_data())
        .map_err(|e| PyValueError::new_err(e.to_string()))?;
    Ok(PyArray2::from_owned_array(py, arr))
}

/// Labels connected foreground components of a 2D boolean mask.
///
/// `connectivity` is 4 (edge-sharing) or 8 (edge- and corner-sharing).
/// Returns `(labels, count)` where labels is a `uint32` array with background 0.
#[pyfunction]
#[pyo3(signature = (mask, connectivity=8))]
fn label_components<'py>(
    py: Python<'py>,
    mask: PyReadonlyArray2<'py, bool>,
    connectivity: u8,
) -> PyResult<(Bound<'py, PyArray2<u32>>, u32)> {
    let conn = match connectivity {
        4 => Connectivity2::Four,
        8 => Connectivity2::Eight,
        _ => return Err(PyValueError::new_err("connectivity must be 4 or 8")),
    };
    let img = image_from(&mask)?;
    let labeling = py.detach(|| label::label_2d(&img, conn)).map_err(image_err)?;
    let count = labeling.count;
    Ok((label_array(py, labeling.labels)?, count))
}

/// Marker-controlled watershed of a 2D image.
///
/// `markers` is a `uint32` array of seed labels (0 = unlabelled). `mask`, if
/// given, is a boolean array restricting flooding to `True` pixels.
#[pyfunction]
#[pyo3(signature = (image, markers, mask=None))]
fn watershed_2d<'py>(
    py: Python<'py>,
    image: &Bound<'py, PyAny>,
    markers: PyReadonlyArray2<'py, u32>,
    mask: Option<PyReadonlyArray2<'py, bool>>,
) -> PyResult<Bound<'py, PyArray2<u32>>> {
    let img = image_f64(image)?;
    let mk = image_from(&markers)?;
    let m = mask.as_ref().map(image_from).transpose()?;
    let out = py.detach(|| watershed::watershed(&img, &mk, m.as_ref())).map_err(image_err)?;
    label_array(py, out)
}

/// Converts an image to a 2D NumPy array with the same element type.
fn image_to_py<'py, T: Element>(
    py: Python<'py>,
    img: Image<T>,
) -> PyResult<Bound<'py, PyArray2<T>>> {
    let (w, h, _) = img.shape();
    let arr = Array2::from_shape_vec((h, w), img.into_data())
        .map_err(|e| PyValueError::new_err(e.to_string()))?;
    Ok(PyArray2::from_owned_array(py, arr))
}

/// Picks one channel of a decoded file. Only single-plane channels are accepted.
fn pick_channel<T>(mut channels: Vec<Image<T>>, channel: usize) -> PyResult<Image<T>> {
    let count = channels.len();
    if channel >= count {
        return Err(PyValueError::new_err(format!(
            "channel {channel} is out of range; the file has {count} channel(s)"
        )));
    }
    let img = channels.swap_remove(channel);
    if img.depth() != 1 {
        return Err(PyValueError::new_err(format!(
            "the file is a volume (depth {}); only 2D images are supported",
            img.depth()
        )));
    }
    Ok(img)
}

fn read_planes(path: &str) -> PyResult<Planes> {
    let source = FileSource::open(Path::new(path))
        .map_err(|e| PyIOError::new_err(format!("cannot open '{path}': {e}")))?;
    let ext = Path::new(path)
        .extension()
        .and_then(|e| e.to_str())
        .map(str::to_ascii_lowercase)
        .unwrap_or_default();
    let decoded = match ext.as_str() {
        "png" => png::read_png(&source),
        "tif" | "tiff" => tiff::read_tiff(&source),
        _ => return Err(PyValueError::new_err(format!("'{path}' is not a PNG or TIFF file"))),
    };
    decoded.map_err(|e| PyValueError::new_err(format!("cannot decode '{path}': {e}")))
}

/// Loads one channel of a PNG or TIFF file as a 2D NumPy array, keeping its sample
/// type (`uint8`, `uint16` or `float32`). Returns `(array, channel_count)`.
#[pyfunction]
#[pyo3(signature = (path, channel=0))]
fn load_image<'py>(
    py: Python<'py>,
    path: &str,
    channel: usize,
) -> PyResult<(Bound<'py, PyAny>, usize)> {
    let planes = py.detach(|| read_planes(path))?;
    let count = planes.channel_count();
    let array: Bound<'py, PyAny> = match planes {
        Planes::U8(v) => image_to_py(py, pick_channel(v, channel)?)?.into_any(),
        Planes::U16(v) => image_to_py(py, pick_channel(v, channel)?)?.into_any(),
        Planes::F32(v) => image_to_py(py, pick_channel(v, channel)?)?.into_any(),
    };
    Ok((array, count))
}

/// Per-object shape and, when `intensity` is given, intensity features.
///
/// `labels` is a 2D `uint32` label image (background 0). Returns a dict mapping
/// each column name to a list with one value per object, in label order.
#[pyfunction]
#[pyo3(signature = (labels, intensity=None))]
fn object_features<'py>(
    py: Python<'py>,
    labels: PyReadonlyArray2<'py, u32>,
    intensity: Option<&Bound<'py, PyAny>>,
) -> PyResult<Bound<'py, PyDict>> {
    let lab = image_from(&labels)?;
    let intensity_img = intensity.map(image_f64).transpose()?;
    let shape = py.detach(|| shape_features_2d(&lab)).map_err(image_err)?;

    let mut columns: Vec<(&str, Vec<f64>)> = vec![
        ("label", shape.iter().map(|s| f64::from(s.label)).collect()),
        ("area", shape.iter().map(|s| s.area as f64).collect()),
        ("perimeter", shape.iter().map(|s| s.perimeter).collect()),
        ("eccentricity", shape.iter().map(|s| s.eccentricity).collect()),
        ("solidity", shape.iter().map(|s| s.solidity).collect()),
    ];
    if let Some(img) = &intensity_img {
        let ints = py.detach(|| intensity_features(&lab, img)).map_err(image_err)?;
        if ints.len() != shape.len() {
            return Err(PyValueError::new_err("shape and intensity measurements disagree"));
        }
        columns.extend([
            ("mean", ints.iter().map(|i| i.mean).collect()),
            ("std", ints.iter().map(|i| i.std).collect()),
            ("min", ints.iter().map(|i| i.min).collect()),
            ("max", ints.iter().map(|i| i.max).collect()),
            ("q25", ints.iter().map(|i| i.q25).collect()),
            ("median", ints.iter().map(|i| i.median).collect()),
            ("q75", ints.iter().map(|i| i.q75).collect()),
            ("integrated", ints.iter().map(|i| i.integrated).collect()),
        ]);
    }
    let dict = PyDict::new(py);
    for (name, values) in columns {
        dict.set_item(name, values)?;
    }
    Ok(dict)
}

/// A pipeline loaded from a TOML definition. Build it with `Pipeline.from_toml(path)`.
#[pyclass(name = "Pipeline", module = "tpt_bio", frozen)]
pub struct PyPipeline {
    inner: TptPipeline,
}

#[pymethods]
impl PyPipeline {
    /// Reads and validates a TOML pipeline definition file.
    ///
    /// Raises OSError if the file cannot be read and ValueError if it is invalid.
    #[staticmethod]
    fn from_toml(path: &str) -> PyResult<Self> {
        let text = std::fs::read_to_string(path)
            .map_err(|e| PyIOError::new_err(format!("cannot read '{path}': {e}")))?;
        let inner =
            TptPipeline::from_toml(&text).map_err(|e| PyValueError::new_err(e.to_string()))?;
        Ok(Self { inner })
    }

    /// Runs every node and returns the outputs in definition order.
    ///
    /// Tables are float64 arrays, labels are uint64 arrays, reads are FastqRecord
    /// lists, BAM records are BamRecord lists, images and masks are 2D arrays, objects
    /// are `(labels, count)` tuples, frames are dicts of lists, and files are paths.
    ///
    /// Raises RuntimeError if a node fails.
    fn run<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyList>> {
        let inner = &self.inner;
        let outputs =
            py.detach(|| inner.run()).map_err(|e| PyRuntimeError::new_err(e.to_string()))?;
        let list = PyList::empty(py);
        for value in outputs {
            list.append(value_to_py(py, value)?)?;
        }
        Ok(list)
    }
}

fn value_to_py(py: Python<'_>, value: PipelineValue) -> PyResult<Bound<'_, PyAny>> {
    Ok(match value {
        PipelineValue::Table(m) => {
            let (rows, cols) = (m.rows(), m.cols());
            let data: Vec<f64> = (0..rows).flat_map(|r| m.row(r).iter().copied()).collect();
            let arr = Array2::from_shape_vec((rows, cols), data)
                .map_err(|e| PyValueError::new_err(e.to_string()))?;
            PyArray2::from_owned_array(py, arr).into_any()
        }
        PipelineValue::Labels(labels) => {
            let v: Vec<u64> = labels.into_iter().map(|x| x as u64).collect();
            PyArray1::from_owned_array(py, Array1::from(v)).into_any()
        }
        PipelineValue::Fastq(reads) => {
            let items = reads
                .iter()
                .map(|buf| Py::new(py, to_py_record(py, buf)))
                .collect::<PyResult<Vec<_>>>()?;
            PyList::new(py, items)?.into_any()
        }
        PipelineValue::Records(records) => {
            let items = records
                .records
                .iter()
                .map(|r| Py::new(py, to_py_bam_record(py, r)))
                .collect::<PyResult<Vec<_>>>()?;
            PyList::new(py, items)?.into_any()
        }
        PipelineValue::File(path) => path.into_pyobject(py)?.into_any(),
        PipelineValue::Image(img) => image_to_py(py, img)?.into_any(),
        PipelineValue::Mask(mask) => image_to_py(py, mask)?.into_any(),
        PipelineValue::Objects(objects) => {
            let count = objects.count;
            let labels = image_to_py(py, objects.labels)?;
            (labels, count).into_pyobject(py)?.into_any()
        }
        PipelineValue::Frame(frame) => {
            let dict = PyDict::new(py);
            for (i, name) in frame.columns.iter().enumerate() {
                let column: Vec<f64> = frame.rows.iter().map(|row| row[i]).collect();
                dict.set_item(name, column)?;
            }
            dict.into_any()
        }
    })
}

macro_rules! opaque_debug {
    ($($t:ty),*) => {$(
        impl std::fmt::Debug for $t {
            fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
                f.debug_struct(stringify!($t)).finish_non_exhaustive()
            }
        }
    )*};
}

opaque_debug!(PyFastqRecord, PyFastqReader, PyBamRecord, PyIndexedBam, PyPipeline);

#[pymodule]
fn tpt_bio(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add("__version__", env!("CARGO_PKG_VERSION"))?;
    m.add_function(wrap_pyfunction!(reverse_complement, m)?)?;
    m.add_function(wrap_pyfunction!(gc_content, m)?)?;
    m.add_function(wrap_pyfunction!(canonical_kmers, m)?)?;
    m.add_class::<PyFastqRecord>()?;
    m.add_class::<PyFastqReader>()?;
    m.add_class::<PyBamRecord>()?;
    m.add_class::<PyIndexedBam>()?;
    m.add_function(wrap_pyfunction!(build_bai, m)?)?;
    m.add_function(wrap_pyfunction!(otsu_threshold, m)?)?;
    m.add_function(wrap_pyfunction!(apply_threshold, m)?)?;
    m.add_function(wrap_pyfunction!(label_components, m)?)?;
    m.add_function(wrap_pyfunction!(watershed_2d, m)?)?;
    m.add_function(wrap_pyfunction!(load_image, m)?)?;
    m.add_function(wrap_pyfunction!(object_features, m)?)?;
    m.add_class::<PyPipeline>()?;
    Ok(())
}
