//! `PyO3` bindings for tpt-bio.
//!
//! The module is `tpt_bio`. The functions release the GIL while they work, as
//! the Phase 5 plan requires for heavy operations. This first slice exposes
//! sequence utilities only.

use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;

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

#[pymodule]
fn tpt_bio(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add("__version__", env!("CARGO_PKG_VERSION"))?;
    m.add_function(wrap_pyfunction!(reverse_complement, m)?)?;
    m.add_function(wrap_pyfunction!(gc_content, m)?)?;
    m.add_function(wrap_pyfunction!(canonical_kmers, m)?)?;
    Ok(())
}
