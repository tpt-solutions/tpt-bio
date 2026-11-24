//! Runs a validated [`Pipeline`] whose nodes are the built-in statistics kinds.
//!
//! Each kind is evaluated by a function in `tpt-bio-core::stats`:
//!
//! * `table_source` builds a [`Matrix`] from its `values` parameter.
//! * `normalize` applies [`normalize_columns`] with z-score or robust scaling.
//! * `pca` fits a [`Pca`] and projects the table onto its components.
//! * `kmeans` clusters the rows with [`kmeans`] and returns one label per row.
//!
//! Values flow between nodes through the [`crate::Dag`] executor, so a node receives its
//! predecessors' outputs in the order of its input slots.

use std::fmt;

use tpt_bio_core::stats::normalize::{Method, normalize_columns};
use tpt_bio_core::stats::pca::Pca;
use tpt_bio_core::stats::{Matrix, cluster::kmeans};

use tpt_bio_seq::{FastqBuf, Header, Record};
use tpt_bio_vision::image::Image;
use tpt_bio_vision::label::Labeling;

use crate::definition::{
    FASTQ, FILE, FRAME, IMAGE, LABELS, MASK, Node, OBJECTS, Pipeline, RECORDS, TABLE, parse_values,
};
use crate::{ExecError, NodeId};

/// Alignment records read from a BAM file, with the header they came with.
#[derive(Debug, Clone, PartialEq)]
pub struct Records {
    /// The BAM header, used for reference names when writing.
    pub header: Header,
    /// The records in file order.
    pub records: Vec<Record>,
}

/// Named numeric columns with one row per object. Rows may be empty.
#[derive(Debug, Clone, PartialEq)]
pub struct Frame {
    /// Column names, in order.
    pub columns: Vec<String>,
    /// Row-major values. Every row has `columns.len()` entries.
    pub rows: Vec<Vec<f64>>,
}

/// A value produced by a node.
#[derive(Debug, Clone, PartialEq)]
pub enum Value {
    /// A numeric table, one row per sample.
    Table(Matrix),
    /// One cluster label per row.
    Labels(Vec<usize>),
    /// Reads parsed from a FASTQ file.
    Fastq(Vec<FastqBuf>),
    /// Alignment records from a BAM file.
    Records(Records),
    /// A path to a file written by a sink node.
    File(String),
    /// A 2D grey-level image.
    Image(Image<f64>),
    /// A 2D boolean foreground mask.
    Mask(Image<bool>),
    /// A labelled 2D image with its object count.
    Objects(Labeling),
    /// Per-object features.
    Frame(Frame),
}

impl Value {
    /// The port type this value satisfies.
    #[must_use]
    pub const fn type_tag(&self) -> crate::TypeTag {
        match self {
            Self::Table(_) => TABLE,
            Self::Labels(_) => LABELS,
            Self::Fastq(_) => FASTQ,
            Self::Records(_) => RECORDS,
            Self::File(_) => FILE,
            Self::Image(_) => IMAGE,
            Self::Mask(_) => MASK,
            Self::Objects(_) => OBJECTS,
            Self::Frame(_) => FRAME,
        }
    }
}

/// Errors from running a pipeline.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum RunError {
    /// A node failed. `node` is its name.
    Node {
        /// Name of the failing node.
        node: String,
        /// What went wrong.
        message: String,
    },
    /// The executor itself failed, such as on an invalid graph.
    Graph(String),
}

impl fmt::Display for RunError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Node { node, message } => write!(f, "node '{node}' failed: {message}"),
            Self::Graph(msg) => write!(f, "pipeline did not run: {msg}"),
        }
    }
}

impl std::error::Error for RunError {}

impl Pipeline {
    /// Runs every node and returns the outputs in definition order.
    ///
    /// Nodes run one after another in dependency order.
    ///
    /// # Errors
    /// Returns [`RunError::Node`] when a node's computation fails, such as a PCA with
    /// more components than columns, or [`RunError::Graph`] for executor failures.
    pub fn run(&self) -> Result<Vec<Value>, RunError> {
        let cancel = std::sync::atomic::AtomicBool::new(false);
        let graph = self.graph();
        graph
            .dag()
            .execute(&cancel, |id: NodeId, node: &Node, inputs: &[&Value]| {
                evaluate(node, inputs).map_err(|message| RunError::Node {
                    node: graph.name(id).unwrap_or_default().to_string(),
                    message,
                })
            })
            .map_err(|e| match e {
                ExecError::Node { error, .. } => error,
                other => RunError::Graph(other.to_string()),
            })
    }
}

fn evaluate(node: &Node, inputs: &[&Value]) -> Result<Value, String> {
    if let Some(result) = crate::seq_nodes::evaluate(node, inputs) {
        return result;
    }
    if let Some(result) = crate::vision_nodes::evaluate(node, inputs) {
        return result;
    }
    match node.kind.as_str() {
        "frame_to_table" => frame_to_table(node, inputs),
        "table_source" => {
            let rows = parse_values(&node.params)?;
            Matrix::from_rows(&rows).map(Value::Table).map_err(|e| e.to_string())
        }
        "normalize" => {
            let data = table_input(inputs)?;
            let method = match node.params.get("method").and_then(toml::Value::as_str) {
                Some("robust") => Method::Robust,
                _ => Method::ZScore,
            };
            normalize_columns(data, method).map(Value::Table).map_err(|e| e.to_string())
        }
        "pca" => {
            let data = table_input(inputs)?;
            let components = int_param(node, "components")?;
            let pca = Pca::fit(data, components).map_err(|e| e.to_string())?;
            pca.transform(data).map(Value::Table).map_err(|e| e.to_string())
        }
        "kmeans" => {
            let data = table_input(inputs)?;
            let k = int_param(node, "k")?;
            let seed = u64::try_from(int_param_or(node, "seed", 0)?)
                .map_err(|_| "seed must be non-negative".to_string())?;
            let max_iter = usize::try_from(int_param_or(node, "max_iter", 100)?)
                .map_err(|_| "max_iter must be positive".to_string())?;
            let fit = kmeans(data, k, seed, max_iter).map_err(|e| e.to_string())?;
            Ok(Value::Labels(fit.labels))
        }
        other => Err(format!("kind '{other}' has no executor")),
    }
}

fn frame_to_table(node: &Node, inputs: &[&Value]) -> Result<Value, String> {
    let frame = match inputs.first() {
        Some(Value::Frame(f)) => f,
        Some(other) => return Err(format!("expected a frame input, got {}", other.type_tag())),
        None => return Err("missing frame input".to_string()),
    };
    let wanted = crate::definition::str_list_param(&node.params, "columns")?;
    let picked: Vec<usize> = match &wanted {
        None => (0..frame.columns.len()).collect(),
        Some(names) => names
            .iter()
            .map(|name| {
                frame
                    .columns
                    .iter()
                    .position(|c| c == name)
                    .ok_or_else(|| format!("frame has no column '{name}'"))
            })
            .collect::<Result<_, _>>()?,
    };
    let rows: Vec<Vec<f64>> =
        frame.rows.iter().map(|row| picked.iter().map(|&c| row[c]).collect()).collect();
    Matrix::from_rows(&rows).map(Value::Table).map_err(|e| e.to_string())
}

pub(crate) fn table_input<'a>(inputs: &[&'a Value]) -> Result<&'a Matrix, String> {
    match inputs.first() {
        Some(Value::Table(m)) => Ok(m),
        Some(other) => Err(format!("expected a table input, got {}", other.type_tag())),
        None => Err("missing table input".to_string()),
    }
}

pub(crate) fn int_param(node: &Node, key: &str) -> Result<usize, String> {
    let v = int_param_or(node, key, -1)?;
    usize::try_from(v).map_err(|_| format!("parameter '{key}' is required"))
}

pub(crate) fn int_param_or(node: &Node, key: &str, default: i64) -> Result<i64, String> {
    match node.params.get(key) {
        None => Ok(default),
        Some(toml::Value::Integer(v)) => Ok(*v),
        Some(_) => Err(format!("parameter '{key}' must be an integer")),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const TEXT: &str = r#"
[[nodes]]
name = "data"
type = "table_source"
params = { values = [[0.0, 0.0], [0.0, 1.0], [10.0, 10.0], [10.0, 11.0]] }

[[nodes]]
name = "scaled"
type = "normalize"
params = { method = "zscore" }
inputs = ["data"]

[[nodes]]
name = "reduced"
type = "pca"
params = { components = 1 }
inputs = ["scaled"]

[[nodes]]
name = "clusters"
type = "kmeans"
params = { k = 2, seed = 1 }
inputs = ["reduced"]
"#;

    #[test]
    fn runs_source_normalize_pca_kmeans() {
        let pipeline = Pipeline::from_toml(TEXT).unwrap();
        let out = pipeline.run().unwrap();
        assert_eq!(out.len(), 4);
        assert!(matches!(out[0], Value::Table(_)));
        assert!(matches!(out[3], Value::Labels(_)));

        // Z-scored columns have mean zero.
        let Value::Table(scaled) = &out[1] else { panic!("expected table") };
        let mean0: f64 = (0..4).map(|r| scaled.get(r, 0)).sum::<f64>() / 4.0;
        assert!(mean0.abs() < 1e-12, "mean={mean0}");

        // The two pairs of rows form the two clusters.
        let Value::Labels(labels) = &out[3] else { panic!("expected labels") };
        assert_eq!(labels[0], labels[1]);
        assert_eq!(labels[2], labels[3]);
        assert_ne!(labels[0], labels[2]);
    }

    #[test]
    fn node_failure_names_the_node() {
        // Asking PCA for more components than the table has columns fails at run time.
        let text = TEXT.replace("components = 1", "components = 5");
        let pipeline = Pipeline::from_toml(&text).unwrap();
        match pipeline.run() {
            Err(RunError::Node { node, .. }) => assert_eq!(node, "reduced"),
            other => panic!("expected a node error, got {other:?}"),
        }
    }
}
