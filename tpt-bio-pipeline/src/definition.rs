//! TOML pipeline definitions.
//!
//! A definition lists nodes. Each node has a unique `name`, a `type` naming one of the
//! kinds in a registry, an optional `params` table, and an optional `inputs` list of
//! node names. Each name in `inputs` becomes an edge into the next input slot of the
//! node, in order. Edges are therefore written once, on the consuming node.
//!
//! ```
//! use tpt_bio_pipeline::Pipeline;
//!
//! let text = r#"
//! [[nodes]]
//! name = "data"
//! type = "table_source"
//! params = { values = [[1.0, 2.0], [3.0, 5.0], [4.0, 1.0]] }
//!
//! [[nodes]]
//! name = "scaled"
//! type = "normalize"
//! params = { method = "robust" }
//! inputs = ["data"]
//! "#;
//! let pipeline = Pipeline::from_toml(text).unwrap();
//! assert_eq!(pipeline.graph().len(), 2);
//! ```
//!
//! Every definition is checked when it is loaded. Errors name the node at fault.

use core::fmt;

use serde::{Deserialize, Serialize};
use toml::{Table, Value};

use crate::typed::{TypeTag, TypedDag, TypedError};
use crate::{DagError, NodeId};

/// A table of finite numbers, produced by sources and stats nodes.
pub const TABLE: TypeTag = TypeTag::new("table");
/// One integer label per row, produced by clustering nodes.
pub const LABELS: TypeTag = TypeTag::new("labels");
/// Reads from a FASTQ file.
pub const FASTQ: TypeTag = TypeTag::new("fastq");
/// Alignment records with their BAM header.
pub const RECORDS: TypeTag = TypeTag::new("records");
/// A file written by a sink node.
pub const FILE: TypeTag = TypeTag::new("file");
/// A 2D grey-level image, one `f64` per pixel.
pub const IMAGE: TypeTag = TypeTag::new("image");
/// A 2D boolean foreground mask.
pub const MASK: TypeTag = TypeTag::new("mask");
/// A labelled 2D image: background 0, objects `1..=count`.
pub const OBJECTS: TypeTag = TypeTag::new("objects");
/// Named numeric columns, one row per object.
pub const FRAME: TypeTag = TypeTag::new("frame");

/// A node type known to the definition format.
#[derive(Debug, Clone, Copy)]
pub struct NodeKind {
    /// The name used in the `type` field.
    pub name: &'static str,
    /// Input port types, in slot order.
    pub inputs: &'static [TypeTag],
    /// Output port type.
    pub output: TypeTag,
    /// The parameter keys this kind accepts. Any other key is an error.
    pub params: &'static [&'static str],
    /// Checks the parameter values. Returns a message describing the first problem.
    pub check: fn(&Table) -> Result<(), String>,
}

/// The node types the crate provides.
///
/// | type                 | inputs                 | output  | params                                          |
/// |----------------------|------------------------|---------|-------------------------------------------------|
/// | `table_source`       | none                   | table   | `values` (rectangular array of numbers)         |
/// | `normalize`          | table                  | table   | `method` = `"zscore"` or `"robust"`             |
/// | `pca`                | table                  | table   | `components` (integer >= 1)                     |
/// | `kmeans`             | table                  | labels  | `k` (>= 1), `seed` (>= 0), `max_iter` (>= 1)    |
/// | `frame_to_table`     | frame                  | table   | `columns` (optional list of column names)       |
/// | `read_fastq`        | none                   | fastq   | `path`                                          |
/// | `read_bam`           | none                   | records | `path`                                          |
/// | `write_records_csv`  | records                | file    | `path`                                          |
/// | `write_frame_csv`    | frame                  | file    | `path`                                          |
/// | `write_labels_csv`   | labels                 | file    | `path`                                          |
/// | `count_reads`        | fastq                  | table   | none (1x1 table holding the count)              |
/// | `count_records`      | records                | table   | none (1x1 table holding the count)              |
/// | `count_objects`      | objects                | table   | none (1x1 table holding the count)              |
/// | `load_image`         | none                   | image   | `path` (PNG or TIFF), `channel` (default 0)     |
/// | `gaussian`           | image                  | image   | `sigma` (> 0)                                   |
/// | `subtract_background`| image                  | image   | `radius` (>= 0)                                 |
/// | `otsu`               | image                  | mask    | `bins` (>= 2, default 256)                      |
/// | `label`              | mask                   | objects | `connectivity` = 4 or 8 (default 8)             |
/// | `filter_objects`     | objects                | objects | `min_area`, `max_area`, `remove_border`         |
/// | `log_seeds`          | image                  | objects | `sigma`, `min_distance`, `min_response`         |
/// | `watershed`          | image, objects, mask   | objects | none                                            |
/// | `measure`            | objects, image         | frame   | none (shape and intensity columns)              |
pub const BUILTIN_KINDS: &[NodeKind] = &[
    NodeKind {
        name: "table_source",
        inputs: &[],
        output: TABLE,
        params: &["values"],
        check: check_table_source,
    },
    NodeKind {
        name: "normalize",
        inputs: &[TABLE],
        output: TABLE,
        params: &["method"],
        check: check_normalize,
    },
    NodeKind {
        name: "pca",
        inputs: &[TABLE],
        output: TABLE,
        params: &["components"],
        check: check_pca,
    },
    NodeKind {
        name: "kmeans",
        inputs: &[TABLE],
        output: LABELS,
        params: &["k", "seed", "max_iter"],
        check: check_kmeans,
    },
    NodeKind {
        name: "frame_to_table",
        inputs: &[FRAME],
        output: TABLE,
        params: &["columns"],
        check: check_frame_to_table,
    },
    NodeKind {
        name: "read_fastq",
        inputs: &[],
        output: FASTQ,
        params: &["path"],
        check: check_path,
    },
    NodeKind {
        name: "read_bam",
        inputs: &[],
        output: RECORDS,
        params: &["path"],
        check: check_path,
    },
    NodeKind {
        name: "write_records_csv",
        inputs: &[RECORDS],
        output: FILE,
        params: &["path"],
        check: check_path,
    },
    NodeKind {
        name: "write_frame_csv",
        inputs: &[FRAME],
        output: FILE,
        params: &["path"],
        check: check_path,
    },
    NodeKind {
        name: "write_labels_csv",
        inputs: &[LABELS],
        output: FILE,
        params: &["path"],
        check: check_path,
    },
    NodeKind {
        name: "count_reads",
        inputs: &[FASTQ],
        output: TABLE,
        params: &[],
        check: check_none,
    },
    NodeKind {
        name: "count_records",
        inputs: &[RECORDS],
        output: TABLE,
        params: &[],
        check: check_none,
    },
    NodeKind {
        name: "count_objects",
        inputs: &[OBJECTS],
        output: TABLE,
        params: &[],
        check: check_none,
    },
    NodeKind {
        name: "load_image",
        inputs: &[],
        output: IMAGE,
        params: &["path", "channel"],
        check: check_load_image,
    },
    NodeKind {
        name: "gaussian",
        inputs: &[IMAGE],
        output: IMAGE,
        params: &["sigma"],
        check: check_gaussian,
    },
    NodeKind {
        name: "subtract_background",
        inputs: &[IMAGE],
        output: IMAGE,
        params: &["radius"],
        check: check_subtract_background,
    },
    NodeKind { name: "otsu", inputs: &[IMAGE], output: MASK, params: &["bins"], check: check_otsu },
    NodeKind {
        name: "label",
        inputs: &[MASK],
        output: OBJECTS,
        params: &["connectivity"],
        check: check_label,
    },
    NodeKind {
        name: "filter_objects",
        inputs: &[OBJECTS],
        output: OBJECTS,
        params: &["min_area", "max_area", "remove_border"],
        check: check_filter_objects,
    },
    NodeKind {
        name: "log_seeds",
        inputs: &[IMAGE],
        output: OBJECTS,
        params: &["sigma", "min_distance", "min_response"],
        check: check_log_seeds,
    },
    NodeKind {
        name: "watershed",
        inputs: &[IMAGE, OBJECTS, MASK],
        output: OBJECTS,
        params: &[],
        check: check_none,
    },
    NodeKind {
        name: "measure",
        inputs: &[OBJECTS, IMAGE],
        output: FRAME,
        params: &[],
        check: check_none,
    },
];

/// The built-in node kinds, for passing to [`Pipeline::from_toml_with`].
#[must_use]
pub const fn builtin_kinds() -> &'static [NodeKind] {
    BUILTIN_KINDS
}

/// One node as written in a definition file.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct Node {
    /// Unique name, used by other nodes to refer to this one.
    pub name: String,
    /// The node kind, written as `type` in TOML.
    #[serde(rename = "type")]
    pub kind: String,
    /// Kind-specific parameters.
    #[serde(default, skip_serializing_if = "Table::is_empty")]
    pub params: Table,
    /// Names of the nodes feeding this node's input slots, in order.
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub inputs: Vec<String>,
}

#[derive(Debug, Serialize, Deserialize)]
struct Document {
    #[serde(default)]
    nodes: Vec<Node>,
}

/// A validated pipeline: a typed graph of [`Node`] definitions.
#[derive(Debug, Clone)]
pub struct Pipeline {
    graph: TypedDag<Node>,
}

impl Pipeline {
    /// Parses and validates a definition using the built-in node kinds.
    ///
    /// # Errors
    /// Returns [`DefinitionError`] for malformed TOML, unknown kinds, bad parameters,
    /// unknown or duplicated inputs, type mismatches, missing inputs and cycles.
    pub fn from_toml(text: &str) -> Result<Self, DefinitionError> {
        Self::from_toml_with(text, BUILTIN_KINDS)
    }

    /// Parses and validates a definition against the given kinds.
    ///
    /// # Errors
    /// As for [`Pipeline::from_toml`], with `kinds` as the registry.
    pub fn from_toml_with(text: &str, kinds: &[NodeKind]) -> Result<Self, DefinitionError> {
        let doc: Document =
            toml::from_str(text).map_err(|e| DefinitionError::Parse(e.to_string()))?;
        build(&doc.nodes, kinds)
    }

    /// Serialises the pipeline back to TOML. Nodes appear in definition order.
    ///
    /// # Errors
    /// Returns the serialiser's error if a parameter cannot be written as TOML.
    pub fn to_toml(&self) -> Result<String, toml::ser::Error> {
        let nodes = self.graph.node_ids().filter_map(|id| self.graph.payload(id).cloned());
        let doc = Document { nodes: nodes.collect() };
        toml::to_string(&doc)
    }

    /// The validated graph. Payloads are the [`Node`] definitions.
    #[must_use]
    pub const fn graph(&self) -> &TypedDag<Node> {
        &self.graph
    }
}

fn build(nodes: &[Node], kinds: &[NodeKind]) -> Result<Pipeline, DefinitionError> {
    let mut graph = TypedDag::new();
    // Pass 1: register every node, checking its kind, parameters and name.
    for node in nodes {
        let kind = kinds.iter().find(|k| k.name == node.kind).ok_or_else(|| {
            DefinitionError::UnknownKind { node: node.name.clone(), kind: node.kind.clone() }
        })?;
        check_keys(&node.params, kind.params)
            .and_then(|()| (kind.check)(&node.params))
            .map_err(|message| DefinitionError::BadParam { node: node.name.clone(), message })?;
        graph.add_node(&node.name, kind.inputs, kind.output, node.clone()).map_err(
            |e| match e {
                TypedError::InvalidName(_) => DefinitionError::InvalidName(node.name.clone()),
                TypedError::DuplicateName(_) => DefinitionError::DuplicateName(node.name.clone()),
                other => DefinitionError::Other(other.to_string()),
            },
        )?;
    }

    // Pass 2: wire the inputs, in order, into each node's slots.
    for node in nodes {
        let to = graph.node_by_name(&node.name).ok_or_else(|| {
            DefinitionError::Other(format!("node '{}' was not registered", node.name))
        })?;
        let declared = graph.inputs(to).map_or(0, <[TypeTag]>::len);
        if node.inputs.len() > declared {
            return Err(DefinitionError::TooManyInputs {
                node: node.name.clone(),
                declared,
                given: node.inputs.len(),
            });
        }
        for (slot, input) in node.inputs.iter().enumerate() {
            if node.inputs[..slot].contains(input) {
                return Err(DefinitionError::DuplicateInput {
                    node: node.name.clone(),
                    input: input.clone(),
                });
            }
            let from = graph.node_by_name(input).ok_or_else(|| DefinitionError::UnknownInput {
                node: node.name.clone(),
                input: input.clone(),
            })?;
            graph.add_edge(from, to).map_err(|e| edge_error(&graph, e, &node.name, input))?;
        }
    }

    graph.validate().map_err(|e| match e {
        TypedError::MissingInput { node, slot, expected } => DefinitionError::MissingInput {
            node: name_of(&graph, node),
            slot,
            expected: expected.name(),
        },
        other => DefinitionError::Other(other.to_string()),
    })?;
    Ok(Pipeline { graph })
}

fn edge_error(graph: &TypedDag<Node>, e: TypedError, node: &str, input: &str) -> DefinitionError {
    match e {
        TypedError::TypeMismatch { from, expected, found, .. } => DefinitionError::TypeMismatch {
            node: node.to_string(),
            input: input.to_string(),
            producer: name_of(graph, from),
            expected: expected.name(),
            found: found.name(),
        },
        TypedError::Dag(DagError::Cycle { from, to }) => {
            DefinitionError::Cycle { from: name_of(graph, from), to: name_of(graph, to) }
        }
        other => DefinitionError::Other(other.to_string()),
    }
}

fn name_of(graph: &TypedDag<Node>, id: NodeId) -> String {
    graph.name(id).map_or_else(|| format!("#{}", id.index()), ToString::to_string)
}

/// Errors from loading a pipeline definition. Each variant names the node at fault.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum DefinitionError {
    /// The text is not valid TOML or does not match the schema.
    Parse(String),
    /// A node name is empty or contains characters outside `[A-Za-z0-9_.-]`.
    InvalidName(String),
    /// Two nodes share a name.
    DuplicateName(String),
    /// A node's `type` is not in the registry.
    UnknownKind {
        /// The node with the unknown type.
        node: String,
        /// The unknown type name.
        kind: String,
    },
    /// A parameter is missing, has the wrong type, is unknown, or is out of range.
    BadParam {
        /// The node with the bad parameter.
        node: String,
        /// What is wrong with it.
        message: String,
    },
    /// An input names a node that is not in the definition.
    UnknownInput {
        /// The node whose input is unknown.
        node: String,
        /// The missing node name.
        input: String,
    },
    /// The same producer is listed twice in one node's inputs.
    DuplicateInput {
        /// The node with the repeated input.
        node: String,
        /// The repeated producer name.
        input: String,
    },
    /// A node lists more inputs than its kind declares.
    TooManyInputs {
        /// The node with too many inputs.
        node: String,
        /// Number of slots the kind declares.
        declared: usize,
        /// Number of inputs listed.
        given: usize,
    },
    /// A producer's output type does not match the consumer's input slot.
    TypeMismatch {
        /// The consuming node.
        node: String,
        /// The producer named in `inputs`.
        input: String,
        /// The producer's name, as stored in the graph.
        producer: String,
        /// The type the consumer's slot expects.
        expected: &'static str,
        /// The type the producer outputs.
        found: &'static str,
    },
    /// A declared input slot has no input listed.
    MissingInput {
        /// The node with the unconnected slot.
        node: String,
        /// Zero-based slot index.
        slot: usize,
        /// The type the slot expects.
        expected: &'static str,
    },
    /// An edge would close a cycle.
    Cycle {
        /// Producer of the rejected edge.
        from: String,
        /// Consumer of the rejected edge.
        to: String,
    },
    /// An error from the graph layer that has no specific variant.
    Other(String),
}

impl fmt::Display for DefinitionError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Parse(msg) => write!(f, "invalid pipeline definition: {msg}"),
            Self::InvalidName(n) => write!(f, "node name {n:?} is not filesystem-safe"),
            Self::DuplicateName(n) => write!(f, "node '{n}' is defined more than once"),
            Self::UnknownKind { node, kind } => {
                write!(f, "node '{node}' has unknown type '{kind}'")
            }
            Self::BadParam { node, message } => write!(f, "node '{node}': {message}"),
            Self::UnknownInput { node, input } => {
                write!(f, "node '{node}' takes input '{input}', which is not defined")
            }
            Self::DuplicateInput { node, input } => {
                write!(f, "node '{node}' lists input '{input}' more than once")
            }
            Self::TooManyInputs { node, declared, given } => {
                write!(f, "node '{node}' has {given} inputs but its type declares {declared}")
            }
            Self::TypeMismatch { node, input, producer, expected, found } => write!(
                f,
                "node '{node}' expects {expected} from input '{input}', \
                 but '{producer}' produces {found}"
            ),
            Self::MissingInput { node, slot, expected } => {
                write!(f, "node '{node}' input {slot} of type {expected} is not connected")
            }
            Self::Cycle { from, to } => {
                write!(f, "input from '{from}' to '{to}' would create a cycle")
            }
            Self::Other(msg) => f.write_str(msg),
        }
    }
}

impl core::error::Error for DefinitionError {}

fn check_keys(params: &Table, allowed: &[&str]) -> Result<(), String> {
    for key in params.keys() {
        if !allowed.contains(&key.as_str()) {
            return Err(format!("unknown parameter '{key}'"));
        }
    }
    Ok(())
}

pub(crate) fn int_param(params: &Table, key: &str, min: i64) -> Result<Option<i64>, String> {
    match params.get(key) {
        None => Ok(None),
        Some(Value::Integer(v)) if *v >= min => Ok(Some(*v)),
        Some(Value::Integer(v)) => Err(format!("parameter '{key}' must be >= {min}, got {v}")),
        Some(other) => {
            Err(format!("parameter '{key}' must be an integer, found {}", other.type_str()))
        }
    }
}

/// Reads the `values` parameter as a rectangular, non-empty matrix of finite numbers.
///
/// # Errors
/// Returns a message describing the first problem found.
pub fn parse_values(params: &Table) -> Result<Vec<Vec<f64>>, String> {
    let Some(outer) = params.get("values") else {
        return Err("missing required parameter 'values'".to_string());
    };
    let Value::Array(rows) = outer else {
        return Err("parameter 'values' must be an array of rows".to_string());
    };
    if rows.is_empty() {
        return Err("parameter 'values' has no rows".to_string());
    }
    let mut out: Vec<Vec<f64>> = Vec::with_capacity(rows.len());
    for (r, row) in rows.iter().enumerate() {
        let Value::Array(cells) = row else {
            return Err(format!("row {r} of 'values' is not an array"));
        };
        let mut parsed = Vec::with_capacity(cells.len());
        for cell in cells {
            let v = match cell {
                Value::Float(v) => *v,
                Value::Integer(v) => *v as f64,
                _ => return Err(format!("row {r} of 'values' has a non-numeric cell")),
            };
            if !v.is_finite() {
                return Err(format!("row {r} of 'values' has a non-finite cell"));
            }
            parsed.push(v);
        }
        if parsed.is_empty() {
            return Err(format!("row {r} of 'values' is empty"));
        }
        if out.first().is_some_and(|first| first.len() != parsed.len()) {
            return Err(format!("row {r} of 'values' has a different length from row 0"));
        }
        out.push(parsed);
    }
    Ok(out)
}

fn check_table_source(params: &Table) -> Result<(), String> {
    parse_values(params).map(|_| ())
}

fn check_normalize(params: &Table) -> Result<(), String> {
    match params.get("method") {
        None => Ok(()),
        Some(Value::String(m)) if m == "zscore" || m == "robust" => Ok(()),
        Some(Value::String(m)) => Err(format!("unknown method '{m}'")),
        Some(other) => {
            Err(format!("parameter 'method' must be a string, found {}", other.type_str()))
        }
    }
}

fn check_pca(params: &Table) -> Result<(), String> {
    match int_param(params, "components", 1)? {
        Some(_) => Ok(()),
        None => Err("missing required parameter 'components'".to_string()),
    }
}

fn check_kmeans(params: &Table) -> Result<(), String> {
    if int_param(params, "k", 1)?.is_none() {
        return Err("missing required parameter 'k'".to_string());
    }
    int_param(params, "seed", 0)?;
    int_param(params, "max_iter", 1)?;
    Ok(())
}

pub(crate) fn str_param<'a>(params: &'a Table, key: &str) -> Result<Option<&'a str>, String> {
    match params.get(key) {
        None => Ok(None),
        Some(Value::String(s)) => Ok(Some(s.as_str())),
        Some(other) => {
            Err(format!("parameter '{key}' must be a string, found {}", other.type_str()))
        }
    }
}

pub(crate) fn float_param(params: &Table, key: &str, min: f64) -> Result<Option<f64>, String> {
    let v = match params.get(key) {
        None => return Ok(None),
        Some(Value::Float(v)) => *v,
        Some(Value::Integer(v)) => *v as f64,
        Some(other) => {
            return Err(format!("parameter '{key}' must be a number, found {}", other.type_str()));
        }
    };
    if !v.is_finite() || v < min {
        return Err(format!("parameter '{key}' must be a finite number >= {min}, got {v}"));
    }
    Ok(Some(v))
}

pub(crate) fn bool_param(params: &Table, key: &str) -> Result<Option<bool>, String> {
    match params.get(key) {
        None => Ok(None),
        Some(Value::Boolean(b)) => Ok(Some(*b)),
        Some(other) => {
            Err(format!("parameter '{key}' must be a boolean, found {}", other.type_str()))
        }
    }
}

pub(crate) fn str_list_param(params: &Table, key: &str) -> Result<Option<Vec<String>>, String> {
    match params.get(key) {
        None => Ok(None),
        Some(Value::Array(items)) => items
            .iter()
            .map(|v| match v {
                Value::String(s) => Ok(s.clone()),
                _ => Err(format!("parameter '{key}' must be a list of strings")),
            })
            .collect::<Result<Vec<_>, _>>()
            .map(Some),
        Some(_) => Err(format!("parameter '{key}' must be a list of strings")),
    }
}

fn check_none(params: &Table) -> Result<(), String> {
    check_keys(params, &[])
}

fn check_path(params: &Table) -> Result<(), String> {
    match str_param(params, "path")? {
        Some(p) if !p.is_empty() => Ok(()),
        Some(_) => Err("parameter 'path' is empty".to_string()),
        None => Err("missing required parameter 'path'".to_string()),
    }
}

fn check_frame_to_table(params: &Table) -> Result<(), String> {
    str_list_param(params, "columns").map(|_| ())
}

fn check_load_image(params: &Table) -> Result<(), String> {
    check_path(params)?;
    int_param(params, "channel", 0).map(|_| ())
}

fn check_gaussian(params: &Table) -> Result<(), String> {
    match float_param(params, "sigma", 0.0)? {
        Some(s) if s > 0.0 => Ok(()),
        Some(_) => Err("parameter 'sigma' must be greater than 0".to_string()),
        None => Err("missing required parameter 'sigma'".to_string()),
    }
}

fn check_subtract_background(params: &Table) -> Result<(), String> {
    match int_param(params, "radius", 0)? {
        Some(_) => Ok(()),
        None => Err("missing required parameter 'radius'".to_string()),
    }
}

fn check_otsu(params: &Table) -> Result<(), String> {
    int_param(params, "bins", 2).map(|_| ())
}

fn check_label(params: &Table) -> Result<(), String> {
    match int_param(params, "connectivity", 0)? {
        None | Some(4 | 8) => Ok(()),
        Some(other) => Err(format!("parameter 'connectivity' must be 4 or 8, got {other}")),
    }
}

fn check_filter_objects(params: &Table) -> Result<(), String> {
    int_param(params, "min_area", 0)?;
    int_param(params, "max_area", 0)?;
    bool_param(params, "remove_border").map(|_| ())
}

fn check_log_seeds(params: &Table) -> Result<(), String> {
    if let Some(s) = float_param(params, "sigma", 0.0)?
        && s <= 0.0
    {
        return Err("parameter 'sigma' must be greater than 0".to_string());
    }
    int_param(params, "min_distance", 0)?;
    float_param(params, "min_response", f64::NEG_INFINITY).map(|_| ())
}

#[cfg(test)]
mod tests {
    use super::*;

    const GOOD: &str = r#"
[[nodes]]
name = "data"
type = "table_source"
params = { values = [[1.0, 2.0], [3.0, 5.0], [4.0, 1.0]] }

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
params = { k = 2, seed = 7 }
inputs = ["reduced"]
"#;

    fn err(text: &str) -> DefinitionError {
        Pipeline::from_toml(text).unwrap_err()
    }

    #[test]
    fn loads_valid_definition() {
        let p = Pipeline::from_toml(GOOD).unwrap();
        assert_eq!(p.graph().len(), 4);
        let clusters = p.graph().node_by_name("clusters").unwrap();
        assert_eq!(p.graph().output(clusters), Some(LABELS));
    }

    #[test]
    fn round_trip_preserves_definition() {
        let first = Pipeline::from_toml(GOOD).unwrap();
        let text = first.to_toml().unwrap();
        let second = Pipeline::from_toml(&text).unwrap();
        assert_eq!(text, second.to_toml().unwrap());
        assert_eq!(second.graph().len(), 4);
        // Definition order is kept.
        let names: Vec<&str> =
            second.graph().node_ids().map(|id| second.graph().name(id).unwrap()).collect();
        assert_eq!(names, ["data", "scaled", "reduced", "clusters"]);
    }

    #[test]
    fn unknown_type_names_the_node() {
        let text = GOOD.replace("type = \"pca\"", "type = \"svd\"");
        let e = err(&text);
        assert_eq!(e, DefinitionError::UnknownKind { node: "reduced".into(), kind: "svd".into() });
        assert!(e.to_string().contains("'reduced'"));
    }

    #[test]
    fn missing_input_names_the_node() {
        let text = GOOD.replace("inputs = [\"scaled\"]", "");
        let e = err(&text);
        assert_eq!(
            e,
            DefinitionError::MissingInput { node: "reduced".into(), slot: 0, expected: "table" }
        );
        assert!(e.to_string().contains("'reduced'"));
    }

    #[test]
    fn unknown_input_name_is_reported() {
        let text = GOOD.replace("inputs = [\"scaled\"]", "inputs = [\"nope\"]");
        assert_eq!(
            err(&text),
            DefinitionError::UnknownInput { node: "reduced".into(), input: "nope".into() }
        );
    }

    #[test]
    fn type_mismatch_names_both_nodes() {
        // kmeans outputs labels, and pca expects a table.
        let text = GOOD.replace("inputs = [\"reduced\"]\n", "")
            + "\n[[nodes]]\nname = \"bad\"\ntype = \"pca\"\nparams = { components = 1 }\ninputs = [\"clusters\"]\n";
        let e = err(&text);
        assert_eq!(
            e,
            DefinitionError::TypeMismatch {
                node: "bad".into(),
                input: "clusters".into(),
                producer: "clusters".into(),
                expected: "table",
                found: "labels",
            }
        );
    }

    #[test]
    fn cycle_is_reported() {
        let text = r#"
[[nodes]]
name = "a"
type = "normalize"
inputs = ["b"]

[[nodes]]
name = "b"
type = "normalize"
inputs = ["a"]
"#;
        assert!(matches!(err(text), DefinitionError::Cycle { .. }));
    }

    #[test]
    fn self_loop_is_a_cycle() {
        let text = "[[nodes]]\nname = \"a\"\ntype = \"normalize\"\ninputs = [\"a\"]\n";
        assert!(matches!(err(text), DefinitionError::Cycle { .. }));
    }

    #[test]
    fn bad_parameters_name_the_node() {
        let text = GOOD.replace("components = 1", "components = 0");
        let e = err(&text);
        assert!(matches!(&e, DefinitionError::BadParam { node, .. } if node == "reduced"));

        let text = GOOD.replace("method = \"zscore\"", "colour = \"red\"");
        assert!(err(&text).to_string().contains("unknown parameter 'colour'"));

        let text = GOOD.replace("k = 2, seed = 7", "seed = 7");
        assert!(err(&text).to_string().contains("'clusters'"));
    }

    #[test]
    fn ragged_values_are_rejected() {
        let text = GOOD.replace("[[1.0, 2.0], [3.0, 5.0], [4.0, 1.0]]", "[[1.0, 2.0], [3.0]]");
        assert!(matches!(err(&text), DefinitionError::BadParam { node, .. } if node == "data"));
    }

    #[test]
    fn duplicate_names_and_inputs_are_rejected() {
        let dup = GOOD.replace("name = \"scaled\"", "name = \"data\"");
        assert_eq!(err(&dup), DefinitionError::DuplicateName("data".into()));

        let source = "[[nodes]]
name = \"a\"
type = \"table_source\"
params = { values = [[1.0]] }
";
        // pca has one slot, so listing the same producer twice is a slot-count error.
        let twice_one_slot = format!(
            "{source}[[nodes]]
name = \"b\"
type = \"pca\"
params = {{ components = 1 }}
inputs = [\"a\", \"a\"]
"
        );
        assert!(
            matches!(err(&twice_one_slot), DefinitionError::TooManyInputs { node, .. } if node == "b")
        );
        // A two-slot kind reaches the duplicate check.
        let mut kinds = BUILTIN_KINDS.to_vec();
        kinds.push(NodeKind {
            name: "join",
            inputs: &[TABLE, TABLE],
            output: TABLE,
            params: &[],
            check: |_| Ok(()),
        });
        let twice = format!(
            "{source}[[nodes]]
name = \"b\"
type = \"join\"
inputs = [\"a\", \"a\"]
"
        );
        let e = Pipeline::from_toml_with(&twice, &kinds).unwrap_err();
        assert!(matches!(e, DefinitionError::DuplicateInput { node, .. } if node == "b"));
    }

    #[test]
    fn invalid_syntax_is_a_parse_error() {
        assert!(matches!(err("[[nodes]\nname = 1"), DefinitionError::Parse(_)));
    }
}
