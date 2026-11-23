//! Typed ports for [`Dag`] nodes.
//!
//! Each node declares an ordered list of input port types and one output port type.
//! An edge `from -> to` fills the next free input slot of `to`, and its type must
//! equal that slot's declared type. Mismatches are rejected when the edge is added,
//! so a [`TypedDag`] cannot hold a wiring error. Declared inputs that are never
//! connected are reported by [`TypedDag::validate`].
//!
//! Node names are also kept here. They identify each node's checkpoint file, so
//! they must be filesystem-safe.

use std::fmt;

use crate::{Dag, DagError, NodeId};

/// A named port type. Two ports are compatible when their names are equal.
///
/// Names are `&'static str` so that type tags can be declared as constants:
///
/// ```
/// use tpt_bio_pipeline::TypeTag;
///
/// const READS: TypeTag = TypeTag::new("seq.reads");
/// assert_eq!(READS.name(), "seq.reads");
/// ```
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct TypeTag(&'static str);

impl TypeTag {
    /// A type tag with the given name.
    #[must_use]
    pub const fn new(name: &'static str) -> Self {
        Self(name)
    }

    /// The name this tag was created with.
    #[must_use]
    pub const fn name(self) -> &'static str {
        self.0
    }
}

impl fmt::Display for TypeTag {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(self.0)
    }
}

/// Errors from building a [`TypedDag`].
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum TypedError {
    /// An underlying graph error, such as a cycle or an unknown node.
    Dag(DagError),
    /// The node name is empty, starts with `.`, or has characters outside `[A-Za-z0-9_.-]`.
    InvalidName(String),
    /// Another node already has this name.
    DuplicateName(String),
    /// The target node has no declared input slot left for this edge.
    NoInputSlot {
        /// The node that received the extra edge.
        node: NodeId,
    },
    /// The source's output type does not match the target's input slot.
    TypeMismatch {
        /// Source of the rejected edge.
        from: NodeId,
        /// Target of the rejected edge.
        to: NodeId,
        /// Input slot of `to` that the edge would fill.
        slot: usize,
        /// Type the slot declares.
        expected: TypeTag,
        /// Type `from` produces.
        found: TypeTag,
    },
    /// A declared input slot has no incoming edge.
    MissingInput {
        /// The node with an unconnected input.
        node: NodeId,
        /// The unconnected slot.
        slot: usize,
        /// The type that slot expects.
        expected: TypeTag,
    },
}

impl fmt::Display for TypedError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Dag(e) => e.fmt(f),
            Self::InvalidName(n) => write!(f, "invalid node name {n:?}"),
            Self::DuplicateName(n) => write!(f, "duplicate node name {n:?}"),
            Self::NoInputSlot { node } => {
                write!(f, "node {} declares no further inputs", node.index())
            }
            Self::TypeMismatch { from, to, slot, expected, found } => write!(
                f,
                "edge {} -> {} fills input {slot} of type {expected}, but produces {found}",
                from.index(),
                to.index()
            ),
            Self::MissingInput { node, slot, expected } => {
                write!(f, "node {} input {slot} of type {expected} is not connected", node.index())
            }
        }
    }
}

impl std::error::Error for TypedError {}

impl From<DagError> for TypedError {
    fn from(e: DagError) -> Self {
        Self::Dag(e)
    }
}

/// Whether `name` is usable as a node name and as a checkpoint file stem.
pub(crate) fn is_valid_name(name: &str) -> bool {
    !name.is_empty()
        && !name.starts_with('.')
        && name.bytes().all(|b| b.is_ascii_alphanumeric() || matches!(b, b'_' | b'-' | b'.'))
}

#[derive(Debug, Clone)]
struct Ports {
    name: String,
    inputs: Vec<TypeTag>,
    output: TypeTag,
}

/// A [`Dag`] whose nodes declare typed ports and unique names.
#[derive(Debug, Clone)]
pub struct TypedDag<T> {
    dag: Dag<T>,
    ports: Vec<Ports>,
}

impl<T> Default for TypedDag<T> {
    fn default() -> Self {
        Self { dag: Dag::new(), ports: Vec::new() }
    }
}

impl<T> TypedDag<T> {
    /// An empty graph.
    #[must_use]
    pub fn new() -> Self {
        Self::default()
    }

    /// Adds a node with the given name, ordered input types, output type and payload.
    ///
    /// # Errors
    /// Returns [`TypedError::InvalidName`] or [`TypedError::DuplicateName`].
    pub fn add_node(
        &mut self,
        name: &str,
        inputs: &[TypeTag],
        output: TypeTag,
        payload: T,
    ) -> Result<NodeId, TypedError> {
        if !is_valid_name(name) {
            return Err(TypedError::InvalidName(name.to_owned()));
        }
        if self.node_by_name(name).is_some() {
            return Err(TypedError::DuplicateName(name.to_owned()));
        }
        let id = self.dag.add_node(payload);
        self.ports.push(Ports { name: name.to_owned(), inputs: inputs.to_vec(), output });
        Ok(id)
    }

    /// Adds an edge `from -> to`, filling the next input slot of `to`.
    ///
    /// Re-adding an existing edge succeeds and changes nothing.
    ///
    /// # Errors
    /// Returns [`TypedError::NoInputSlot`] when `to` has no declared slot left,
    /// [`TypedError::TypeMismatch`] when the types differ, and [`TypedError::Dag`]
    /// for unknown nodes or cycles.
    pub fn add_edge(&mut self, from: NodeId, to: NodeId) -> Result<(), TypedError> {
        let out = self.port(from)?.output;
        let target = self.port(to)?;
        if self.dag.predecessors(to).contains(&from) {
            return Ok(());
        }
        let slot = self.dag.predecessors(to).len();
        let Some(&expected) = target.inputs.get(slot) else {
            return Err(TypedError::NoInputSlot { node: to });
        };
        if expected != out {
            return Err(TypedError::TypeMismatch { from, to, slot, expected, found: out });
        }
        self.dag.add_edge(from, to)?;
        Ok(())
    }

    fn port(&self, id: NodeId) -> Result<&Ports, DagError> {
        self.ports.get(id.index()).ok_or(DagError::UnknownNode(id))
    }

    /// Checks that every declared input slot is connected.
    ///
    /// # Errors
    /// Returns [`TypedError::MissingInput`] for the first unconnected slot, in node order.
    pub fn validate(&self) -> Result<(), TypedError> {
        for (i, ports) in self.ports.iter().enumerate() {
            let id = NodeId::from_index(i);
            let connected = self.dag.predecessors(id).len();
            if let Some(&expected) = ports.inputs.get(connected..).and_then(<[TypeTag]>::first) {
                return Err(TypedError::MissingInput { node: id, slot: connected, expected });
            }
        }
        Ok(())
    }

    /// The underlying graph, for running it with the [`Dag`] executors.
    #[must_use]
    pub const fn dag(&self) -> &Dag<T> {
        &self.dag
    }

    /// The payload of a node.
    #[must_use]
    pub fn payload(&self, id: NodeId) -> Option<&T> {
        self.dag.payload(id)
    }

    /// The name a node was added with.
    #[must_use]
    pub fn name(&self, id: NodeId) -> Option<&str> {
        self.ports.get(id.index()).map(|p| p.name.as_str())
    }

    /// The declared input types of a node, in slot order.
    #[must_use]
    pub fn inputs(&self, id: NodeId) -> Option<&[TypeTag]> {
        self.ports.get(id.index()).map(|p| p.inputs.as_slice())
    }

    /// The declared output type of a node.
    #[must_use]
    pub fn output(&self, id: NodeId) -> Option<TypeTag> {
        self.ports.get(id.index()).map(|p| p.output)
    }

    /// Finds a node by name.
    #[must_use]
    pub fn node_by_name(&self, name: &str) -> Option<NodeId> {
        self.ports.iter().position(|p| p.name == name).map(NodeId::from_index)
    }

    /// The ids of all nodes, in insertion order.
    pub fn node_ids(&self) -> impl Iterator<Item = NodeId> + use<T> {
        (0..self.ports.len()).map(NodeId::from_index)
    }

    /// Number of nodes.
    #[must_use]
    pub fn len(&self) -> usize {
        self.ports.len()
    }

    /// Whether the graph has no nodes.
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.ports.is_empty()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const READS: TypeTag = TypeTag::new("seq.reads");
    const COUNTS: TypeTag = TypeTag::new("stats.counts");
    const TABLE: TypeTag = TypeTag::new("stats.table");

    fn sample() -> (TypedDag<()>, NodeId, NodeId, NodeId) {
        let mut g = TypedDag::new();
        let src = g.add_node("load", &[], READS, ()).unwrap();
        let count = g.add_node("count", &[READS], COUNTS, ()).unwrap();
        let report = g.add_node("report", &[COUNTS], TABLE, ()).unwrap();
        g.add_edge(src, count).unwrap();
        g.add_edge(count, report).unwrap();
        (g, src, count, report)
    }

    #[test]
    fn matching_types_connect() {
        let (g, src, count, report) = sample();
        assert_eq!(g.dag().predecessors(count), &[src]);
        assert_eq!(g.dag().predecessors(report).len(), 1);
        assert!(g.validate().is_ok());
        assert_eq!(g.output(src), Some(READS));
        assert_eq!(g.node_by_name("count"), Some(count));
    }

    #[test]
    fn mismatched_type_is_rejected_at_construction() {
        let mut g = TypedDag::new();
        let src = g.add_node("load", &[], READS, ()).unwrap();
        let report = g.add_node("report", &[COUNTS], TABLE, ()).unwrap();
        let err = g.add_edge(src, report).unwrap_err();
        assert_eq!(
            err,
            TypedError::TypeMismatch {
                from: src,
                to: report,
                slot: 0,
                expected: COUNTS,
                found: READS,
            }
        );
        // The rejected edge is not in the graph.
        assert!(g.dag().predecessors(report).is_empty());
    }

    #[test]
    fn extra_edges_have_no_slot() {
        let (mut g, src, count, _) = sample();
        let extra = g.add_node("extra", &[], READS, ()).unwrap();
        // `count` declares one input and already has one edge.
        assert_eq!(g.add_edge(extra, count), Err(TypedError::NoInputSlot { node: count }));
        assert!(g.add_edge(src, count).is_ok(), "re-adding an edge is a no-op");
    }

    #[test]
    fn slots_fill_in_declaration_order() {
        let mut g = TypedDag::new();
        let a = g.add_node("a", &[], READS, ()).unwrap();
        let b = g.add_node("b", &[], COUNTS, ()).unwrap();
        let join = g.add_node("join", &[READS, COUNTS], TABLE, ()).unwrap();
        // Connecting b first would fill slot 0 with COUNTS and fail.
        assert!(matches!(g.add_edge(b, join), Err(TypedError::TypeMismatch { slot: 0, .. })));
        g.add_edge(a, join).unwrap();
        g.add_edge(b, join).unwrap();
        assert!(g.validate().is_ok());
    }

    #[test]
    fn unconnected_inputs_fail_validation() {
        let mut g = TypedDag::new();
        let _a = g.add_node("a", &[], READS, ()).unwrap();
        let _join = g.add_node("join", &[READS], TABLE, ()).unwrap();
        assert!(matches!(g.validate(), Err(TypedError::MissingInput { slot: 0, .. })));
    }

    #[test]
    fn cycles_are_still_rejected() {
        let mut g = TypedDag::new();
        let a = g.add_node("a", &[TABLE], TABLE, ()).unwrap();
        let b = g.add_node("b", &[TABLE], TABLE, ()).unwrap();
        g.add_edge(a, b).unwrap();
        assert_eq!(g.add_edge(b, a), Err(TypedError::Dag(DagError::Cycle { from: b, to: a })));
    }

    #[test]
    fn names_are_validated_and_unique() {
        let mut g: TypedDag<()> = TypedDag::new();
        assert!(matches!(g.add_node("", &[], READS, ()), Err(TypedError::InvalidName(_))));
        assert!(matches!(g.add_node("../x", &[], READS, ()), Err(TypedError::InvalidName(_))));
        assert!(matches!(g.add_node(".hidden", &[], READS, ()), Err(TypedError::InvalidName(_))));
        g.add_node("ok-name_1.v2", &[], READS, ()).unwrap();
        assert!(matches!(
            g.add_node("ok-name_1.v2", &[], READS, ()),
            Err(TypedError::DuplicateName(_))
        ));
    }
}
