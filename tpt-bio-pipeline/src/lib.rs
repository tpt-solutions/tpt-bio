//! DAG-based workflow model and executor for tpt-bio.
//!
//! A [`Dag`] holds payloads at nodes and directed edges between them. Edges are
//! checked for cycles when they are added, so every `Dag` value is acyclic.
//! [`Dag::execute`] runs the nodes in dependency order, passing each node the
//! results of its predecessors. [`Dag::execute_parallel`] runs every level of
//! independent nodes on the rayon thread pool.
//!
//! This is the in-memory executor. Streaming between nodes, checkpointing and
//! pipeline files are not implemented yet; see `todo.md`.

#![forbid(unsafe_code)]
#![warn(missing_docs)]

use std::cmp::Reverse;
use std::collections::BinaryHeap;
use std::fmt;
use std::sync::atomic::{AtomicBool, Ordering};

use rayon::prelude::*;

/// Identifier of a node within one [`Dag`].
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct NodeId(usize);

impl NodeId {
    /// The node's index in insertion order.
    #[must_use]
    pub const fn index(self) -> usize {
        self.0
    }
}

/// Errors from building or ordering a graph.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum DagError {
    /// The node ID does not belong to this graph.
    UnknownNode(NodeId),
    /// Adding the edge would create a cycle.
    Cycle {
        /// Source of the rejected edge.
        from: NodeId,
        /// Target of the rejected edge.
        to: NodeId,
    },
}

impl fmt::Display for DagError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::UnknownNode(n) => write!(f, "node {} is not in this graph", n.0),
            Self::Cycle { from, to } => {
                write!(f, "edge {} -> {} would create a cycle", from.0, to.0)
            }
        }
    }
}

impl std::error::Error for DagError {}

/// Errors from running a graph.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum ExecError<E> {
    /// The graph itself is invalid.
    Dag(DagError),
    /// A node's function returned an error.
    Node {
        /// The failing node.
        node: NodeId,
        /// The error it returned.
        error: E,
    },
    /// The run was cancelled before every node finished.
    Cancelled,
}

impl<E: fmt::Display> fmt::Display for ExecError<E> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Dag(e) => e.fmt(f),
            Self::Node { node, error } => write!(f, "node {} failed: {error}", node.0),
            Self::Cancelled => f.write_str("pipeline cancelled"),
        }
    }
}

impl<E: fmt::Debug + fmt::Display> std::error::Error for ExecError<E> {}

/// A directed acyclic graph with a payload `T` at each node.
#[derive(Debug, Clone, Default)]
pub struct Dag<T> {
    payloads: Vec<T>,
    preds: Vec<Vec<NodeId>>,
    succs: Vec<Vec<NodeId>>,
}

impl<T> Dag<T> {
    /// An empty graph.
    #[must_use]
    pub fn new() -> Self {
        Self { payloads: Vec::new(), preds: Vec::new(), succs: Vec::new() }
    }

    /// Adds a node and returns its ID.
    pub fn add_node(&mut self, payload: T) -> NodeId {
        let id = NodeId(self.payloads.len());
        self.payloads.push(payload);
        self.preds.push(Vec::new());
        self.succs.push(Vec::new());
        id
    }

    /// Number of nodes.
    #[must_use]
    pub fn len(&self) -> usize {
        self.payloads.len()
    }

    /// Whether the graph has no nodes.
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.payloads.is_empty()
    }

    /// The payload of a node.
    #[must_use]
    pub fn payload(&self, id: NodeId) -> Option<&T> {
        self.payloads.get(id.0)
    }

    /// The direct predecessors of a node, in the order their edges were added.
    #[must_use]
    pub fn predecessors(&self, id: NodeId) -> &[NodeId] {
        self.preds.get(id.0).map_or(&[], Vec::as_slice)
    }

    fn check(&self, id: NodeId) -> Result<(), DagError> {
        if id.0 < self.payloads.len() { Ok(()) } else { Err(DagError::UnknownNode(id)) }
    }

    /// Adds an edge `from -> to`, meaning `to` runs after `from`.
    ///
    /// Duplicate edges are ignored.
    ///
    /// # Errors
    /// Returns [`DagError::UnknownNode`] for an ID outside this graph, and
    /// [`DagError::Cycle`] if the edge would close a cycle (including a self-loop).
    pub fn add_edge(&mut self, from: NodeId, to: NodeId) -> Result<(), DagError> {
        self.check(from)?;
        self.check(to)?;
        if from == to || self.reaches(to, from) {
            return Err(DagError::Cycle { from, to });
        }
        if !self.succs[from.0].contains(&to) {
            self.succs[from.0].push(to);
            self.preds[to.0].push(from);
        }
        Ok(())
    }

    /// Whether `start` can reach `target` by following edges.
    fn reaches(&self, start: NodeId, target: NodeId) -> bool {
        let mut seen = vec![false; self.payloads.len()];
        let mut stack = vec![start];
        while let Some(n) = stack.pop() {
            if n == target {
                return true;
            }
            if std::mem::replace(&mut seen[n.0], true) {
                continue;
            }
            stack.extend(self.succs[n.0].iter().copied());
        }
        false
    }

    /// A dependency-respecting order. Among ready nodes, the lowest ID comes first,
    /// so the order is deterministic.
    #[must_use]
    pub fn topological_order(&self) -> Vec<NodeId> {
        let mut indegree: Vec<usize> = self.preds.iter().map(Vec::len).collect();
        let mut ready: BinaryHeap<Reverse<NodeId>> =
            (0..self.len()).filter(|&i| indegree[i] == 0).map(|i| Reverse(NodeId(i))).collect();
        let mut order = Vec::with_capacity(self.len());
        while let Some(Reverse(n)) = ready.pop() {
            order.push(n);
            for &s in &self.succs[n.0] {
                indegree[s.0] -= 1;
                if indegree[s.0] == 0 {
                    ready.push(Reverse(s));
                }
            }
        }
        // Edges are checked on insertion, so every node is ordered.
        debug_assert_eq!(order.len(), self.len());
        order
    }

    /// Groups nodes into levels. Every node's predecessors lie in earlier levels,
    /// so the nodes within one level are independent of each other.
    #[must_use]
    pub fn levels(&self) -> Vec<Vec<NodeId>> {
        let mut level_of = vec![0usize; self.len()];
        for n in self.topological_order() {
            let lvl = self.preds[n.0].iter().map(|p| level_of[p.0] + 1).max().unwrap_or(0);
            level_of[n.0] = lvl;
        }
        let depth = level_of.iter().copied().max().map_or(0, |d| d + 1);
        let mut levels = vec![Vec::new(); depth];
        for n in 0..self.len() {
            levels[level_of[n]].push(NodeId(n));
        }
        levels
    }

    /// Runs every node in dependency order on the current thread.
    ///
    /// `f` receives the node ID, its payload, and the results of its predecessors
    /// in edge-insertion order. `cancel` is checked before each node.
    ///
    /// # Errors
    /// Returns the first node error, or [`ExecError::Cancelled`].
    pub fn execute<R, E, F>(&self, cancel: &AtomicBool, f: F) -> Result<Vec<R>, ExecError<E>>
    where
        F: Fn(NodeId, &T, &[&R]) -> Result<R, E>,
    {
        let mut results: Vec<Option<R>> = (0..self.len()).map(|_| None).collect();
        for n in self.topological_order() {
            if cancel.load(Ordering::Relaxed) {
                return Err(ExecError::Cancelled);
            }
            let inputs: Vec<&R> = self.preds[n.0]
                .iter()
                .map(|p| results[p.0].as_ref().expect("predecessors run first"))
                .collect();
            let value = f(n, &self.payloads[n.0], &inputs)
                .map_err(|error| ExecError::Node { node: n, error })?;
            results[n.0] = Some(value);
        }
        Ok(results.into_iter().flatten().collect())
    }
}

impl<T: Sync> Dag<T> {
    /// Runs each level of independent nodes in parallel on the rayon thread pool.
    ///
    /// Results are the same as [`Dag::execute`] for pure functions. If several
    /// nodes in a level fail, the one with the lowest ID is reported.
    ///
    /// # Errors
    /// Returns the first node error by level, or [`ExecError::Cancelled`].
    pub fn execute_parallel<R, E, F>(
        &self,
        cancel: &AtomicBool,
        f: F,
    ) -> Result<Vec<R>, ExecError<E>>
    where
        R: Send + Sync,
        E: Send,
        F: Fn(NodeId, &T, &[&R]) -> Result<R, E> + Sync,
    {
        let mut results: Vec<Option<R>> = (0..self.len()).map(|_| None).collect();
        for level in self.levels() {
            if cancel.load(Ordering::Relaxed) {
                return Err(ExecError::Cancelled);
            }
            let computed: Vec<(NodeId, Result<R, E>)> = level
                .par_iter()
                .map(|&n| {
                    let inputs: Vec<&R> = self.preds[n.0]
                        .iter()
                        .map(|p| results[p.0].as_ref().expect("earlier levels ran first"))
                        .collect();
                    (n, f(n, &self.payloads[n.0], &inputs))
                })
                .collect();
            let mut failed: Option<(NodeId, E)> = None;
            for (n, outcome) in computed {
                match outcome {
                    Ok(v) => results[n.0] = Some(v),
                    Err(e) => {
                        if failed.as_ref().is_none_or(|(prev, _)| n < *prev) {
                            failed = Some((n, e));
                        }
                    }
                }
            }
            if let Some((node, error)) = failed {
                return Err(ExecError::Node { node, error });
            }
        }
        Ok(results.into_iter().flatten().collect())
    }
}

#[cfg(test)]
#[allow(clippy::many_single_char_names)] // Diamond graphs use node names a, b, c, d.
mod tests {
    use super::*;
    use std::convert::Infallible;

    fn diamond() -> (Dag<u32>, [NodeId; 4]) {
        let mut g = Dag::new();
        let a = g.add_node(1);
        let b = g.add_node(2);
        let c = g.add_node(3);
        let d = g.add_node(4);
        g.add_edge(a, b).unwrap();
        g.add_edge(a, c).unwrap();
        g.add_edge(b, d).unwrap();
        g.add_edge(c, d).unwrap();
        (g, [a, b, c, d])
    }

    #[test]
    fn rejects_cycles_and_self_loops() {
        let (mut g, [a, _, _, d]) = diamond();
        assert_eq!(g.add_edge(d, a), Err(DagError::Cycle { from: d, to: a }));
        assert_eq!(g.add_edge(a, a), Err(DagError::Cycle { from: a, to: a }));
        assert_eq!(g.add_edge(a, NodeId(99)), Err(DagError::UnknownNode(NodeId(99))));
    }

    #[test]
    fn duplicate_edges_are_ignored() {
        let (mut g, [a, b, _, _]) = diamond();
        g.add_edge(a, b).unwrap();
        assert_eq!(g.predecessors(b), &[a]);
    }

    #[test]
    fn topological_order_respects_edges() {
        let (g, [a, b, c, d]) = diamond();
        assert_eq!(g.topological_order(), vec![a, b, c, d]);
        assert_eq!(g.levels(), vec![vec![a], vec![b, c], vec![d]]);
    }

    #[test]
    fn sequential_and_parallel_agree() {
        let (g, _) = diamond();
        let cancel = AtomicBool::new(false);
        let f = |_: NodeId, v: &u32, inputs: &[&u64]| -> Result<u64, Infallible> {
            Ok(u64::from(*v) + inputs.iter().map(|x| **x).sum::<u64>())
        };
        let seq = g.execute(&cancel, f).unwrap();
        let par = g.execute_parallel(&cancel, f).unwrap();
        assert_eq!(seq, par);
        assert_eq!(seq[3], 4 + (2 + 1) + (3 + 1));
    }

    #[test]
    fn node_error_stops_the_run() {
        let (g, [_, b, _, _]) = diamond();
        let cancel = AtomicBool::new(false);
        let err = g
            .execute(&cancel, |n, v, _: &[&u32]| if n == b { Err("boom") } else { Ok(*v) })
            .unwrap_err();
        assert_eq!(err, ExecError::Node { node: b, error: "boom" });
    }

    #[test]
    fn cancellation_is_honored() {
        let (g, _) = diamond();
        let cancel = AtomicBool::new(true);
        let err = g.execute(&cancel, |_, v, _: &[&u32]| Ok::<_, Infallible>(*v)).unwrap_err();
        assert_eq!(err, ExecError::Cancelled);
    }
}
