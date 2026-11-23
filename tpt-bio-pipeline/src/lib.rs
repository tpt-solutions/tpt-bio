//! DAG-based workflow model and executor for tpt-bio.
//!
//! A [`Dag`] holds payloads at nodes and directed edges between them. Edges are
//! checked for cycles when they are added, so every `Dag` value is acyclic.
//! [`Dag::execute`] runs the nodes in dependency order, passing each node the
//! results of its predecessors. [`Dag::execute_parallel`] runs every level of
//! independent nodes on the rayon thread pool.
//!
//! The [`Dag`] executors run nodes in memory. [`Progress`] receives a callback when
//! each node starts and finishes. The [`typed`] module adds port types, and
//! [`checkpoint`] persists node outputs so a rerun can resume. [`stream`] provides
//! a bounded, backpressured map over an iterator. The [`definition`] module parses
//! TOML pipeline files, and [`execute`] runs the built-in nodes: statistics, sequence I/O and 2D vision steps.

#![forbid(unsafe_code)]
#![warn(missing_docs)]

pub mod checkpoint;
pub mod definition;
pub mod execute;
pub mod stream;
pub mod typed;

mod seq_nodes;
mod vision_nodes;

use std::cmp::Reverse;
use std::collections::BinaryHeap;
use std::fmt;
use std::sync::atomic::{AtomicBool, AtomicUsize, Ordering};

use rayon::prelude::*;

pub use checkpoint::{CheckpointError, CheckpointStore, Codec};
pub use definition::{DefinitionError, NodeKind, Pipeline};
pub use execute::{Frame, Records, RunError, Value};
pub use typed::{TypeTag, TypedDag, TypedError};

/// Receives a callback as each node starts and finishes.
///
/// Parallel executors call these methods from rayon worker threads, so the
/// implementation must be `Sync`. The counts are approximate under parallelism:
/// `completed` is read when the event fires, not at a single synchronized instant.
pub trait Progress: Sync {
    /// Called just before a node's function runs. `completed` is the number of
    /// nodes that had finished when this node started, and `total` is the node count.
    fn node_started(&self, node: NodeId, completed: usize, total: usize);

    /// Called after a node's function returns. `completed` includes this node.
    /// `ok` is false when the function returned an error.
    fn node_finished(&self, node: NodeId, completed: usize, total: usize, ok: bool);
}

/// A [`Progress`] implementation that ignores every event.
#[derive(Debug, Clone, Copy, Default)]
pub struct NoProgress;

impl Progress for NoProgress {
    fn node_started(&self, _: NodeId, _: usize, _: usize) {}
    fn node_finished(&self, _: NodeId, _: usize, _: usize, _: bool) {}
}

/// Identifier of a node within one [`Dag`].
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct NodeId(usize);

impl NodeId {
    /// The node's index in insertion order.
    #[must_use]
    pub const fn index(self) -> usize {
        self.0
    }

    pub(crate) const fn from_index(index: usize) -> Self {
        Self(index)
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
        self.execute_observed(cancel, &NoProgress, f)
    }

    /// Like [`Dag::execute`], reporting each node's start and finish to `progress`.
    ///
    /// # Errors
    /// Same as [`Dag::execute`].
    pub fn execute_observed<R, E, F, P>(
        &self,
        cancel: &AtomicBool,
        progress: &P,
        f: F,
    ) -> Result<Vec<R>, ExecError<E>>
    where
        F: Fn(NodeId, &T, &[&R]) -> Result<R, E>,
        P: Progress,
    {
        let total = self.len();
        let mut results: Vec<Option<R>> = (0..total).map(|_| None).collect();
        for (completed, n) in self.topological_order().into_iter().enumerate() {
            if cancel.load(Ordering::Relaxed) {
                return Err(ExecError::Cancelled);
            }
            let inputs: Vec<&R> = self.preds[n.0]
                .iter()
                .map(|p| results[p.0].as_ref().expect("predecessors run first"))
                .collect();
            progress.node_started(n, completed, total);
            let outcome = f(n, &self.payloads[n.0], &inputs);
            progress.node_finished(n, completed + 1, total, outcome.is_ok());
            let value = outcome.map_err(|error| ExecError::Node { node: n, error })?;
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
        self.execute_parallel_observed(cancel, &NoProgress, f)
    }

    /// Like [`Dag::execute_parallel`], reporting each node's start and finish to `progress`.
    ///
    /// # Errors
    /// Same as [`Dag::execute_parallel`].
    pub fn execute_parallel_observed<R, E, F, P>(
        &self,
        cancel: &AtomicBool,
        progress: &P,
        f: F,
    ) -> Result<Vec<R>, ExecError<E>>
    where
        R: Send + Sync,
        E: Send,
        F: Fn(NodeId, &T, &[&R]) -> Result<R, E> + Sync,
        P: Progress,
    {
        let total = self.len();
        let finished = AtomicUsize::new(0);
        let mut results: Vec<Option<R>> = (0..total).map(|_| None).collect();
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
                    progress.node_started(n, finished.load(Ordering::SeqCst), total);
                    let outcome = f(n, &self.payloads[n.0], &inputs);
                    let done = finished.fetch_add(1, Ordering::SeqCst) + 1;
                    progress.node_finished(n, done, total, outcome.is_ok());
                    (n, outcome)
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

#[cfg(test)]
mod progress_tests {
    use super::*;
    use std::convert::Infallible;
    use std::sync::Mutex;

    /// Records `(event, node index, completed, ok)` for every callback.
    #[derive(Default)]
    struct Recorder(Mutex<Vec<(&'static str, usize, usize, bool)>>);

    impl Progress for Recorder {
        fn node_started(&self, node: NodeId, completed: usize, _total: usize) {
            self.0.lock().unwrap().push(("start", node.index(), completed, true));
        }
        fn node_finished(&self, node: NodeId, completed: usize, _total: usize, ok: bool) {
            self.0.lock().unwrap().push(("finish", node.index(), completed, ok));
        }
    }

    fn chain(n: usize) -> Dag<u32> {
        let mut g = Dag::new();
        let ids: Vec<NodeId> = (0..n).map(|i| g.add_node(i as u32)).collect();
        for w in ids.windows(2) {
            g.add_edge(w[0], w[1]).unwrap();
        }
        g
    }

    #[test]
    fn sequential_reports_each_node_in_order() {
        let g = chain(3);
        let rec = Recorder::default();
        let cancel = AtomicBool::new(false);
        let out =
            g.execute_observed(&cancel, &rec, |_, v, _: &[&u32]| Ok::<_, Infallible>(*v)).unwrap();
        assert_eq!(out, vec![0, 1, 2]);
        let events = rec.0.lock().unwrap().clone();
        assert_eq!(
            events,
            vec![
                ("start", 0, 0, true),
                ("finish", 0, 1, true),
                ("start", 1, 1, true),
                ("finish", 1, 2, true),
                ("start", 2, 2, true),
                ("finish", 2, 3, true),
            ]
        );
    }

    #[test]
    fn parallel_reports_every_node_once() {
        let g = chain(4);
        let rec = Recorder::default();
        let cancel = AtomicBool::new(false);
        g.execute_parallel_observed(&cancel, &rec, |_, v, _: &[&u32]| Ok::<_, Infallible>(*v))
            .unwrap();
        let events = rec.0.lock().unwrap().clone();
        assert_eq!(events.len(), 8);
        let mut finished: Vec<usize> =
            events.iter().filter(|e| e.0 == "finish").map(|e| e.1).collect();
        finished.sort_unstable();
        assert_eq!(finished, vec![0, 1, 2, 3]);
        assert_eq!(events.iter().filter(|e| e.0 == "finish").map(|e| e.2).max(), Some(4));
    }

    #[test]
    fn failing_node_reports_not_ok() {
        let g = chain(2);
        let rec = Recorder::default();
        let cancel = AtomicBool::new(false);
        let err = g
            .execute_observed(&cancel, &rec, |_, _, _: &[&u32]| Err::<u32, _>("nope"))
            .unwrap_err();
        assert!(matches!(err, ExecError::Node { error: "nope", .. }));
        let events = rec.0.lock().unwrap().clone();
        assert_eq!(events, vec![("start", 0, 0, true), ("finish", 0, 1, false)]);
    }
}
