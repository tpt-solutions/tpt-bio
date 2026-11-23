//! Two-label foreground/background segmentation by min-cut (graph cut).
//!
//! The energy of a labelling `L` (1 = foreground, 0 = background) is
//! `sum_p cost(p, L_p) + smoothness * #{neighbouring pairs (p, q) with L_p != L_q}`.
//! Each unary cost is shifted by the per-pixel minimum so that all arcs have
//! non-negative capacity, which leaves the minimiser unchanged. The minimiser is
//! found exactly with a max-flow / min-cut on a sparse graph.
//!
//! The max-flow is Dinic's algorithm on an adjacency-list residual graph. The sparse
//! graph is implemented here, because `tpt-math` is not a workspace dependency.

use std::collections::VecDeque;
use std::fmt;

use crate::image::{FACE_OFFSETS, Image, ImageError, OFFSETS_4};

/// Errors from graph-cut segmentation.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum GraphCutError {
    /// The unary cost images do not share one shape.
    Shape(ImageError),
    /// A unary cost is NaN or infinite, or the smoothness weight is.
    NonFinite,
    /// The smoothness weight is negative, which would make the cut non-submodular.
    NegativeSmoothness,
}

impl fmt::Display for GraphCutError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Shape(e) => write!(f, "{e}"),
            Self::NonFinite => f.write_str("graph-cut costs must be finite"),
            Self::NegativeSmoothness => f.write_str("smoothness weight must be non-negative"),
        }
    }
}

impl std::error::Error for GraphCutError {}

/// Pairwise neighbourhood for the smoothness term.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Neighbourhood {
    /// Four in-plane neighbours (2D). On a volume, only the `x` and `y` neighbours are used.
    Four,
    /// Six face neighbours (3D). On a 2D image this is the same as [`Neighbourhood::Four`].
    Six,
}

/// A sparse flow network stored as an adjacency list of paired arcs.
///
/// Arc `e` and its reverse `e ^ 1` are stored together, so the residual capacity of the
/// reverse arc is updated in O(1) when flow is pushed along `e`.
#[derive(Debug, Clone)]
pub struct FlowNetwork {
    adj: Vec<Vec<usize>>,
    to: Vec<usize>,
    cap: Vec<f64>,
}

impl FlowNetwork {
    /// An empty network with `nodes` nodes and no arcs.
    #[must_use]
    pub fn new(nodes: usize) -> Self {
        Self { adj: vec![Vec::new(); nodes], to: Vec::new(), cap: Vec::new() }
    }

    /// Number of nodes.
    #[must_use]
    pub fn node_count(&self) -> usize {
        self.adj.len()
    }

    /// Adds an arc `u -> v` with capacity `forward` and an arc `v -> u` with capacity
    /// `backward`. Passing the same value for both gives an undirected edge.
    pub fn add_edge(&mut self, u: usize, v: usize, forward: f64, backward: f64) {
        self.adj[u].push(self.to.len());
        self.to.push(v);
        self.cap.push(forward);
        self.adj[v].push(self.to.len());
        self.to.push(u);
        self.cap.push(backward);
    }

    /// Maximum flow from `s` to `t`, by Dinic's algorithm. After the call the residual
    /// capacities describe the minimum cut.
    pub fn max_flow(&mut self, s: usize, t: usize) -> f64 {
        let n = self.node_count();
        if s == t {
            return 0.0;
        }
        let mut total = 0.0;
        let mut level = vec![-1i32; n];
        let mut iter = vec![0usize; n];
        while self.bfs_levels(s, t, &mut level) {
            iter.fill(0);
            total += self.blocking_flow(s, t, &level, &mut iter);
        }
        total
    }

    /// Nodes reachable from `s` through arcs with positive residual capacity. After
    /// [`FlowNetwork::max_flow`], this is the source side of a minimum cut.
    #[must_use]
    pub fn reachable_from(&self, s: usize) -> Vec<bool> {
        let mut seen = vec![false; self.node_count()];
        let mut queue = VecDeque::new();
        seen[s] = true;
        queue.push_back(s);
        while let Some(u) = queue.pop_front() {
            for &e in &self.adj[u] {
                let v = self.to[e];
                if self.cap[e] > 0.0 && !seen[v] {
                    seen[v] = true;
                    queue.push_back(v);
                }
            }
        }
        seen
    }

    /// Builds the BFS level graph. Returns true if `t` is reachable.
    fn bfs_levels(&self, s: usize, t: usize, level: &mut [i32]) -> bool {
        level.fill(-1);
        let mut queue = VecDeque::new();
        level[s] = 0;
        queue.push_back(s);
        while let Some(u) = queue.pop_front() {
            for &e in &self.adj[u] {
                let v = self.to[e];
                if self.cap[e] > 0.0 && level[v] < 0 {
                    level[v] = level[u] + 1;
                    queue.push_back(v);
                }
            }
        }
        level[t] >= 0
    }

    /// Pushes a blocking flow along the level graph, without recursion.
    fn blocking_flow(&mut self, s: usize, t: usize, level: &[i32], iter: &mut [usize]) -> f64 {
        let mut total = 0.0;
        let mut path: Vec<usize> = Vec::new();
        let mut u = s;
        loop {
            if u == t {
                let bottleneck = path.iter().map(|&e| self.cap[e]).fold(f64::INFINITY, f64::min);
                for &e in &path {
                    self.cap[e] -= bottleneck;
                    self.cap[e ^ 1] += bottleneck;
                }
                total += bottleneck;
                path.clear();
                u = s;
                continue;
            }
            // Advance along the next admissible arc out of `u`.
            let mut advanced = false;
            while iter[u] < self.adj[u].len() {
                let e = self.adj[u][iter[u]];
                let v = self.to[e];
                if self.cap[e] > 0.0 && level[v] == level[u] + 1 {
                    path.push(e);
                    u = v;
                    advanced = true;
                    break;
                }
                iter[u] += 1;
            }
            if advanced {
                continue;
            }
            // Dead end: retreat, and skip the arc that led here.
            if u == s {
                return total;
            }
            let e = path.pop().expect("a non-source dead end has an arc behind it");
            u = self.to[e ^ 1];
            iter[u] += 1;
        }
    }
}

/// Segments a 2D or 3D image into foreground (`1`) and background (`0`) by min-cut.
///
/// `fg_cost[p]` is the cost of labelling pixel `p` foreground, and `bg_cost[p]` the cost
/// of labelling it background. `smoothness` is paid once for every neighbouring pair
/// with different labels, with neighbours chosen by `neighbourhood`. The result is a
/// label image with `1` for foreground and `0` for background, which is a global
/// minimiser of the energy.
///
/// # Errors
/// Returns [`GraphCutError::Shape`] if the two cost images differ in shape,
/// [`GraphCutError::NonFinite`] for NaN or infinite values, and
/// [`GraphCutError::NegativeSmoothness`] for a negative `smoothness`.
pub fn min_cut_segment(
    fg_cost: &Image<f64>,
    bg_cost: &Image<f64>,
    smoothness: f64,
    neighbourhood: Neighbourhood,
) -> Result<Image<u32>, GraphCutError> {
    fg_cost.check_same_shape(bg_cost).map_err(GraphCutError::Shape)?;
    if !smoothness.is_finite()
        || fg_cost.data().iter().chain(bg_cost.data()).any(|v| !v.is_finite())
    {
        return Err(GraphCutError::NonFinite);
    }
    if smoothness < 0.0 {
        return Err(GraphCutError::NegativeSmoothness);
    }
    let n = fg_cost.len();
    let (source, sink) = (n, n + 1);
    let mut net = FlowNetwork::new(n + 2);
    for p in 0..n {
        let (cf, cb) = (fg_cost.data()[p], bg_cost.data()[p]);
        let shift = cf.min(cb);
        // Source side means foreground. Cutting `s -> p` puts `p` on the sink side
        // (background), so it costs `cb`. Cutting `p -> t` puts `p` on the source side
        // (foreground), so it costs `cf`.
        net.add_edge(source, p, cb - shift, 0.0);
        net.add_edge(p, sink, cf - shift, 0.0);
    }
    let offsets: &[(isize, isize, isize)] = match neighbourhood {
        Neighbourhood::Four => &OFFSETS_4,
        Neighbourhood::Six => &FACE_OFFSETS,
    };
    for p in 0..n {
        for &d in offsets {
            // Visit each undirected pair once, from its smaller index.
            if let Some(q) = fg_cost.offset(p, d)
                && q > p
            {
                net.add_edge(p, q, smoothness, smoothness);
            }
        }
    }
    net.max_flow(source, sink);
    let fg = net.reachable_from(source);
    let labels: Vec<u32> = (0..n).map(|p| u32::from(fg[p])).collect();
    Image::new_3d(fg_cost.width(), fg_cost.height(), fg_cost.depth(), labels)
        .map_err(GraphCutError::Shape)
}

#[allow(clippy::float_cmp)]
#[cfg(test)]
mod tests {
    use super::*;

    /// Energy of a labelling, with the same convention as [`min_cut_segment`].
    fn energy(
        fg: &Image<f64>,
        bg: &Image<f64>,
        labels: &[u32],
        smoothness: f64,
        nb: Neighbourhood,
    ) -> f64 {
        let offsets: &[(isize, isize, isize)] = match nb {
            Neighbourhood::Four => &OFFSETS_4,
            Neighbourhood::Six => &FACE_OFFSETS,
        };
        let mut e = 0.0;
        for p in 0..labels.len() {
            e += if labels[p] == 1 { fg.data()[p] } else { bg.data()[p] };
            for &d in offsets {
                if let Some(q) = fg.offset(p, d)
                    && q > p
                    && labels[p] != labels[q]
                {
                    e += smoothness;
                }
            }
        }
        e
    }

    #[test]
    fn two_pixel_chain_cut_is_known_by_hand() {
        // 5x1 strip: the first two pixels look foreground, the last three background.
        // The unique optimum cuts the single edge between pixels 1 and 2, with energy 1.
        let fg = Image::new_2d(5, 1, vec![0.0, 0.0, 10.0, 10.0, 10.0]).unwrap();
        let bg = Image::new_2d(5, 1, vec![10.0, 10.0, 0.0, 0.0, 0.0]).unwrap();
        let seg = min_cut_segment(&fg, &bg, 1.0, Neighbourhood::Four).unwrap();
        assert_eq!(seg.data(), &[1, 1, 0, 0, 0]);
    }

    #[test]
    fn strong_smoothness_prefers_a_single_label() {
        // With a huge smoothness weight, changing any label costs more than the unaries gain.
        let fg = Image::new_2d(4, 1, vec![0.0, 1.0, 0.0, 1.0]).unwrap();
        let bg = Image::new_2d(4, 1, vec![1.0, 0.0, 1.0, 0.0]).unwrap();
        let seg = min_cut_segment(&fg, &bg, 100.0, Neighbourhood::Four).unwrap();
        assert!(seg.data().iter().all(|&v| v == seg.data()[0]));
    }

    #[test]
    fn three_dimensional_slab_cut_is_known_by_hand() {
        // 2x2x2 volume: the z = 0 layer is foreground, z = 1 is background.
        // The optimum cuts the four z-edges, with energy 4 at smoothness 1.
        let mut fg = Vec::new();
        let mut bg = Vec::new();
        for z in 0..2 {
            for _ in 0..4 {
                fg.push(if z == 0 { 0.0 } else { 10.0 });
                bg.push(if z == 0 { 10.0 } else { 0.0 });
            }
        }
        let fg = Image::new_3d(2, 2, 2, fg).unwrap();
        let bg = Image::new_3d(2, 2, 2, bg).unwrap();
        let seg = min_cut_segment(&fg, &bg, 1.0, Neighbourhood::Six).unwrap();
        assert_eq!(seg.data(), &[1, 1, 1, 1, 0, 0, 0, 0]);
    }

    #[test]
    fn matches_exhaustive_minimum_on_small_random_grids() {
        // Pseudo-random unary costs; the min-cut labelling must match brute force over 2^n.
        for &(w, h) in &[(2usize, 2usize), (3, 2), (4, 2)] {
            for seed in 0..6u64 {
                let n = w * h;
                let mut state = seed.wrapping_mul(0x9E37_79B9_7F4A_7C15) | 1;
                let mut next = || {
                    state ^= state << 13;
                    state ^= state >> 7;
                    state ^= state << 17;
                    (state % 1000) as f64 / 100.0
                };
                let fg: Vec<f64> = (0..n).map(|_| next()).collect();
                let bg: Vec<f64> = (0..n).map(|_| next()).collect();
                let fg = Image::new_2d(w, h, fg).unwrap();
                let bg = Image::new_2d(w, h, bg).unwrap();
                let lambda = 1.5;
                let seg = min_cut_segment(&fg, &bg, lambda, Neighbourhood::Four).unwrap();
                let got = energy(&fg, &bg, seg.data(), lambda, Neighbourhood::Four);
                let mut best = f64::INFINITY;
                for mask in 0..(1u32 << n) {
                    let labels: Vec<u32> = (0..n).map(|p| (mask >> p) & 1).collect();
                    best = best.min(energy(&fg, &bg, &labels, lambda, Neighbourhood::Four));
                }
                assert!((got - best).abs() < 1e-9, "{w}x{h} seed {seed}: {got} vs {best}");
            }
        }
    }

    #[test]
    fn rejects_bad_inputs() {
        let a = Image::new_2d(2, 1, vec![0.0, 0.0]).unwrap();
        let b = Image::new_2d(3, 1, vec![0.0, 0.0, 0.0]).unwrap();
        assert!(matches!(
            min_cut_segment(&a, &b, 1.0, Neighbourhood::Four),
            Err(GraphCutError::Shape(_))
        ));
        assert_eq!(
            min_cut_segment(&a, &a, -1.0, Neighbourhood::Four),
            Err(GraphCutError::NegativeSmoothness)
        );
        let nan = Image::new_2d(2, 1, vec![f64::NAN, 0.0]).unwrap();
        assert_eq!(
            min_cut_segment(&nan, &a, 1.0, Neighbourhood::Four),
            Err(GraphCutError::NonFinite)
        );
    }
}
