//! Isosurface extraction by marching tetrahedra on a Freudenthal (Kuhn) split.
//!
//! The full 256-case marching-cubes tables are not reproduced here. Instead each grid cell
//! is split into six tetrahedra that all share the cell's main diagonal, one per axis order.
//! That split is conforming: neighbouring cells agree on every shared face, so the mesh has
//! no cracks. Each tetrahedron is cut by the level set with one or two triangles.
//!
//! A grid point is inside when its value is `>= iso`. Triangles are oriented so their normals
//! point out of the inside region, that is towards lower field values. Crossing vertices are
//! shared per grid edge, so the mesh is watertight wherever the level set lies inside the
//! grid. Every edge of a closed surface is used by exactly two triangles.
//!
//! Grid point `(i, j, k)` sits at world position `((i + 0.5) * sx, (j + 0.5) * sy, (k + 0.5) * sz)`,
//! the same convention as [`super::volume::Volume`].

use std::collections::HashMap;

use super::RenderError;
use crate::image::Image;

/// A triangle mesh in world space.
#[derive(Debug, Clone, PartialEq, Default)]
pub struct Mesh {
    /// Vertex positions.
    pub vertices: Vec<[f32; 3]>,
    /// Triangles as counter-clockwise (seen from the outside) vertex indices.
    pub triangles: Vec<[u32; 3]>,
}

impl Mesh {
    /// Enclosed volume, from the signed tetrahedra sum. Exact for closed meshes with
    /// consistently oriented triangles. Returns the absolute value.
    #[must_use]
    pub fn enclosed_volume(&self) -> f64 {
        let mut total = 0.0f64;
        for tri in &self.triangles {
            let a = to_f64(self.vertices[tri[0] as usize]);
            let b = to_f64(self.vertices[tri[1] as usize]);
            let c = to_f64(self.vertices[tri[2] as usize]);
            total += dot(a, cross(b, c)) / 6.0;
        }
        total.abs()
    }
}

/// Extracts the isosurface `field = iso` as a triangle mesh.
///
/// `spacing` is the physical voxel size. A field with an extent below 2 along any axis has
/// no cells and yields an empty mesh.
///
/// # Errors
/// Returns [`RenderError::InvalidParameter`] when `iso` or a spacing component is not finite,
/// or when a spacing component is not positive.
pub fn isosurface(field: &Image<f32>, spacing: [f64; 3], iso: f32) -> Result<Mesh, RenderError> {
    if !iso.is_finite() {
        return Err(RenderError::InvalidParameter("iso value must be finite"));
    }
    if spacing.iter().any(|s| !s.is_finite() || *s <= 0.0) {
        return Err(RenderError::InvalidParameter("voxel spacing must be finite and positive"));
    }
    let (w, h, d) = field.shape();
    let mut mesher = Mesher { field, spacing, iso, mesh: Mesh::default(), edges: HashMap::new() };
    if w >= 2 && h >= 2 && d >= 2 {
        for z in 0..d - 1 {
            for y in 0..h - 1 {
                for x in 0..w - 1 {
                    mesher.cell(x, y, z);
                }
            }
        }
    }
    Ok(mesher.mesh)
}

/// Builds an indicator field for one label: `1.0` where `labels == id`, else `0.0`.
/// Use with `iso = 0.5` to extract the surface of that label.
#[must_use]
pub fn label_indicator(labels: &Image<u32>, id: u32) -> Image<f32> {
    labels.map(|&l| if l == id { 1.0 } else { 0.0 })
}

/// Corner of a cell with its global index, world position and field value.
#[derive(Clone, Copy)]
struct Corner {
    id: usize,
    pos: [f64; 3],
    val: f32,
}

/// Axis orders: each gives one tetrahedron on the path `0 -> e_a -> e_a + e_b -> 7`.
const AXIS_ORDERS: [[usize; 3]; 6] =
    [[0, 1, 2], [0, 2, 1], [1, 0, 2], [1, 2, 0], [2, 0, 1], [2, 1, 0]];

struct Mesher<'a> {
    field: &'a Image<f32>,
    spacing: [f64; 3],
    iso: f32,
    mesh: Mesh,
    /// Crossing vertex per grid edge, keyed by the ordered pair of global indices.
    edges: HashMap<(usize, usize), u32>,
}

impl Mesher<'_> {
    fn cell(&mut self, x: usize, y: usize, z: usize) {
        let data = self.field.data();
        let corners: [Corner; 8] = std::array::from_fn(|c| {
            let (xi, yi, zi) = (x + (c & 1), y + ((c >> 1) & 1), z + ((c >> 2) & 1));
            let id = self.field.index(xi, yi, zi);
            Corner {
                id,
                pos: [
                    (xi as f64 + 0.5) * self.spacing[0],
                    (yi as f64 + 0.5) * self.spacing[1],
                    (zi as f64 + 0.5) * self.spacing[2],
                ],
                val: data[id],
            }
        });
        for order in AXIS_ORDERS {
            let c1 = 1usize << order[0];
            let c2 = c1 | (1usize << order[1]);
            self.tetrahedron([corners[0], corners[c1], corners[c2], corners[7]]);
        }
    }

    fn tetrahedron(&mut self, t: [Corner; 4]) {
        let ins: Vec<Corner> = t.iter().copied().filter(|c| c.val >= self.iso).collect();
        let outs: Vec<Corner> = t.iter().copied().filter(|c| c.val < self.iso).collect();
        match (ins.len(), outs.len()) {
            (1, 3) => {
                let a = ins[0];
                let tri = [
                    self.crossing(a, outs[0]),
                    self.crossing(a, outs[1]),
                    self.crossing(a, outs[2]),
                ];
                self.push(tri, centroid(&ins), centroid(&outs));
            }
            (3, 1) => {
                let o = outs[0];
                let tri =
                    [self.crossing(ins[0], o), self.crossing(ins[1], o), self.crossing(ins[2], o)];
                self.push(tri, centroid(&ins), centroid(&outs));
            }
            (2, 2) => {
                let (a, b) = (ins[0], ins[1]);
                let (c, d) = (outs[0], outs[1]);
                // The four crossings form a quad a-c, a-d, b-d, b-c in cyclic order.
                let q0 = self.crossing(a, c);
                let q1 = self.crossing(a, d);
                let q2 = self.crossing(b, d);
                let q3 = self.crossing(b, c);
                let (ci, co) = (centroid(&ins), centroid(&outs));
                self.push([q0, q1, q2], ci, co);
                self.push([q0, q2, q3], ci, co);
            }
            _ => {}
        }
    }

    /// The vertex where the level set crosses the grid edge between `a` (inside) and `b`
    /// (outside), created on first use.
    fn crossing(&mut self, a: Corner, b: Corner) -> u32 {
        let key = if a.id < b.id { (a.id, b.id) } else { (b.id, a.id) };
        if let Some(&v) = self.edges.get(&key) {
            return v;
        }
        let (va, vb) = (f64::from(a.val), f64::from(b.val));
        let t = (f64::from(self.iso) - va) / (vb - va);
        let p = [0, 1, 2].map(|i| a.pos[i] + (b.pos[i] - a.pos[i]) * t);
        let v = self.mesh.vertices.len() as u32;
        self.mesh.vertices.push(p.map(|c| c as f32));
        self.edges.insert(key, v);
        v
    }

    /// Adds a triangle, flipping it if needed so its normal points from `inside` to `outside`.
    fn push(&mut self, tri: [u32; 3], inside: [f64; 3], outside: [f64; 3]) {
        let a = to_f64(self.mesh.vertices[tri[0] as usize]);
        let b = to_f64(self.mesh.vertices[tri[1] as usize]);
        let c = to_f64(self.mesh.vertices[tri[2] as usize]);
        let normal = cross(sub(b, a), sub(c, a));
        let dir = sub(outside, inside);
        let tri = if dot(normal, dir) < 0.0 { [tri[0], tri[2], tri[1]] } else { tri };
        self.mesh.triangles.push(tri);
    }
}

fn centroid(points: &[Corner]) -> [f64; 3] {
    let n = points.len() as f64;
    let mut out = [0.0; 3];
    for p in points {
        for (o, c) in out.iter_mut().zip(p.pos) {
            *o += c / n;
        }
    }
    out
}

fn to_f64(p: [f32; 3]) -> [f64; 3] {
    p.map(f64::from)
}

fn sub(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]]
}

#[allow(clippy::float_cmp)]
#[cfg(test)]
mod tests {
    use super::*;
    use std::collections::HashMap;

    /// Field `r - distance` from the grid centre: positive inside a sphere of radius `r`.
    fn sphere_field(n: usize, centre: f64, r: f64) -> Image<f32> {
        let mut data = Vec::with_capacity(n * n * n);
        for z in 0..n {
            for y in 0..n {
                for x in 0..n {
                    let d = ((x as f64 - centre).powi(2)
                        + (y as f64 - centre).powi(2)
                        + (z as f64 - centre).powi(2))
                    .sqrt();
                    data.push((r - d) as f32);
                }
            }
        }
        Image::new_3d(n, n, n, data).unwrap()
    }

    /// Counts triangles per undirected edge.
    fn edge_use(mesh: &Mesh) -> HashMap<(u32, u32), usize> {
        let mut uses = HashMap::new();
        for t in &mesh.triangles {
            for k in 0..3 {
                let (a, b) = (t[k], t[(k + 1) % 3]);
                let key = if a < b { (a, b) } else { (b, a) };
                *uses.entry(key).or_insert(0) += 1;
            }
        }
        uses
    }

    #[test]
    fn render3d_isosurface_sphere_is_closed() {
        let field = sphere_field(32, 15.5, 8.0);
        let mesh = isosurface(&field, [1.0; 3], 0.0).unwrap();
        assert!(!mesh.triangles.is_empty());
        let uses = edge_use(&mesh);
        assert!(uses.values().all(|&n| n == 2), "some edge is not shared by exactly two triangles");
    }

    #[test]
    fn render3d_isosurface_sphere_volume_matches_formula() {
        let r = 10.0f64;
        let field = sphere_field(40, 19.5, r);
        let mesh = isosurface(&field, [1.0; 3], 0.0).unwrap();
        let expected = 4.0 / 3.0 * std::f64::consts::PI * r.powi(3);
        let got = mesh.enclosed_volume();
        let rel = (got - expected).abs() / expected;
        assert!(rel < 0.03, "volume {got}, expected {expected}, relative error {rel}");
    }

    #[test]
    fn render3d_isosurface_respects_spacing() {
        let field = sphere_field(40, 19.5, 10.0);
        let unit = isosurface(&field, [1.0; 3], 0.0).unwrap().enclosed_volume();
        let scaled = isosurface(&field, [2.0, 2.0, 2.0], 0.0).unwrap().enclosed_volume();
        assert!((scaled / unit - 8.0).abs() < 1e-3);
    }

    #[test]
    fn render3d_isosurface_label_indicator_sphere() {
        let n = 30usize;
        let mut labels = vec![0u32; n * n * n];
        for z in 0..n {
            for y in 0..n {
                for x in 0..n {
                    let d2 = (x as f64 - 14.5).powi(2)
                        + (y as f64 - 14.5).powi(2)
                        + (z as f64 - 14.5).powi(2);
                    if d2 <= 36.0 {
                        labels[(z * n + y) * n + x] = 5;
                    }
                }
            }
        }
        let labels = Image::new_3d(n, n, n, labels).unwrap();
        let field = label_indicator(&labels, 5);
        let mesh = isosurface(&field, [1.0; 3], 0.5).unwrap();
        assert!(!mesh.triangles.is_empty());
        assert!(edge_use(&mesh).values().all(|&n| n == 2));
    }

    #[test]
    fn render3d_isosurface_flat_field_is_empty() {
        let field = Image::new_3d(4, 4, 4, vec![1.0f32; 64]).unwrap();
        let mesh = isosurface(&field, [1.0; 3], 0.5).unwrap();
        assert!(mesh.triangles.is_empty());
        let thin = Image::new_3d(4, 4, 1, vec![1.0f32; 16]).unwrap();
        assert!(isosurface(&thin, [1.0; 3], 0.5).unwrap().triangles.is_empty());
    }
}
