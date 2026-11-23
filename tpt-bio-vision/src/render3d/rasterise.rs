//! Flat-shaded software rasteriser for [`Mesh`] with a z-buffer.
//!
//! Each triangle gets one colour from its world-space normal and a fixed light direction
//! (Lambert, two-sided, with an ambient floor). Triangles are filled by testing pixel centres
//! with edge functions. Depth is interpolated linearly in screen space, which is exact for
//! orthographic cameras and approximate for perspective ones. There is no near-plane
//! clipping: a triangle with any vertex at or behind a perspective eye is skipped.

use super::camera::Camera;
use super::isosurface::Mesh;
use super::{RenderError, Rgba8, cross, dot, normalise, sub};
use crate::image::Image;

/// Ambient light fraction, so faces turned away from the light stay visible.
const AMBIENT: f64 = 0.2;

/// Renders `mesh` seen through `camera` into an RGBA image of the viewport.
///
/// `colour` is the surface albedo and `background` fills pixels no triangle covers.
///
/// # Errors
/// Returns [`RenderError::InvalidParameter`] for an invalid camera or a triangle index
/// that is out of range.
pub fn rasterise_mesh(
    mesh: &Mesh,
    camera: &Camera,
    colour: [u8; 3],
    background: [u8; 3],
) -> Result<Image<Rgba8>, RenderError> {
    camera.validate()?;
    let n = mesh.vertices.len();
    if mesh.triangles.iter().any(|t| t.iter().any(|&i| i as usize >= n)) {
        return Err(RenderError::InvalidParameter("triangle index out of range"));
    }
    let (w, h) = camera.viewport;
    let mut color_buf = vec![[background[0], background[1], background[2], 255u8]; w * h];
    let mut zbuf = vec![f64::INFINITY; w * h];
    let light = normalise([0.4, -0.5, 0.8]);
    let projected: Vec<Option<[f64; 3]>> =
        mesh.vertices.iter().map(|v| camera.project(v.map(f64::from))).collect();

    for tri in &mesh.triangles {
        let idx = tri.map(|i| i as usize);
        let [Some(pa), Some(pb), Some(pc)] = idx.map(|i| projected[i]) else {
            continue;
        };
        let world = idx.map(|i| mesh.vertices[i].map(f64::from));
        let normal = normalise(cross(sub(world[1], world[0]), sub(world[2], world[0])));
        if length_is_zero(normal) {
            continue;
        }
        let shade = AMBIENT + (1.0 - AMBIENT) * dot(normal, light).abs();
        let pixel = [
            (f64::from(colour[0]) * shade).round().min(255.0) as u8,
            (f64::from(colour[1]) * shade).round().min(255.0) as u8,
            (f64::from(colour[2]) * shade).round().min(255.0) as u8,
            255,
        ];
        fill_triangle([pa, pb, pc], w, h, |x, y, z| {
            let i = y * w + x;
            if z < zbuf[i] {
                zbuf[i] = z;
                color_buf[i] = pixel;
            }
        });
    }

    Ok(Image::new_2d(w, h, color_buf)?)
}

fn length_is_zero(v: [f64; 3]) -> bool {
    dot(v, v) == 0.0
}

/// Calls `emit(x, y, z)` for every pixel whose centre lies inside the screen triangle.
/// `z` is the depth interpolated at the pixel centre.
fn fill_triangle(p: [[f64; 3]; 3], w: usize, h: usize, mut emit: impl FnMut(usize, usize, f64)) {
    let area = edge(p[0], p[1], p[2]);
    if area.abs() < 1e-12 {
        return;
    }
    let min_x = p.iter().map(|q| q[0]).fold(f64::INFINITY, f64::min);
    let max_x = p.iter().map(|q| q[0]).fold(f64::NEG_INFINITY, f64::max);
    let min_y = p.iter().map(|q| q[1]).fold(f64::INFINITY, f64::min);
    let max_y = p.iter().map(|q| q[1]).fold(f64::NEG_INFINITY, f64::max);
    if max_x < 0.0 || max_y < 0.0 || min_x > w as f64 || min_y > h as f64 {
        return;
    }
    let x0 = min_x.floor().max(0.0) as usize;
    let x1 = (max_x.ceil().max(0.0) as usize).min(w - 1);
    let y0 = min_y.floor().max(0.0) as usize;
    let y1 = (max_y.ceil().max(0.0) as usize).min(h - 1);
    for y in y0..=y1 {
        for x in x0..=x1 {
            let c = [x as f64 + 0.5, y as f64 + 0.5, 0.0];
            let w0 = edge(p[1], p[2], c) / area;
            let w1 = edge(p[2], p[0], c) / area;
            let w2 = edge(p[0], p[1], c) / area;
            if w0 >= 0.0 && w1 >= 0.0 && w2 >= 0.0 {
                let z = w0 * p[0][2] + w1 * p[1][2] + w2 * p[2][2];
                emit(x, y, z);
            }
        }
    }
}

/// Twice the signed area of the screen triangle `(a, b, p)`, using only `x` and `y`.
fn edge(a: [f64; 3], b: [f64; 3], p: [f64; 3]) -> f64 {
    (b[0] - a[0]) * (p[1] - a[1]) - (b[1] - a[1]) * (p[0] - a[0])
}

#[allow(clippy::float_cmp)]
#[cfg(test)]
mod tests {
    use super::*;
    use crate::render3d::camera::Projection;
    use crate::render3d::isosurface::isosurface;

    /// Field `r - distance` on an `n^3` grid, centred at `centre` (grid index units).
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

    #[test]
    fn render3d_rasterise_sphere_covers_expected_disk() {
        // Grid index centre 19.5 is world 20.0 (grid point i sits at world i + 0.5).
        let field = sphere_field(40, 19.5, 10.0);
        let mesh = isosurface(&field, [1.0; 3], 0.0).unwrap();
        // Orthographic, 64 px for 32 units: 2 px per unit, so the disk radius is 20 px.
        let cam = Camera {
            target: [20.0, 20.0, 20.0],
            yaw: 0.0,
            pitch: 0.0,
            distance: 100.0,
            projection: Projection::Orthographic { height: 32.0 },
            viewport: (64, 64),
        };
        let img = rasterise_mesh(&mesh, &cam, [200, 100, 50], [0, 0, 0]).unwrap();
        let covered = img.data().iter().filter(|p| p[3] == 255 && **p != [0, 0, 0, 255]).count();
        let radius_px = (covered as f64 / std::f64::consts::PI).sqrt();
        assert!((radius_px - 20.0).abs() < 1.0, "estimated radius {radius_px} px");
        // The centre pixel is covered and the corner is background.
        assert_ne!(*img.get(32, 32, 0).unwrap(), [0, 0, 0, 255]);
        assert_eq!(*img.get(0, 0, 0).unwrap(), [0, 0, 0, 255]);
    }

    #[test]
    fn render3d_rasterise_rejects_bad_index() {
        let mesh = Mesh { vertices: vec![[0.0; 3]], triangles: vec![[0, 1, 2]] };
        let cam = Camera {
            target: [0.0; 3],
            yaw: 0.0,
            pitch: 0.0,
            distance: 10.0,
            projection: Projection::Orthographic { height: 4.0 },
            viewport: (8, 8),
        };
        assert!(rasterise_mesh(&mesh, &cam, [255; 3], [0; 3]).is_err());
    }
}
