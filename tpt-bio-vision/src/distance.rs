//! Exact Euclidean distance transform (Felzenszwalb and Huttenlocher, 2012).
//!
//! For each foreground pixel, the output is the distance to the nearest background
//! pixel (`mask == false`). Background pixels get `0`. The transform is separable:
//! a 1D squared-distance lower envelope is applied along each axis in turn.
//! If the mask has no background at all, every distance is `f64::INFINITY`.

use crate::image::Image;

/// Stand-in for "no background seen yet" in the squared-distance pass. Large enough
/// that it never wins a lower-envelope comparison, small enough that `f64` keeps
/// integer precision when it is combined with squared offsets.
const INF: f64 = 1e15;

/// Exact squared Euclidean distance from each foreground voxel to the nearest background
/// voxel. Background voxels get `0`. Distances are in voxel units, measured between
/// voxel centres.
#[must_use]
pub fn squared_distance(mask: &Image<bool>) -> Image<f64> {
    let (w, h, d) = mask.shape();
    let mut f: Vec<f64> = mask.data().iter().map(|&m| if m { INF } else { 0.0 }).collect();
    let mut scratch = Scratch::default();
    for axis in 0..3 {
        for line in axis_lines(w, h, d, axis) {
            scratch.transform_line(&mut f, &line);
        }
    }

    let data = f.into_iter().map(|v| if v >= INF / 2.0 { f64::INFINITY } else { v }).collect();
    Image::new_3d(w, h, d, data).expect("distance buffer matches mask shape")
}

/// Exact Euclidean distance from each foreground voxel to the nearest background voxel.
///
/// Works for 2D (`depth == 1`) and 3D images.
#[must_use]
pub fn distance_transform(mask: &Image<bool>) -> Image<f64> {
    squared_distance(mask).map(|&v| v.sqrt())
}

/// Index lists for every 1D line along `axis` (0 = x, 1 = y, 2 = z).
fn axis_lines(w: usize, h: usize, d: usize, axis: usize) -> Vec<Vec<usize>> {
    let (len, stride) = match axis {
        0 => (w, 1),
        1 => (h, w),
        _ => (d, w * h),
    };
    let mut out = Vec::new();
    for z in 0..d {
        for y in 0..h {
            for x in 0..w {
                let coord = [x, y, z][axis];
                if coord == 0 {
                    let base = (z * h + y) * w + x;
                    out.push((0..len).map(|k| base + k * stride).collect());
                }
            }
        }
    }
    out
}

#[derive(Default)]
struct Scratch {
    v: Vec<usize>,
    z: Vec<f64>,
    line: Vec<f64>,
    out: Vec<f64>,
}

impl Scratch {
    fn transform_line(&mut self, f: &mut [f64], idx: &[usize]) {
        let n = idx.len();
        if n == 0 {
            return;
        }
        self.line.clear();
        self.line.extend(idx.iter().map(|&i| f[i]));
        self.v.clear();
        self.v.resize(n, 0);
        self.z.clear();
        self.z.resize(n + 1, 0.0);
        self.out.clear();
        self.out.resize(n, 0.0);
        lower_envelope(&self.line, &mut self.v, &mut self.z, &mut self.out);
        for (&i, &val) in idx.iter().zip(&self.out) {
            f[i] = val;
        }
    }
}

/// 1D squared distance transform of `g` (Felzenszwalb and Huttenlocher), written to `out`.
fn lower_envelope(g: &[f64], v: &mut [usize], z: &mut [f64], out: &mut [f64]) {
    let n = g.len();
    let mut k = 0usize;
    v[0] = 0;
    z[0] = f64::NEG_INFINITY;
    z[1] = f64::INFINITY;
    for q in 1..n {
        let mut s = intersection(g, v[k], q);
        while s <= z[k] {
            k -= 1;
            s = intersection(g, v[k], q);
        }
        k += 1;
        v[k] = q;
        z[k] = s;
        z[k + 1] = f64::INFINITY;
    }
    k = 0;
    for (q, slot) in out.iter_mut().enumerate() {
        while z[k + 1] < q as f64 {
            k += 1;
        }
        let dq = q as f64 - v[k] as f64;
        *slot = dq * dq + g[v[k]];
    }
}

/// Abscissa where the parabolas rooted at `p` and `q` intersect.
fn intersection(g: &[f64], p: usize, q: usize) -> f64 {
    let (pf, qf) = (p as f64, q as f64);
    ((g[q] + qf * qf) - (g[p] + pf * pf)) / (2.0 * qf - 2.0 * pf)
}

#[allow(clippy::float_cmp)]
#[cfg(test)]
mod tests {
    use super::*;

    /// Brute-force reference: squared distance to the nearest background voxel.
    fn brute(mask: &Image<bool>) -> Vec<f64> {
        let bg: Vec<(usize, usize, usize)> =
            (0..mask.len()).filter(|&i| !mask.data()[i]).map(|i| mask.coords(i)).collect();
        (0..mask.len())
            .map(|i| {
                if !mask.data()[i] {
                    return 0.0;
                }
                let (x, y, z) = mask.coords(i);
                bg.iter()
                    .map(|&(bx, by, bz)| {
                        let dx = x as f64 - bx as f64;
                        let dy = y as f64 - by as f64;
                        let dz = z as f64 - bz as f64;
                        dx * dx + dy * dy + dz * dz
                    })
                    .fold(f64::INFINITY, f64::min)
            })
            .collect()
    }

    #[test]
    fn single_background_pixel_gives_radial_distances() {
        // 5x5 foreground with the centre pixel as the only background.
        let mut data = vec![true; 25];
        data[12] = false;
        let m = Image::new_2d(5, 5, data).unwrap();
        let sq = squared_distance(&m);
        assert_eq!(*sq.get(2, 2, 0).unwrap(), 0.0);
        assert_eq!(*sq.get(3, 2, 0).unwrap(), 1.0);
        assert_eq!(*sq.get(3, 3, 0).unwrap(), 2.0);
        assert_eq!(*sq.get(0, 0, 0).unwrap(), 8.0);
        let d = distance_transform(&m);
        assert!((d.get(4, 4, 0).unwrap() - 8f64.sqrt()).abs() < 1e-12);
    }

    #[test]
    fn full_foreground_has_infinite_distance() {
        let m = Image::new_2d(2, 2, vec![true; 4]).unwrap();
        assert!(distance_transform(&m).data().iter().all(|v| v.is_infinite()));
    }

    #[test]
    fn matches_brute_force_on_irregular_2d_and_3d() {
        for &(w, h, d) in &[(1usize, 1usize, 1usize), (7, 5, 1), (9, 4, 3)] {
            let data: Vec<bool> = (0..w * h * d).map(|i| (i * 37 + 11) % 5 != 0).collect();
            let m = Image::new_3d(w, h, d, data).unwrap();
            let got = squared_distance(&m);
            let want = brute(&m);
            for (g, r) in got.data().iter().zip(&want) {
                if r.is_infinite() {
                    assert!(g.is_infinite(), "{w}x{h}x{d}: {g} vs {r}");
                } else {
                    assert!((g - r).abs() < 1e-9, "{w}x{h}x{d}: {g} vs {r}");
                }
            }
        }
    }

    #[test]
    fn three_d_distance_along_axis() {
        let mut data = vec![true; 5];
        data[0] = false;
        let m = Image::new_3d(1, 1, 5, data).unwrap();
        let d = distance_transform(&m);
        assert_eq!(d.data(), &[0.0, 1.0, 2.0, 3.0, 4.0]);
    }
}
