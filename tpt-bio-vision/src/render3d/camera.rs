//! Orbit camera with perspective or orthographic projection.
//!
//! World space is right-handed with `z` up. The camera orbits `target` at `distance`: the eye
//! sits at `target + distance * (cos(pitch) cos(yaw), cos(pitch) sin(yaw), sin(pitch))` and
//! looks at the target. Pitch is clamped just short of plus or minus 90 degrees so the view
//! basis stays defined. Pixel coordinates run from `(0, 0)` at the top-left corner to
//! `(width, height)` at the bottom-right, and the target projects to the viewport centre.

use std::f64::consts::{FRAC_PI_2, PI};

use super::{RenderError, add, cross, dot, length, normalise, scale, sub};

const PITCH_LIMIT: f64 = FRAC_PI_2 - 1e-6;
const WORLD_UP: [f64; 3] = [0.0, 0.0, 1.0];

/// Projection model.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum Projection {
    /// Pinhole projection.
    Perspective {
        /// Vertical field of view in radians, in `(0, pi)`.
        fov_y: f64,
    },
    /// Parallel projection.
    Orthographic {
        /// World units visible across the viewport height.
        height: f64,
    },
}

/// A ray with a unit direction.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Ray {
    /// Ray origin in world space.
    pub origin: [f64; 3],
    /// Unit direction.
    pub dir: [f64; 3],
}

/// An orbit camera and its viewport.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Camera {
    /// World point the camera looks at (and orbits).
    pub target: [f64; 3],
    /// Rotation about the `z` axis, in radians. Zero looks along `-x`.
    pub yaw: f64,
    /// Elevation above the `xy` plane, in radians.
    pub pitch: f64,
    /// Distance from the eye to the target, in world units.
    pub distance: f64,
    /// Projection model.
    pub projection: Projection,
    /// Viewport size `(width, height)` in pixels.
    pub viewport: (usize, usize),
}

#[derive(Debug, Clone, Copy)]
struct Basis {
    eye: [f64; 3],
    forward: [f64; 3],
    right: [f64; 3],
    up: [f64; 3],
}

impl Camera {
    /// Checks that every parameter is in range.
    ///
    /// # Errors
    /// Returns [`RenderError::InvalidParameter`] describing the first bad value.
    pub fn validate(&self) -> Result<(), RenderError> {
        let finite = self.target.iter().all(|c| c.is_finite())
            && self.yaw.is_finite()
            && self.pitch.is_finite();
        if !finite {
            return Err(RenderError::InvalidParameter("camera target and angles must be finite"));
        }
        if !(self.distance.is_finite() && self.distance > 0.0) {
            return Err(RenderError::InvalidParameter("camera distance must be positive"));
        }
        if self.viewport.0 == 0 || self.viewport.1 == 0 {
            return Err(RenderError::InvalidParameter("viewport must be non-empty"));
        }
        let ok = match self.projection {
            Projection::Perspective { fov_y } => fov_y > 0.0 && fov_y < PI,
            Projection::Orthographic { height } => height.is_finite() && height > 0.0,
        };
        if !ok {
            return Err(RenderError::InvalidParameter("projection parameter out of range"));
        }
        Ok(())
    }

    /// World position of the eye.
    #[must_use]
    pub fn eye(&self) -> [f64; 3] {
        self.basis().eye
    }

    /// Projects a world point to `[px, py, depth]`, where `(px, py)` are pixel coordinates
    /// and `depth` is the distance along the viewing direction.
    ///
    /// Returns `None` for perspective points at or behind the eye. Orthographic projection
    /// has no such limit.
    #[must_use]
    pub fn project(&self, p: [f64; 3]) -> Option<[f64; 3]> {
        let b = self.basis();
        let d = sub(p, b.eye);
        let (xc, yc, zc) = (dot(d, b.right), dot(d, b.up), dot(d, b.forward));
        let (w, h) = (self.viewport.0 as f64, self.viewport.1 as f64);
        let aspect = w / h;
        let (nx, ny) = match self.projection {
            Projection::Perspective { fov_y } => {
                if zc <= 0.0 {
                    return None;
                }
                let t = (fov_y / 2.0).tan();
                (xc / (zc * t * aspect), yc / (zc * t))
            }
            Projection::Orthographic { height } => {
                let half_h = height / 2.0;
                (xc / (half_h * aspect), yc / half_h)
            }
        };
        Some([(nx + 1.0) * w / 2.0, (1.0 - ny) * h / 2.0, zc])
    }

    /// The ray through pixel coordinate `(px, py)`, which may be a pixel centre such as
    /// `x + 0.5`.
    ///
    /// For perspective cameras the origin is the eye. For orthographic cameras the origin
    /// lies on the eye's image plane, moved `ortho_back` units away from the viewing
    /// direction so that rays start behind the volume. Pass a value larger than the
    /// volume's depth as seen from the eye; perspective cameras ignore it.
    #[must_use]
    pub fn ray(&self, px: f64, py: f64, ortho_back: f64) -> Ray {
        let b = self.basis();
        let (w, h) = (self.viewport.0 as f64, self.viewport.1 as f64);
        let aspect = w / h;
        let nx = 2.0 * px / w - 1.0;
        let ny = 1.0 - 2.0 * py / h;
        match self.projection {
            Projection::Perspective { fov_y } => {
                let t = (fov_y / 2.0).tan();
                let dir = normalise(add(
                    b.forward,
                    add(scale(b.right, nx * t * aspect), scale(b.up, ny * t)),
                ));
                Ray { origin: b.eye, dir }
            }
            Projection::Orthographic { height } => {
                let half_h = height / 2.0;
                let plane =
                    add(b.eye, add(scale(b.right, nx * half_h * aspect), scale(b.up, ny * half_h)));
                Ray { origin: sub(plane, scale(b.forward, ortho_back)), dir: b.forward }
            }
        }
    }

    /// World units covered by one pixel at `depth` (perspective) or everywhere (orthographic).
    ///
    /// Used to turn a physical length into a screen footprint in pixels.
    #[must_use]
    pub fn pixels_per_unit(&self, depth: f64) -> f64 {
        let h = self.viewport.1 as f64;
        match self.projection {
            Projection::Perspective { fov_y } => h / (2.0 * depth * (fov_y / 2.0).tan()),
            Projection::Orthographic { height } => h / height,
        }
    }

    fn basis(&self) -> Basis {
        let pitch = self.pitch.clamp(-PITCH_LIMIT, PITCH_LIMIT);
        let (sin_y, cos_y) = self.yaw.sin_cos();
        let (sin_p, cos_p) = pitch.sin_cos();
        let offset = [cos_p * cos_y, cos_p * sin_y, sin_p];
        let eye = add(self.target, scale(offset, self.distance));
        let forward = normalise(sub(self.target, eye));
        let right = normalise(cross(forward, WORLD_UP));
        let up = cross(right, forward);
        debug_assert!(length(up) > 0.5);
        Basis { eye, forward, right, up }
    }
}

#[allow(clippy::float_cmp)]
#[cfg(test)]
mod tests {
    use super::*;

    fn persp(yaw: f64, pitch: f64) -> Camera {
        Camera {
            target: [5.0, -3.0, 2.0],
            yaw,
            pitch,
            distance: 40.0,
            projection: Projection::Perspective { fov_y: 0.8 },
            viewport: (101, 61),
        }
    }

    #[test]
    fn render3d_camera_target_projects_to_centre() {
        for (yaw, pitch) in [(0.0, 0.0), (0.7, 0.3), (-2.1, -0.9)] {
            let cam = persp(yaw, pitch);
            let [px, py, depth] = cam.project(cam.target).unwrap();
            assert!((px - 50.5).abs() < 1e-9, "px {px}");
            assert!((py - 30.5).abs() < 1e-9, "py {py}");
            assert!((depth - 40.0).abs() < 1e-9);
        }
        let ortho =
            Camera { projection: Projection::Orthographic { height: 20.0 }, ..persp(0.4, 0.2) };
        let [px, py, _] = ortho.project(ortho.target).unwrap();
        assert!((px - 50.5).abs() < 1e-9 && (py - 30.5).abs() < 1e-9);
    }

    #[test]
    fn render3d_camera_centre_ray_points_at_target() {
        let cam = persp(0.3, 0.2);
        let ray = cam.ray(50.5, 30.5, 0.0);
        let to_target = normalise(sub(cam.target, ray.origin));
        assert!(dot(ray.dir, to_target) > 1.0 - 1e-12);
    }

    #[test]
    fn render3d_camera_rejects_bad_parameters() {
        let mut cam = persp(0.0, 0.0);
        assert!(cam.validate().is_ok());
        cam.distance = 0.0;
        assert!(cam.validate().is_err());
        cam.distance = 1.0;
        cam.viewport = (0, 4);
        assert!(cam.validate().is_err());
    }

    #[test]
    fn render3d_camera_ortho_scale_is_exact() {
        let cam = Camera {
            target: [0.0; 3],
            yaw: 0.0,
            pitch: 0.0,
            distance: 100.0,
            projection: Projection::Orthographic { height: 10.0 },
            viewport: (40, 20),
        };
        assert!((cam.pixels_per_unit(100.0) - 2.0).abs() < 1e-12);
        // One world unit along the screen's up axis (+z at yaw 0) moves two pixels.
        let [_, py, _] = cam.project([0.0, 0.0, 1.0]).unwrap();
        assert!((py - 8.0).abs() < 1e-9, "py {py}");
    }
}
