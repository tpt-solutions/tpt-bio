//! Transfer functions: windowed, piecewise-linear colour and opacity per channel, and a
//! clipping plane.
//!
//! A sample value is first normalised through its channel's window `(min, max)` to
//! `t` in `[0, 1]` (values outside the window clamp). The colour and opacity knots are then
//! evaluated at `t`. Colour is multiplied by the channel's `gain` and clamped to `[0, 1]`.
//!
//! Several channels combine at one sample point: opacity composes as
//! `1 - prod(1 - a_i)`, and colour is the opacity-weighted mean of the channel colours.

use super::{RenderError, dot};

/// Colour and opacity mapping for one channel.
#[derive(Debug, Clone, PartialEq)]
pub struct ChannelTransfer {
    /// Window `(min, max)`. Values map linearly to `t` in `[0, 1]` across it.
    pub window: (f32, f32),
    /// Multiplier applied to the colour, clamped to `[0, 1]` afterwards.
    pub gain: f32,
    /// Colour knots `(t, rgb)` sorted by `t`, each `rgb` component in `[0, 1]`.
    pub colour: Vec<(f32, [f32; 3])>,
    /// Opacity knots `(t, alpha)` sorted by `t`, alpha in `[0, 1]`.
    pub opacity: Vec<(f32, f32)>,
}

impl ChannelTransfer {
    /// A channel with a constant colour, unit gain, and opacity rising linearly from 0 at
    /// the window minimum to 1 at the maximum.
    #[must_use]
    pub fn new(colour: [f32; 3], window: (f32, f32)) -> Self {
        Self {
            window,
            gain: 1.0,
            colour: vec![(0.0, colour), (1.0, colour)],
            opacity: vec![(0.0, 0.0), (1.0, 1.0)],
        }
    }

    /// Maps a raw sample value to `t` in `[0, 1]` through the window.
    #[must_use]
    pub fn normalise(&self, value: f32) -> f32 {
        let (lo, hi) = self.window;
        ((value - lo) / (hi - lo)).clamp(0.0, 1.0)
    }

    /// The colour (after gain) and opacity for a raw sample value.
    ///
    /// Non-finite values map to transparent black.
    #[must_use]
    pub fn evaluate(&self, value: f32) -> ([f32; 3], f32) {
        if !value.is_finite() {
            return ([0.0; 3], 0.0);
        }
        let t = self.normalise(value);
        let mut colour =
            piecewise(&self.colour, t, |a, b, f| [0, 1, 2].map(|i| a[i] + (b[i] - a[i]) * f));
        for c in &mut colour {
            *c = (*c * self.gain).clamp(0.0, 1.0);
        }
        let alpha = piecewise(&self.opacity, t, |a, b, f| a + (b - a) * f).clamp(0.0, 1.0);
        (colour, alpha)
    }
}

/// A plane `dot(normal, p) + offset = 0`. Points with a value below zero are on the
/// discarded side.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ClipPlane {
    /// Plane normal in world space. Need not be unit length.
    pub normal: [f64; 3],
    /// Signed offset of the plane.
    pub offset: f64,
}

impl ClipPlane {
    /// True when `p` lies on the kept side, where `dot(normal, p) + offset >= 0`.
    #[must_use]
    pub fn keeps(&self, p: [f64; 3]) -> bool {
        dot(self.normal, p) + self.offset >= 0.0
    }
}

/// Per-channel transfer functions plus an optional clipping plane.
#[derive(Debug, Clone, PartialEq)]
pub struct TransferFunction {
    /// One entry per volume channel, in channel order.
    pub channels: Vec<ChannelTransfer>,
    /// Samples on the discarded side of this plane contribute nothing.
    pub clip: Option<ClipPlane>,
}

impl TransferFunction {
    /// Checks windows, gains and knot lists.
    ///
    /// # Errors
    /// Returns [`RenderError::InvalidParameter`] for the first invalid entry.
    pub fn validate(&self) -> Result<(), RenderError> {
        if self.channels.is_empty() {
            return Err(RenderError::InvalidParameter("transfer function has no channels"));
        }
        for ch in &self.channels {
            let (lo, hi) = ch.window;
            if !(lo.is_finite() && hi.is_finite() && hi > lo) {
                return Err(RenderError::InvalidParameter("window must be finite with max > min"));
            }
            if !(ch.gain.is_finite() && ch.gain >= 0.0) {
                return Err(RenderError::InvalidParameter("gain must be finite and non-negative"));
            }
            if ch.colour.is_empty() || ch.opacity.is_empty() {
                return Err(RenderError::InvalidParameter("knot lists must not be empty"));
            }
            if !sorted(ch.colour.iter().map(|k| k.0)) || !sorted(ch.opacity.iter().map(|k| k.0)) {
                return Err(RenderError::InvalidParameter("knots must be sorted by t"));
            }
        }
        if let Some(c) = self.clip
            && !(c.offset.is_finite() && c.normal.iter().all(|v| v.is_finite()))
        {
            return Err(RenderError::InvalidParameter("clip plane must be finite"));
        }
        Ok(())
    }

    /// Combines one value per channel into a single `(rgb, alpha)` sample.
    ///
    /// `values` must have one entry per channel; extra entries are ignored.
    #[must_use]
    pub fn combine(&self, values: &[f32]) -> ([f32; 3], f32) {
        let mut transparency = 1.0f32;
        let mut weighted = [0.0f32; 3];
        let mut alpha_sum = 0.0f32;
        for (ch, &v) in self.channels.iter().zip(values) {
            let (c, a) = ch.evaluate(v);
            for (acc, cv) in weighted.iter_mut().zip(c) {
                *acc += cv * a;
            }
            alpha_sum += a;
            transparency *= 1.0 - a;
        }
        if alpha_sum > 0.0 {
            let colour = weighted.map(|w| w / alpha_sum);
            (colour, 1.0 - transparency)
        } else {
            ([0.0; 3], 0.0)
        }
    }
}

/// Evaluates a sorted knot list at `t` with clamping outside the knot range.
fn piecewise<T: Copy>(knots: &[(f32, T)], t: f32, mix: impl Fn(T, T, f32) -> T) -> T {
    let first = knots[0];
    let last = knots[knots.len() - 1];
    if t <= first.0 {
        return first.1;
    }
    if t >= last.0 {
        return last.1;
    }
    for pair in knots.windows(2) {
        if let [a, b] = pair
            && t <= b.0
        {
            let span = b.0 - a.0;
            if span <= 0.0 {
                return b.1;
            }
            return mix(a.1, b.1, (t - a.0) / span);
        }
    }
    last.1
}

fn sorted(mut ts: impl Iterator<Item = f32>) -> bool {
    let Some(mut prev) = ts.next() else { return true };
    for t in ts {
        if t.is_nan() || t < prev {
            return false;
        }
        prev = t;
    }
    true
}

#[allow(clippy::float_cmp)]
#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn render3d_transfer_maps_knots_exactly() {
        let ch = ChannelTransfer {
            window: (100.0, 300.0),
            gain: 1.0,
            colour: vec![(0.0, [0.0, 0.0, 0.0]), (0.5, [1.0, 0.0, 0.0]), (1.0, [1.0, 1.0, 1.0])],
            opacity: vec![(0.0, 0.0), (0.5, 0.75), (1.0, 0.25)],
        };
        // Window position 100 -> t 0, 200 -> t 0.5, 300 -> t 1.
        assert_eq!(ch.evaluate(100.0), ([0.0, 0.0, 0.0], 0.0));
        assert_eq!(ch.evaluate(200.0), ([1.0, 0.0, 0.0], 0.75));
        assert_eq!(ch.evaluate(300.0), ([1.0, 1.0, 1.0], 0.25));
        // Values outside the window clamp to the end knots.
        assert_eq!(ch.evaluate(0.0).1, 0.0);
        assert_eq!(ch.evaluate(1e6).1, 0.25);
        // Midway between the first two knots: opacity 0.375.
        assert!((ch.evaluate(150.0).1 - 0.375).abs() < 1e-6);
        assert_eq!(ch.evaluate(f32::NAN), ([0.0; 3], 0.0));
    }

    #[test]
    fn render3d_transfer_gain_scales_colour_with_clamp() {
        let mut ch = ChannelTransfer::new([0.4, 0.5, 0.6], (0.0, 1.0));
        ch.gain = 2.0;
        let (c, a) = ch.evaluate(1.0);
        assert_eq!(c, [0.8, 1.0, 1.0]);
        assert_eq!(a, 1.0);
    }

    #[test]
    fn render3d_transfer_combines_channels() {
        let red = ChannelTransfer::new([1.0, 0.0, 0.0], (0.0, 1.0));
        let blue = ChannelTransfer::new([0.0, 0.0, 1.0], (0.0, 1.0));
        let tf = TransferFunction { channels: vec![red, blue], clip: None };
        assert!(tf.validate().is_ok());
        // Equal opacity 0.5 in each channel: alpha 0.75, colour the even mix.
        let (c, a) = tf.combine(&[0.5, 0.5]);
        assert!((a - 0.75).abs() < 1e-6);
        assert!((c[0] - 0.5).abs() < 1e-6 && (c[2] - 0.5).abs() < 1e-6);
    }

    #[test]
    fn render3d_transfer_clip_plane_discards_negative_side() {
        let plane = ClipPlane { normal: [1.0, 0.0, 0.0], offset: -2.0 };
        assert!(!plane.keeps([1.0, 0.0, 0.0]));
        assert!(plane.keeps([2.0, 5.0, 5.0]));
        assert!(plane.keeps([3.0, 0.0, 0.0]));
    }

    #[test]
    fn render3d_transfer_validate_rejects_bad_knots() {
        let mut ch = ChannelTransfer::new([1.0; 3], (0.0, 1.0));
        ch.opacity = vec![(0.8, 0.0), (0.2, 1.0)];
        let tf = TransferFunction { channels: vec![ch], clip: None };
        assert!(tf.validate().is_err());
        let flat = ChannelTransfer::new([1.0; 3], (1.0, 1.0));
        assert!(TransferFunction { channels: vec![flat], clip: None }.validate().is_err());
    }
}
