// Front-to-back ray marching for `render3d::gpu`. One invocation per pixel.
//
// Mirrors `render3d::raycast::march` on the CPU: the ray arrives already clipped to the
// volume box (t0 clamped to >= 0), marching starts at t0 + step / 2, samples are trilinear,
// per-channel transfer functions are read from a lookup table and combined as in
// `TransferFunction::combine`, and compositing stops at the opacity cutoff.

struct Params {
    // width, height, depth, channel count
    dims: vec4<u32>,
    // spacing x, y, z, march step
    spacing_step: vec4<f32>,
    // lookup table length, clip plane enabled (0 or 1), viewport width, viewport height
    counts: vec4<u32>,
    // clip plane normal x, y, z and offset: a point is kept when dot(n, p) + offset >= 0
    clip: vec4<f32>,
    // opacity cutoff in x
    cutoff: vec4<f32>,
};

@group(0) @binding(0) var<uniform> P: Params;
// Channel-major voxel data: channel c, voxel index (z * h + y) * w + x.
@group(0) @binding(1) var<storage, read> vol: array<f32>;
// Two entries per pixel: (origin, t0) then (direction, t1).
@group(0) @binding(2) var<storage, read> rays: array<vec4<f32>>;
// Per channel: lookup table of (rgb after gain, alpha) over the window, LUT length entries each.
@group(0) @binding(3) var<storage, read> lut: array<vec4<f32>>;
// Per channel window (min, max).
@group(0) @binding(4) var<storage, read> windows: array<vec2<f32>>;
// Premultiplied colour and opacity per pixel.
@group(0) @binding(5) var<storage, read_write> out_px: array<vec4<f32>>;

fn sample_channel(c: u32, p: vec3<f32>) -> f32 {
    let w = P.dims.x;
    let h = P.dims.y;
    let d = P.dims.z;
    let n = w * h * d;
    let last = vec3<f32>(f32(w) - 1.0, f32(h) - 1.0, f32(d) - 1.0);
    let u = clamp(p / P.spacing_step.xyz - vec3<f32>(0.5), vec3<f32>(0.0), last);
    let lo = floor(u);
    let fr = u - lo;
    let loi = vec3<u32>(lo);
    let hii = min(loi + vec3<u32>(1u), vec3<u32>(w - 1u, h - 1u, d - 1u));
    var out = 0.0;
    for (var k = 0u; k < 8u; k = k + 1u) {
        let hx = (k & 1u) != 0u;
        let hy = (k & 2u) != 0u;
        let hz = (k & 4u) != 0u;
        let x = select(loi.x, hii.x, hx);
        let y = select(loi.y, hii.y, hy);
        let z = select(loi.z, hii.z, hz);
        let wt = select(1.0 - fr.x, fr.x, hx)
            * select(1.0 - fr.y, fr.y, hy)
            * select(1.0 - fr.z, fr.z, hz);
        out = out + wt * vol[c * n + (z * h + y) * w + x];
    }
    return out;
}

// Colour (after gain) and opacity for a raw value, from the lookup table. Linear in the
// table index, so it approximates `ChannelTransfer::evaluate` to within the table spacing.
fn lut_entry(c: u32, v: f32) -> vec4<f32> {
    let win = windows[c];
    let t = clamp((v - win.x) / (win.y - win.x), 0.0, 1.0);
    let n = P.counts.x;
    let u = t * f32(n - 1u);
    let i = min(u32(floor(u)), n - 2u);
    let f = u - f32(i);
    let base = c * n + i;
    return mix(lut[base], lut[base + 1u], f);
}

@compute @workgroup_size(8, 8, 1)
fn main(@builtin(global_invocation_id) gid: vec3<u32>) {
    let vw = P.counts.z;
    let vh = P.counts.w;
    if (gid.x >= vw || gid.y >= vh) {
        return;
    }
    let pix = gid.y * vw + gid.x;
    let r0 = rays[2u * pix];
    let r1 = rays[2u * pix + 1u];
    let origin = r0.xyz;
    let t0 = r0.w;
    let dir = r1.xyz;
    let t1 = r1.w;
    let step = P.spacing_step.w;
    let channels = P.dims.w;
    let clip_on = P.counts.y != 0u;
    let cutoff = P.cutoff.x;

    var colour = vec3<f32>(0.0);
    var alpha = 0.0;
    var t = t0 + 0.5 * step;
    loop {
        if (t >= t1 || alpha >= cutoff) {
            break;
        }
        let p = origin + dir * t;
        let kept = !clip_on || (dot(P.clip.xyz, p) + P.clip.w >= 0.0);
        if (kept) {
            var weighted = vec3<f32>(0.0);
            var alpha_sum = 0.0;
            var transparency = 1.0;
            for (var c = 0u; c < channels; c = c + 1u) {
                let v = sample_channel(c, p);
                var e = vec4<f32>(0.0);
                // Non-finite samples are transparent, as on the CPU.
                if (abs(v) <= 3.4028235e38) {
                    e = lut_entry(c, v);
                }
                weighted = weighted + e.rgb * e.a;
                alpha_sum = alpha_sum + e.a;
                transparency = transparency * (1.0 - e.a);
            }
            if (alpha_sum > 0.0) {
                let s_colour = weighted / alpha_sum;
                let k = (1.0 - alpha) * (1.0 - transparency);
                colour = colour + k * s_colour;
                alpha = alpha + k;
            }
        }
        t = t + step;
    }
    out_px[pix] = vec4<f32>(colour, alpha);
}
