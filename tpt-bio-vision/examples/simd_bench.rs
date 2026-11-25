//! Times a separable Gaussian blur, to compare the SIMD backend with the scalar build.
//!
//! Run `cargo run --release --example simd_bench` for the default build, and
//! `cargo run --release --example simd_bench --features scalar-only` for the scalar
//! kernels. Each build prints the backend name, so the two outputs can be compared.

use std::time::{Duration, Instant};

use tpt_bio_vision::filter::gaussian_blur;
use tpt_bio_vision::simd;

const SIZE: usize = 2048;
const SIGMA: f32 = 2.0;
const RUNS: usize = 7;

fn main() {
    let input: Vec<f32> = (0..SIZE * SIZE).map(|i| ((i * 37) % 251) as f32).collect();
    // One warm-up run, so the first timed run does not pay for page faults.
    let _ = gaussian_blur(&input, SIZE, SIZE, SIGMA).expect("blur");
    let mut times: Vec<Duration> = (0..RUNS)
        .map(|_| {
            let t = Instant::now();
            let out = gaussian_blur(&input, SIZE, SIZE, SIGMA).expect("blur");
            std::hint::black_box(&out);
            t.elapsed()
        })
        .collect();
    times.sort();
    let median = times[RUNS / 2];
    println!(
        "backend={} size={SIZE}x{SIZE} sigma={SIGMA} median={:.2} ms (min {:.2} ms, {RUNS} runs)",
        simd::backend(),
        median.as_secs_f64() * 1e3,
        times[0].as_secs_f64() * 1e3,
    );
}
