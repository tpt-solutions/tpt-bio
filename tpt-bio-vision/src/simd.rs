//! CPU feature detection for the `tpt-simd` kernels, with a scalar fallback.
//!
//! `tpt-simd` selects a kernel at run time from the CPU's features (the default
//! `runtime-dispatch` feature). It falls back to scalar code on CPUs without the features it
//! needs. This module reports which features the host has and which backend a build uses.
//! Building with the `scalar-only` feature forces the scalar kernels, so the fallback path can
//! be tested on any machine.

use tpt_simd::detect::Features;

/// Backend names reported by [`backend`].
pub const BACKENDS: [&str; 6] = ["scalar", "sse4.1", "avx2+fma", "avx512", "neon", "sve"];

/// CPU features the host reports at run time.
#[must_use]
pub fn detected_features() -> Features {
    Features::runtime()
}

/// Name of the backend this build uses on the running CPU.
///
/// The name is the highest of the detected ISAs that `tpt-simd` dispatches on. It is `"scalar"`
/// when the crate was built with `scalar-only`, or when the CPU has none of them.
#[must_use]
pub fn backend() -> &'static str {
    if cfg!(feature = "scalar-only") {
        return "scalar";
    }
    let f = detected_features();
    if f.avx512f {
        "avx512"
    } else if f.avx2 && f.fma {
        "avx2+fma"
    } else if f.sse41 {
        "sse4.1"
    } else if f.neon {
        "neon"
    } else if f.sve {
        "sve"
    } else {
        "scalar"
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn simd_backend_is_a_known_name() {
        assert!(BACKENDS.contains(&backend()), "unknown backend {}", backend());
    }

    #[test]
    fn simd_scalar_only_build_reports_scalar() {
        if cfg!(feature = "scalar-only") {
            assert_eq!(backend(), "scalar");
        }
    }
}
