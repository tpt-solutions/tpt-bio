//! CPU feature reporting for the `tpt-simd` kernels.
//!
//! `tpt-simd` chooses its vector kernels when the crate is compiled, from the target
//! features in effect. The default build targets the baseline `x86_64` ISA, so its kernels
//! are the non-vector fallback even on a CPU with AVX2. [`backend`] reports the kernels
//! the binary was built with. [`detected_features`] reports what the running CPU supports,
//! which is useful for choosing a build, but does not change the kernels.
//!
//! To build with vector kernels, set the target features when compiling, for example
//! `RUSTFLAGS="-C target-feature=+avx,+fma"` (or `-C target-cpu=native`). A binary built
//! that way needs a CPU with those features. The `scalar-only` feature forces the fallback
//! kernels on any build.

use tpt_simd::detect::Features;

/// Backend names reported by [`backend`].
pub const BACKENDS: [&str; 6] = ["scalar", "sse4.1", "avx+fma", "avx512", "neon", "sve"];

/// CPU features the host reports at run time.
#[must_use]
pub fn detected_features() -> Features {
    Features::runtime()
}

/// Name of the kernels this binary was compiled to use.
///
/// `"scalar"` when the crate was built with `scalar-only`, or when the target features
/// enable no vector kernel in `tpt-simd`.
#[must_use]
pub fn backend() -> &'static str {
    if cfg!(feature = "scalar-only") {
        "scalar"
    } else if cfg!(all(target_arch = "x86_64", target_feature = "avx512f")) {
        "avx512"
    } else if cfg!(all(target_arch = "x86_64", target_feature = "avx", target_feature = "fma")) {
        "avx+fma"
    } else if cfg!(all(target_arch = "x86_64", target_feature = "sse4.1")) {
        "sse4.1"
    } else if cfg!(target_arch = "aarch64") {
        "neon"
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
