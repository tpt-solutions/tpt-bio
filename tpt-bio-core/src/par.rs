//! Rayon-based parallel helpers. Available with the `std` feature.

use alloc::vec::Vec;

use rayon::prelude::*;

/// Applies `f` to every element in parallel, preserving input order.
pub fn par_map<T, U, F>(items: &[T], f: F) -> Vec<U>
where
    T: Sync,
    U: Send,
    F: Fn(&T) -> U + Sync + Send,
{
    items.par_iter().map(f).collect()
}

/// Applies `f` to fixed-size chunks in parallel, preserving chunk order.
///
/// # Panics
/// Panics if `chunk_size` is zero.
pub fn par_chunk_map<T, U, F>(items: &[T], chunk_size: usize, f: F) -> Vec<U>
where
    T: Sync,
    U: Send,
    F: Fn(&[T]) -> U + Sync + Send,
{
    assert!(chunk_size > 0, "chunk_size must be non-zero");
    items.par_chunks(chunk_size).map(f).collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn par_map_preserves_order() {
        let input: Vec<u32> = (0..1000).collect();
        let out = par_map(&input, |x| x * 2);
        assert!(out.iter().enumerate().all(|(i, &v)| v == 2 * i as u32));
    }

    #[test]
    fn par_chunk_map_sums() {
        let input: Vec<u64> = (1..=100).collect();
        let sums = par_chunk_map(&input, 10, |c| c.iter().sum::<u64>());
        assert_eq!(sums.len(), 10);
        assert_eq!(sums.iter().sum::<u64>(), 5050);
    }
}
