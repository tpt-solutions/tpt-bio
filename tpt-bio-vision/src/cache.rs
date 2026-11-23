//! A block cache for [`ByteSource`]s, aimed at remote reads.
//!
//! [`BlockCache`] splits the wrapped source into aligned fixed-size blocks and keeps
//! recently used blocks in memory. Reads are served block by block, so a small read
//! fetches at most one whole block from the underlying source. Eviction is least
//! recently used, and the capacity is set in bytes.
//!
//! The lock is not held while the underlying source is read, so two threads that
//! miss on the same block may both fetch it. The second copy is dropped.

use std::collections::{BTreeMap, HashMap};
use std::io;
use std::sync::{Arc, Mutex, MutexGuard, PoisonError};

use crate::source::{self, ByteSource};

/// Default block size: 64 KiB.
pub const DEFAULT_BLOCK_SIZE: usize = 64 * 1024;

/// Hit, miss and eviction counters for a [`BlockCache`].
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct CacheStats {
    /// Block lookups served from memory.
    pub hits: u64,
    /// Block lookups that read the underlying source.
    pub misses: u64,
    /// Blocks dropped to make room.
    pub evictions: u64,
    /// Blocks currently held.
    pub resident_blocks: usize,
}

#[derive(Debug)]
struct Block {
    data: Arc<Vec<u8>>,
    tick: u64,
}

#[derive(Debug, Default)]
struct State {
    blocks: HashMap<u64, Block>,
    /// Maps a recency tick to a block index. The smallest tick is the least recent.
    order: BTreeMap<u64, u64>,
    tick: u64,
    stats: CacheStats,
}

impl State {
    fn next_tick(&mut self) -> u64 {
        self.tick += 1;
        self.tick
    }

    /// Marks an resident block as most recently used and returns its data.
    fn touch(&mut self, index: u64) -> Option<Arc<Vec<u8>>> {
        let tick = self.next_tick();
        let block = self.blocks.get_mut(&index)?;
        self.order.remove(&block.tick);
        block.tick = tick;
        self.order.insert(tick, index);
        Some(Arc::clone(&block.data))
    }
}

/// An LRU cache of aligned blocks over another [`ByteSource`].
#[derive(Debug)]
pub struct BlockCache<S: ByteSource> {
    inner: S,
    block_size: usize,
    max_blocks: usize,
    state: Mutex<State>,
}

impl<S: ByteSource> BlockCache<S> {
    /// Wraps `inner` with blocks of `block_size` bytes and room for about
    /// `capacity_bytes` of cached data. A capacity below one block disables caching:
    /// every lookup then reads the underlying source.
    ///
    /// # Errors
    /// Returns [`io::ErrorKind::InvalidInput`] if `block_size` is zero.
    pub fn new(inner: S, block_size: usize, capacity_bytes: u64) -> io::Result<Self> {
        if block_size == 0 {
            return Err(io::Error::new(io::ErrorKind::InvalidInput, "block size must be non-zero"));
        }
        let max_blocks = usize::try_from(capacity_bytes / block_size as u64).unwrap_or(usize::MAX);
        Ok(Self { inner, block_size, max_blocks, state: Mutex::new(State::default()) })
    }

    /// Hit, miss and eviction counters.
    #[must_use]
    pub fn stats(&self) -> CacheStats {
        let state = self.state();
        CacheStats { resident_blocks: state.blocks.len(), ..state.stats }
    }

    /// The wrapped source.
    #[must_use]
    pub fn inner(&self) -> &S {
        &self.inner
    }

    /// Unwraps the cache and returns the source. Cached blocks are dropped.
    #[must_use]
    pub fn into_inner(self) -> S {
        self.inner
    }

    fn state(&self) -> MutexGuard<'_, State> {
        // The cache holds no invariant a panic could leave half-written, so a
        // poisoned lock is still safe to use.
        self.state.lock().unwrap_or_else(PoisonError::into_inner)
    }

    /// Returns block `index`, from memory when possible.
    fn block(&self, index: u64) -> io::Result<Arc<Vec<u8>>> {
        {
            let mut state = self.state();
            if let Some(data) = state.touch(index) {
                state.stats.hits += 1;
                return Ok(data);
            }
            state.stats.misses += 1;
        }

        let start = index * self.block_size as u64;
        let len = (self.inner.len() - start).min(self.block_size as u64);
        let data = Arc::new(self.inner.read_vec_at(start, len as usize)?);
        if self.max_blocks == 0 {
            return Ok(data);
        }

        let mut state = self.state();
        if let Some(existing) = state.touch(index) {
            // Another thread inserted the block while the lock was released.
            return Ok(existing);
        }
        while state.blocks.len() >= self.max_blocks {
            let Some((_, victim)) = state.order.pop_first() else { break };
            state.blocks.remove(&victim);
            state.stats.evictions += 1;
        }
        let tick = state.next_tick();
        state.order.insert(tick, index);
        state.blocks.insert(index, Block { data: Arc::clone(&data), tick });
        Ok(data)
    }
}

impl<S: ByteSource> ByteSource for BlockCache<S> {
    fn len(&self) -> u64 {
        self.inner.len()
    }

    fn read_exact_at(&self, offset: u64, buf: &mut [u8]) -> io::Result<()> {
        source::check_range(offset, buf.len(), self.inner.len())?;
        if buf.is_empty() {
            return Ok(());
        }
        let end = offset + buf.len() as u64;
        let block_size = self.block_size as u64;
        for index in offset / block_size..=(end - 1) / block_size {
            let block = self.block(index)?;
            let block_start = index * block_size;
            let from = offset.max(block_start);
            let to = end.min(block_start + block.len() as u64);
            let dst = &mut buf[(from - offset) as usize..(to - offset) as usize];
            let src = &block[(from - block_start) as usize..(to - block_start) as usize];
            dst.copy_from_slice(src);
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use std::sync::atomic::{AtomicUsize, Ordering};

    use super::*;

    /// An in-memory source that counts the reads it serves.
    #[derive(Debug)]
    struct CountingSource {
        data: Vec<u8>,
        reads: AtomicUsize,
    }

    impl CountingSource {
        fn new(data: Vec<u8>) -> Self {
            Self { data, reads: AtomicUsize::new(0) }
        }

        fn reads(&self) -> usize {
            self.reads.load(Ordering::SeqCst)
        }
    }

    impl ByteSource for CountingSource {
        fn len(&self) -> u64 {
            self.data.len() as u64
        }

        fn read_exact_at(&self, offset: u64, buf: &mut [u8]) -> io::Result<()> {
            self.reads.fetch_add(1, Ordering::SeqCst);
            self.data.as_slice().read_exact_at(offset, buf)
        }
    }

    fn data(len: usize) -> Vec<u8> {
        (0..len).map(|i| (i % 251) as u8).collect()
    }

    #[test]
    fn repeated_reads_are_served_from_memory() {
        let cache = BlockCache::new(CountingSource::new(data(16)), 4, 1024).unwrap();
        let first = cache.read_vec_at(0, 8).unwrap();
        assert_eq!(first, data(16)[0..8]);
        assert_eq!(cache.inner().reads(), 2);

        let second = cache.read_vec_at(2, 4).unwrap();
        assert_eq!(second, data(16)[2..6]);
        assert_eq!(cache.inner().reads(), 2, "cached blocks must not be re-read");

        let stats = cache.stats();
        assert_eq!((stats.hits, stats.misses, stats.evictions), (2, 2, 0));
        assert_eq!(stats.resident_blocks, 2);
    }

    #[test]
    fn partial_last_block_and_unaligned_reads() {
        let cache = BlockCache::new(CountingSource::new(data(10)), 4, 1024).unwrap();
        assert_eq!(cache.read_vec_at(3, 7).unwrap(), data(10)[3..10]);
        assert_eq!(cache.read_vec_at(8, 2).unwrap(), data(10)[8..10]);
        assert!(cache.read_vec_at(9, 2).is_err());
        assert_eq!(cache.len(), 10);
    }

    #[test]
    fn least_recently_used_block_is_evicted() {
        // Room for two 4-byte blocks.
        let cache = BlockCache::new(CountingSource::new(data(16)), 4, 8).unwrap();
        let mut buf = [0u8; 1];
        cache.read_exact_at(0, &mut buf).unwrap(); // load block 0
        cache.read_exact_at(4, &mut buf).unwrap(); // load block 1
        cache.read_exact_at(0, &mut buf).unwrap(); // hit block 0, making block 1 the oldest
        cache.read_exact_at(8, &mut buf).unwrap(); // load block 2, evicting block 1
        assert_eq!(cache.inner().reads(), 3);
        assert_eq!(cache.stats().evictions, 1);

        cache.read_exact_at(0, &mut buf).unwrap(); // block 0 must still be resident
        assert_eq!(cache.inner().reads(), 3);

        cache.read_exact_at(4, &mut buf).unwrap(); // block 1 was evicted, so this misses
        assert_eq!(cache.inner().reads(), 4);
        assert_eq!(cache.stats().resident_blocks, 2);
        assert!(cache.stats().hits >= 2);
    }

    #[test]
    fn capacity_below_one_block_disables_caching() {
        let cache = BlockCache::new(CountingSource::new(data(8)), 4, 3).unwrap();
        for _ in 0..3 {
            assert_eq!(cache.read_vec_at(0, 4).unwrap(), data(8)[0..4]);
        }
        assert_eq!(cache.inner().reads(), 3);
        let stats = cache.stats();
        assert_eq!((stats.hits, stats.misses, stats.resident_blocks), (0, 3, 0));
    }

    #[test]
    fn zero_block_size_is_rejected() {
        assert!(BlockCache::new(CountingSource::new(data(4)), 0, 16).is_err());
    }

    #[test]
    fn default_block_size_is_64_kib() {
        assert_eq!(DEFAULT_BLOCK_SIZE, 64 * 1024);
        let cache = BlockCache::new(data(1000), DEFAULT_BLOCK_SIZE, 1 << 20).unwrap();
        assert_eq!(cache.read_vec_at(999, 1).unwrap(), [data(1000)[999]]);
        let inner: Vec<u8> = cache.into_inner();
        assert_eq!(inner.len(), 1000);
    }
}
