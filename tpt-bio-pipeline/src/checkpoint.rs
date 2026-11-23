//! Checkpointing and resume for [`TypedDag`].
//!
//! When a node finishes, its output is encoded with [`Codec`] and written to
//! `<dir>/<node name>.ckpt`. On a later run with the same store, a node whose file
//! exists and decodes is skipped, and its stored value is used as its result.
//!
//! The format is whatever the [`Codec`] implementation writes. Files are written to
//! a temporary name and renamed into place, so a crash never leaves a partial
//! checkpoint under the final name.
//!
//! Checkpoints are keyed only by node name. Changing a node's function or its
//! input data does not invalidate its old checkpoint. Delete the directory, or call
//! [`CheckpointStore::remove`], to force a node to rerun.

use std::fmt;
use std::fs;
use std::io;
use std::path::{Path, PathBuf};
use std::sync::atomic::AtomicBool;

use crate::typed::{TypedDag, is_valid_name};
use crate::{ExecError, NodeId, Progress};

/// Serialization for values stored in a checkpoint.
pub trait Codec: Sized {
    /// Encodes the value as bytes.
    fn encode(&self) -> Vec<u8>;

    /// Decodes bytes written by [`Codec::encode`]. Returns `None` for anything else,
    /// in which case the node is rerun.
    fn decode(bytes: &[u8]) -> Option<Self>;
}

impl Codec for u64 {
    fn encode(&self) -> Vec<u8> {
        self.to_le_bytes().to_vec()
    }

    fn decode(bytes: &[u8]) -> Option<Self> {
        <[u8; 8]>::try_from(bytes).ok().map(u64::from_le_bytes)
    }
}

impl Codec for Vec<u8> {
    fn encode(&self) -> Vec<u8> {
        self.clone()
    }

    fn decode(bytes: &[u8]) -> Option<Self> {
        Some(bytes.to_vec())
    }
}

impl Codec for String {
    fn encode(&self) -> Vec<u8> {
        self.as_bytes().to_vec()
    }

    fn decode(bytes: &[u8]) -> Option<Self> {
        String::from_utf8(bytes.to_vec()).ok()
    }
}

/// A node failure, or a failure to read or write a checkpoint.
#[derive(Debug)]
pub enum CheckpointError<E> {
    /// The node's own function returned an error.
    Node(E),
    /// Reading or writing a checkpoint file failed.
    Io(io::Error),
}

impl<E: fmt::Display> fmt::Display for CheckpointError<E> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Node(e) => e.fmt(f),
            Self::Io(e) => write!(f, "checkpoint I/O failed: {e}"),
        }
    }
}

impl<E: fmt::Debug + fmt::Display> std::error::Error for CheckpointError<E> {}

/// A directory that holds one checkpoint file per node.
#[derive(Debug, Clone)]
pub struct CheckpointStore {
    dir: PathBuf,
}

impl CheckpointStore {
    /// Opens a store, creating the directory if needed.
    ///
    /// # Errors
    /// Returns the I/O error from creating the directory.
    pub fn open(dir: impl Into<PathBuf>) -> io::Result<Self> {
        let dir = dir.into();
        fs::create_dir_all(&dir)?;
        Ok(Self { dir })
    }

    /// The directory that holds the checkpoint files.
    #[must_use]
    pub fn dir(&self) -> &Path {
        &self.dir
    }

    /// The stored bytes for `name`, or `None` if there is no checkpoint.
    ///
    /// # Errors
    /// Returns an error for an invalid name or an I/O failure other than "not found".
    pub fn load(&self, name: &str) -> io::Result<Option<Vec<u8>>> {
        match fs::read(self.path_for(name)?) {
            Ok(bytes) => Ok(Some(bytes)),
            Err(e) if e.kind() == io::ErrorKind::NotFound => Ok(None),
            Err(e) => Err(e),
        }
    }

    /// Writes the checkpoint for `name`, replacing any previous one.
    ///
    /// # Errors
    /// Returns an error for an invalid name or an I/O failure.
    pub fn save(&self, name: &str, bytes: &[u8]) -> io::Result<()> {
        let path = self.path_for(name)?;
        let tmp = path.with_extension("ckpt.tmp");
        fs::write(&tmp, bytes)?;
        fs::rename(&tmp, &path)
    }

    /// Deletes the checkpoint for `name`. Returns `false` if there was none.
    ///
    /// # Errors
    /// Returns an error for an invalid name or an I/O failure other than "not found".
    pub fn remove(&self, name: &str) -> io::Result<bool> {
        match fs::remove_file(self.path_for(name)?) {
            Ok(()) => Ok(true),
            Err(e) if e.kind() == io::ErrorKind::NotFound => Ok(false),
            Err(e) => Err(e),
        }
    }

    fn path_for(&self, name: &str) -> io::Result<PathBuf> {
        if !is_valid_name(name) {
            return Err(io::Error::new(io::ErrorKind::InvalidInput, "invalid checkpoint name"));
        }
        Ok(self.dir.join(format!("{name}.ckpt")))
    }
}

impl<T: Sync> TypedDag<T> {
    /// Runs the graph in parallel levels, like [`crate::Dag::execute_parallel`], but
    /// loads each node's output from `store` when a valid checkpoint exists and
    /// saves it after the node runs.
    ///
    /// Skipped nodes still produce `started` and `finished` progress events.
    ///
    /// # Errors
    /// Returns [`ExecError::Node`] with [`CheckpointError::Node`] when a node fails,
    /// [`ExecError::Node`] with [`CheckpointError::Io`] when a checkpoint cannot be
    /// written, and [`ExecError::Cancelled`] when `cancel` is set.
    pub fn execute_checkpointed<R, E, F, P>(
        &self,
        cancel: &AtomicBool,
        store: &CheckpointStore,
        progress: &P,
        f: F,
    ) -> Result<Vec<R>, ExecError<CheckpointError<E>>>
    where
        R: Codec + Send + Sync,
        E: Send,
        F: Fn(NodeId, &T, &[&R]) -> Result<R, E> + Sync,
        P: Progress,
    {
        self.dag().execute_parallel_observed(cancel, progress, |n, payload, inputs| {
            let name = self.name(n).unwrap_or_default();
            if let Some(bytes) = store.load(name).map_err(CheckpointError::Io)?
                && let Some(value) = R::decode(&bytes)
            {
                return Ok(value);
            }
            let value = f(n, payload, inputs).map_err(CheckpointError::Node)?;
            store.save(name, &value.encode()).map_err(CheckpointError::Io)?;
            Ok(value)
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::TypeTag;
    use std::convert::Infallible;
    use std::sync::atomic::{AtomicBool, AtomicUsize, Ordering};
    use std::time::{SystemTime, UNIX_EPOCH};

    const NUM: TypeTag = TypeTag::new("num");

    /// A scratch directory under the system temp dir, removed on drop.
    struct ScratchDir(PathBuf);

    impl ScratchDir {
        fn new(tag: &str) -> Self {
            let nanos = SystemTime::now().duration_since(UNIX_EPOCH).unwrap().as_nanos();
            let dir = std::env::temp_dir()
                .join(format!("tpt-bio-pipeline-{tag}-{}-{nanos}", std::process::id()));
            Self(dir)
        }
    }

    impl Drop for ScratchDir {
        fn drop(&mut self) {
            let _ = fs::remove_dir_all(&self.0);
        }
    }

    /// src(7) -> double -> add <- src2(1)  where add = sum of inputs, double = 2 * input.
    fn graph() -> TypedDag<u64> {
        let mut g = TypedDag::new();
        let src = g.add_node("src", &[], NUM, 7).unwrap();
        let src2 = g.add_node("src2", &[], NUM, 1).unwrap();
        let double = g.add_node("double", &[NUM], NUM, 0).unwrap();
        let add = g.add_node("add", &[NUM, NUM], NUM, 0).unwrap();
        g.add_edge(src, double).unwrap();
        g.add_edge(double, add).unwrap();
        g.add_edge(src2, add).unwrap();
        g.validate().unwrap();
        g
    }

    fn run(g: &TypedDag<u64>, store: &CheckpointStore, calls: &AtomicUsize) -> Vec<u64> {
        let cancel = AtomicBool::new(false);
        g.execute_checkpointed(&cancel, store, &crate::NoProgress, |n, payload, inputs| {
            calls.fetch_add(1, Ordering::SeqCst);
            let sum: u64 = inputs.iter().map(|x| **x).sum();
            // Node 2 doubles its single input; node 3 adds its two inputs.
            Ok::<u64, Infallible>(match g.name(n) {
                Some("double") => 2 * sum,
                _ => *payload + sum,
            })
        })
        .unwrap()
    }

    #[test]
    fn first_run_saves_every_node_then_rerun_skips_them() {
        let scratch = ScratchDir::new("ckpt-skip");
        let store = CheckpointStore::open(&scratch.0).unwrap();
        let g = graph();

        let calls = AtomicUsize::new(0);
        let first = run(&g, &store, &calls);
        assert_eq!(calls.load(Ordering::SeqCst), 4);
        // src=7, src2=1, double=14, add=14+1=15.
        assert_eq!(first, vec![7, 1, 14, 15]);
        for name in ["src", "src2", "double", "add"] {
            assert!(store.load(name).unwrap().is_some(), "{name} not saved");
        }

        let calls = AtomicUsize::new(0);
        let second = run(&g, &store, &calls);
        assert_eq!(calls.load(Ordering::SeqCst), 0, "completed nodes must be skipped");
        assert_eq!(second, first);
    }

    #[test]
    fn only_missing_or_corrupt_nodes_rerun() {
        let scratch = ScratchDir::new("ckpt-partial");
        let store = CheckpointStore::open(&scratch.0).unwrap();
        let g = graph();
        let calls = AtomicUsize::new(0);
        let expected = run(&g, &store, &calls);

        assert!(store.remove("double").unwrap());
        // Five bytes cannot decode as u64, so `src` must rerun too.
        store.save("src", &[1, 2, 3, 4, 5]).unwrap();

        let calls = AtomicUsize::new(0);
        let resumed = run(&g, &store, &calls);
        assert_eq!(calls.load(Ordering::SeqCst), 2);
        assert_eq!(resumed, expected);
        assert_eq!(store.load("src").unwrap(), Some(7u64.encode()));
    }

    #[test]
    fn failed_node_is_not_checkpointed() {
        let scratch = ScratchDir::new("ckpt-fail");
        let store = CheckpointStore::open(&scratch.0).unwrap();
        let g = graph();
        let cancel = AtomicBool::new(false);
        let err = g
            .execute_checkpointed(&cancel, &store, &crate::NoProgress, |n, payload, _: &[&u64]| {
                if g.name(n) == Some("double") { Err("bad input") } else { Ok(*payload) }
            })
            .unwrap_err();
        assert!(matches!(err, ExecError::Node { error: CheckpointError::Node("bad input"), .. }));
        assert!(store.load("double").unwrap().is_none());
        assert!(store.load("src").unwrap().is_some());
    }

    #[test]
    fn names_are_checked_before_touching_the_filesystem() {
        let scratch = ScratchDir::new("ckpt-names");
        let store = CheckpointStore::open(&scratch.0).unwrap();
        assert_eq!(store.load("../escape").unwrap_err().kind(), io::ErrorKind::InvalidInput);
        assert!(store.save("", &[]).is_err());
    }

    #[test]
    fn codecs_round_trip() {
        assert_eq!(u64::decode(&42u64.encode()), Some(42));
        assert_eq!(u64::decode(&[1, 2]), None);
        assert_eq!(String::decode(&"héllo".to_string().encode()).as_deref(), Some("héllo"));
        assert_eq!(String::decode(&[0xFF]), None);
        assert_eq!(Vec::<u8>::decode(&vec![9u8, 8].encode()), Some(vec![9, 8]));
    }
}
