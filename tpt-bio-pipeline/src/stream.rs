//! A streaming map with bounded channels and backpressure.
//!
//! [`stream_map`] pulls items from an iterator on a producer thread, applies a
//! function on `workers` threads, and passes results to a sink on the calling
//! thread. Both hand-offs use bounded crossbeam channels of `capacity` slots. When
//! the sink is slow, the channels fill, the workers stop taking items, and the
//! producer blocks. Memory use therefore stays bounded no matter how long the input is.
//!
//! The sink receives results in input order. A small reorder buffer holds results
//! that finish early. The buffer is not bounded by `capacity`, so one very slow item
//! can make it grow while the other workers keep going.

use std::collections::BTreeMap;
use std::sync::atomic::{AtomicBool, Ordering};
use std::thread;

use crossbeam::channel::bounded;

/// Maps `f` over `source` in parallel and feeds results to `sink` in input order.
///
/// The first error returned by `f` stops the run: no new items are started, and the
/// producer and workers exit. Items that were already in flight are discarded. With
/// more than one worker, "first" means the first error to arrive, not the first by
/// input position.
///
/// Returns the number of results passed to `sink` when the run succeeds.
///
/// # Errors
/// Returns the first error from `f`.
///
/// # Panics
/// Panics if `workers` or `capacity` is zero.
pub fn stream_map<I, T, U, E, F, S>(
    source: I,
    workers: usize,
    capacity: usize,
    f: F,
    mut sink: S,
) -> Result<usize, E>
where
    I: IntoIterator<Item = T>,
    I::IntoIter: Send,
    T: Send,
    U: Send,
    E: Send,
    F: Fn(T) -> Result<U, E> + Sync,
    S: FnMut(U),
{
    assert!(workers > 0, "workers must be non-zero");
    assert!(capacity > 0, "capacity must be non-zero");

    let (in_tx, in_rx) = bounded::<(usize, T)>(capacity);
    let (out_tx, out_rx) = bounded::<(usize, Result<U, E>)>(capacity);
    let stop = AtomicBool::new(false);
    let f = &f;
    let stop = &stop;
    let source = source.into_iter();

    thread::scope(|s| {
        s.spawn(move || {
            for (i, item) in source.enumerate() {
                if stop.load(Ordering::Relaxed) || in_tx.send((i, item)).is_err() {
                    break;
                }
            }
        });

        for _ in 0..workers {
            let in_rx = in_rx.clone();
            let out_tx = out_tx.clone();
            s.spawn(move || {
                for (i, item) in in_rx {
                    if stop.load(Ordering::Relaxed) {
                        break;
                    }
                    if out_tx.send((i, f(item))).is_err() {
                        break;
                    }
                }
            });
        }
        drop(in_rx);
        drop(out_tx);

        let mut next = 0usize;
        let mut pending: BTreeMap<usize, U> = BTreeMap::new();
        let mut count = 0usize;
        let mut failure: Option<E> = None;
        while let Ok((i, outcome)) = out_rx.recv() {
            match outcome {
                Ok(value) => {
                    pending.insert(i, value);
                    while let Some(value) = pending.remove(&next) {
                        sink(value);
                        next += 1;
                        count += 1;
                    }
                }
                Err(e) => {
                    failure = Some(e);
                    break;
                }
            }
        }
        // Wake any blocked worker and the producer so the scope can join.
        stop.store(true, Ordering::Relaxed);
        drop(out_rx);
        failure.map_or(Ok(count), Err)
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::sync::atomic::AtomicUsize;
    use std::time::Duration;

    #[test]
    fn results_arrive_in_input_order() {
        let mut out = Vec::new();
        let n = stream_map(0..500u64, 4, 8, |x| Ok::<_, ()>(x * x), |v| out.push(v)).unwrap();
        assert_eq!(n, 500);
        assert_eq!(out, (0..500u64).map(|x| x * x).collect::<Vec<_>>());
    }

    #[test]
    fn empty_input_is_fine() {
        let n = stream_map(std::iter::empty::<u8>(), 2, 2, Ok::<_, ()>, |_| {}).unwrap();
        assert_eq!(n, 0);
    }

    #[test]
    fn first_error_stops_the_run() {
        let produced = AtomicUsize::new(0);
        let source = (0..100_000u32).inspect(|_| {
            produced.fetch_add(1, Ordering::Relaxed);
        });
        let mut sunk = Vec::new();
        let err = stream_map(
            source,
            2,
            4,
            |x| if x == 10 { Err("boom") } else { Ok(x) },
            |v| sunk.push(v),
        )
        .unwrap_err();
        assert_eq!(err, "boom");
        // The sink only ever sees an in-order prefix of the successful results.
        assert!(sunk.len() <= 10);
        assert!(sunk.iter().copied().eq(0..sunk.len() as u32));
        // The run stopped well before the 100k inputs were consumed.
        assert!(produced.load(Ordering::Relaxed) < 100_000);
    }

    #[test]
    fn slow_sink_bounds_the_lead_of_the_producer() {
        // Single worker, capacity 2: at most about 2 + 1 + 2 + 1 items are in flight.
        let produced = AtomicUsize::new(0);
        let consumed = AtomicUsize::new(0);
        let max_lead = AtomicUsize::new(0);
        let source = (0..60u32).inspect(|_| {
            let p = produced.fetch_add(1, Ordering::SeqCst) + 1;
            let lead = p.saturating_sub(consumed.load(Ordering::SeqCst));
            max_lead.fetch_max(lead, Ordering::SeqCst);
        });
        stream_map(source, 1, 2, Ok::<_, ()>, |_| {
            thread::sleep(Duration::from_millis(1));
            consumed.fetch_add(1, Ordering::SeqCst);
        })
        .unwrap();
        let lead = max_lead.load(Ordering::SeqCst);
        assert!(lead <= 8, "producer ran {lead} items ahead of the sink");
        assert_eq!(consumed.load(Ordering::SeqCst), 60);
    }

    #[test]
    #[should_panic(expected = "workers must be non-zero")]
    fn zero_workers_panics() {
        let _ = stream_map(0..1u8, 0, 1, Ok::<_, ()>, |_| {});
    }
}
