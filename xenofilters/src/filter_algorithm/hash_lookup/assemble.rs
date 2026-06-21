//! In-memory fragment accumulator keyed by canonical read name.
//!
//! `NameTable` accumulates records from both streams until a fragment is
//! "complete" — i.e. both streams have contributed at least one primary
//! alignment. Supplementary and secondary records are buffered alongside
//! primaries and follow whatever classification decision is made.

use crate::alignment::{FragmentState, SimpleRec};
use crate::config::StripReadSuffix;
use crate::filter_algorithm::collated::reader::canonical_name;
use smallvec::SmallVec;
use std::collections::HashMap;

/// Per-stream record accumulator inside a pending fragment.
#[derive(Default)]
struct StreamBuf<R> {
    records: SmallVec<[R; 2]>,
    has_primary: bool,
}

impl<R: SimpleRec> StreamBuf<R> {
    fn push(&mut self, r: R) {
        if !r.flags().map(|f| f.is_secondary() || f.is_supplementary()).unwrap_or(false) {
            self.has_primary = true;
        }
        self.records.push(r);
    }

    fn is_ready(&self) -> bool {
        self.has_primary
    }

    fn drain_into_fragment(&mut self, nr: usize) -> Result<FragmentState<R>, anyhow::Error> {
        let mut iter = self.records.drain(..);
        let first = iter.next().ok_or_else(|| anyhow::anyhow!("empty StreamBuf"))?;
        let mut state = FragmentState::from_record(first, nr)?;
        for r in iter {
            state.add_record(r)?;
        }
        Ok(state)
    }
}

/// A fragment in the process of being assembled from both streams.
pub(crate) struct PendingFragment<R> {
    pub(crate) driving: StreamBuf<R>,   // stream 0
    pub(crate) lookup: StreamBuf<R>,    // stream 1
    /// Insertion-order sequence number for output staging.
    pub(crate) seq_nr: u64,
}

impl<R: SimpleRec> PendingFragment<R> {
    fn new(seq_nr: u64) -> Self {
        Self {
            driving: StreamBuf::default(),
            lookup: StreamBuf::default(),
            seq_nr,
        }
    }

    pub(crate) fn push(&mut self, r: R, stream_nr: usize) {
        if stream_nr == 0 {
            self.driving.push(r);
        } else {
            self.lookup.push(r);
        }
    }

    /// True once both streams have at least one primary alignment.
    pub(crate) fn is_complete(&self) -> bool {
        self.driving.is_ready() && self.lookup.is_ready()
    }

    /// Destructure into a pair of `FragmentState`s. Panics if not complete.
    pub(crate) fn into_pair(mut self) -> Result<(FragmentState<R>, FragmentState<R>), anyhow::Error> {
        let a = self.driving.drain_into_fragment(0)?;
        let b = self.lookup.drain_into_fragment(1)?;
        Ok((a, b))
    }

    /// Drain stream-0 records as a fragment (for unmatched emission).
    pub(crate) fn into_driving_fragment(mut self) -> Result<FragmentState<R>, anyhow::Error> {
        self.driving.drain_into_fragment(0)
    }
}

pub(crate) type NameTable<R> = HashMap<Box<[u8]>, PendingFragment<R>>;

/// Insert `rec` from stream `nr` into `table`. Returns the key and `true`
/// if the fragment is now complete and should be extracted for scoring.
pub(crate) fn insert<R: SimpleRec>(
    table: &mut NameTable<R>,
    rec: R,
    nr: usize,
    strip: StripReadSuffix,
    seq_counter: &mut u64,
) -> (Box<[u8]>, bool) {
    let key: Box<[u8]> = match rec.name() {
        Some(n) => canonical_name(n.as_ref(), strip),
        None => {
            // Unnamed records get a unique synthetic key so they don't collide.
            let k: Box<[u8]> = format!("__unnamed_{}", *seq_counter).into_bytes().into();
            *seq_counter += 1;
            k
        }
    };

    let entry = table.entry(key.clone()).or_insert_with(|| {
        let sn = *seq_counter;
        *seq_counter += 1;
        PendingFragment::new(sn)
    });

    entry.push(rec, nr);
    let complete = entry.is_complete();
    (key, complete)
}
