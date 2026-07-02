//!
//! Design principle: always accumulate `RecordKind`s during pass 1, with the
//! cheapest representation that fits the record's role:
//!
//! - `Secondary`        — virtual offset + flags only; never scored, tags along.
//! - `UnmappedPrimary`   — virtual offset + flags + qualities; never NW-scored
//!                          (no CIGAR/MD exists), qualities reserved for future
//!                          quality-weighted ambiguous tie-breaking.
//! - `Mapped`            — boxed full record (CIGAR/MD/qualities/etc.) for
//!                          primary or supplementary mapped alignments; the
//!                          only variant ever passed to NW scoring.
//!
//! Boxing the expensive variant keeps `RecordKind` small (pointer + tag) so
//! `SmallVec<[RecordKind; 2]>` stays cheap even when most fragments are
//! perfect/unmapped and never reach the boxed branch.
//!
//! Early-assignability is evaluated at fragment *completion* time (see
//! `StreamAccumulator::classify`). When neither fast path applies, a
//! per-mate `MateKind` summary is computed alongside the records so that
//! `resolve_fragment` can skip NW scoring for individual mates that are
//! identical-contribution in both streams (Unmapped-Unmapped or
//! Perfect-Perfect), even when the fragment as a whole requires scoring.

pub(crate) mod pending;
pub(crate) mod stream;

pub(crate) use pending::{PendingFragment, insert};
pub(crate) use stream::{EarlyKind, MappedRecord, RecordKind, StreamKind};

use std::collections::HashMap;

pub(crate) type FragmentTable = HashMap<Box<[u8]>, PendingFragment>;
