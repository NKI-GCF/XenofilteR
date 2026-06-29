
//! [`CollatedReader`] extracts complete [`FragmentState`]s from a collated
//! BAM stream — consuming all records that share a canonical read name
//! before returning.

use crate::alignment::{FragmentState, SimpleRec};
use crate::aln_stream::AlignmentStream;
use crate::config::StripReadSuffix;
use crate::variant::StoreTrait;
use noodles::sam::alignment::record_buf::RecordBuf;
use noodles::sam::Header;
use std::sync::Arc;
use crate::Error;

/// Strip the `/1` or `/2` suffix from `raw` according to `mode`.
pub(crate) fn canonical_name(raw: &[u8], mode: StripReadSuffix) -> Box<[u8]> {
    let stripped = match mode {
        StripReadSuffix::True => {
            if raw.len() >= 2 { &raw[..raw.len() - 2] } else { raw }
        }
        StripReadSuffix::Variable => {
            if raw.ends_with(b"/1") || raw.ends_with(b"/2") {
                &raw[..raw.len() - 2]
            } else {
                raw
            }
        }
        _ => raw, // False / Auto: exact match
    };
    stripped.into()
}

pub(crate) struct CollatedReader<R> {
    pub(super) inner: Box<dyn AlignmentStream<R>>,
    peeked: Option<R>,
    strip: StripReadSuffix,
    species_nr: usize,
}

impl<R: SimpleRec> CollatedReader<R> {
    pub(crate) fn new(
        inner: Box<dyn AlignmentStream<R>>,
        strip: StripReadSuffix,
        species_nr: usize,
    ) -> Self {
        Self { inner, peeked: None, strip, species_nr }
    }

    /// Yield the next complete fragment (all records sharing a canonical name),
    /// or `None` when the stream is exhausted.
    pub(crate) fn next_fragment(&mut self) -> Result<Option<FragmentState<R>>, Error> {
        let first = match self.peeked.take() {
            Some(r) => r,
            None => match self.inner.next_rec()? {
                Some(r) => r,
                None => return Ok(None),
            },
        };

        let first_key = match first.name() {
            Some(n) => canonical_name(n.as_ref(), self.strip),
            None => {
                return Ok(Some(FragmentState::from_record(first, self.species_nr)?));
            }
        };

        let mut state = FragmentState::from_record(first, self.species_nr)?;

        loop {
            match self.inner.next_rec()? {
                None => break,
                Some(r) => {
                    let key = match r.name() {
                        Some(n) => canonical_name(n.as_ref(), self.strip),
                        None => {
                            self.peeked = Some(r);
                            break;
                        }
                    };
                    if key == first_key {
                        state.add_record(r)?;
                    } else {
                        self.peeked = Some(r);
                        break;
                    }
                }
            }
        }

        Ok(Some(state))
    }

    pub(crate) fn header(&self) -> &Header {
        self.inner.header()
    }

    pub(crate) fn variant_store(&self) -> Option<Arc<dyn StoreTrait>> {
        self.inner.variant_store()
    }

    pub(crate) fn write_record(&mut self, rec: RecordBuf, state: Option<bool>) -> Result<(), Error> {
        self.inner.write_record(rec, state)
    }
}
