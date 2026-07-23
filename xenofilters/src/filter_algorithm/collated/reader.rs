//! [`CollatedReader`] extracts complete [`FragmentState`]s from a collated
//! BAM stream -- consuming all records that share a canonical read name
//! before returning.

use crate::Error;
use crate::alignment::{FragmentState, SimpleRec};
use crate::aln_stream::AlignmentStream;
use crate::region::{PositiveRegions, ScoreFn};
use crate::variant::StoreTrait;
use noodles::sam::Header;
use noodles::sam::alignment::record_buf::RecordBuf;
use std::sync::Arc;

pub(crate) struct CollatedReader<R> {
    pub(super) inner: Box<dyn AlignmentStream<R>>,
    peeked: Option<R>,
    species_nr: usize,
    bisulfite: bool,
}

impl<R: SimpleRec> CollatedReader<R> {
    pub(crate) fn new(
        inner: Box<dyn AlignmentStream<R>>,
        bisulfite: bool,
        species_nr: usize,
    ) -> Self {
        Self {
            inner,
            peeked: None,
            species_nr,
            bisulfite,
        }
    }
    pub(crate) fn positive_regions(&self) -> Option<(&PositiveRegions, ScoreFn)> {
        self.inner.positive_regions()
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
            Some(n) => {
                let b: &[u8] = n.as_ref();
                b.to_vec()
            }
            None => {
                return Ok(Some(FragmentState::from_record(
                    first,
                    self.species_nr,
                    self.bisulfite,
                )?));
            }
        };

        let mut state = FragmentState::from_record(first, self.species_nr, self.bisulfite)?;

        loop {
            match self.inner.next_rec()? {
                None => break,
                Some(r) => {
                    let key: &[u8] = match r.name() {
                        Some(n) => n.as_ref(),
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

    pub(crate) fn write_record(
        &mut self,
        rec: RecordBuf,
        state: Option<bool>,
    ) -> Result<(), Error> {
        self.inner.write_record(rec, state)
    }
}
