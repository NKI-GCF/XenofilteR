use anyhow::{Result, ensure};
use rust_htslib::bam::record::Record;
use smallvec::{SmallVec, smallvec};
use crate::{Options, Penalties};
use crate::aln_stream::AlnStream;
use crate::fragment::FragmentState;
use crate::alignment::stitched_fragment;
use std::cmp::Ordering;
use std::sync::OnceLock;

// Usually two species compared so two alignments, two fragment states.

type QnameCompareFn = fn(&AlnBuffer, &[u8]) -> Option<bool>;
type AlnBuffer = SmallVec<[FragmentState; 2]>;

pub struct LineByLine {
    aln: SmallVec<[AlnStream; 2]>,
    branch_counters: [u64; 32],
    is_secondary_skipped: fn(&Record) -> bool,
    is_unmapped_skipped: fn(&Record) -> bool,
    compare_qname: fn($AlnBuffer, &[u8]) -> bool,
    penalties: Penalties
}

impl LineByLine {
    pub fn new(config: Config, aln: SmallVec<[AlnStream; 2]>) -> Self {
        let is_unmapped_skipped = match self.discard_unmapped {
            true => |rec: &Record| rec.is_unmapped() && rec.is_mate_unmapped(),
            false => |_| false,
        };

        let is_secondary_skipped = match self.skip_secondary {
            true => |rec: &Record| rec.is_secondary(),
            false => |_| false,
        };
        let compare_qname = match self.strip_read_suffix {
            Some(true) => |best: &AlnBuffer, qname2: &[u8]| {
                best.first()
                    .map(|b| b.first_qname())
                    .map(|qname1| qname1[..qname1.len() - 2] != qname2[..qname2.len() - 2])
            },
            Some(false) => |best: &AlnBuffer, qname2: &[u8]| {
                best.first()
                    .map(|b| b.first_qname())
                    .map(|qname1| qname1 != qname2)
            },
            None => |best: &AlnBuffer, qname2: &[u8]| {
                best.first().map(|b| b.first_qname()).map(|qname1| {
                    if qname1.ends_with(b"/1") || qname1.ends_with(b"/2") {
                        qname1[..qname1.len() - 2] != qname2[..qname2.len() - 2]
                    } else {
                        qname1 != qname2
                    }
               })
           },
        };
        LineByLine {
            aln,
            branch_counters: [0; 32],
            is_secondary_skipped,
            is_unmapped_skipped,
            compare_qname,
            penalties: config.to_penalties()
        }
    }

    fn write_record(&mut self, i: usize, rec: Record, best_state: Option<bool>) -> Result<()> {
        if self.is_skipped_unmapped(&rec) {
            return Ok(());
        }
        match (i, best_state) {
            (i, Some(false)) => self.branch_counters[i << 1] += 1,
            (i, Some(true)) => self.branch_counters[1 + (i << 1)] += 1,
            (i, None) => self.branch_counters[16 + i] += 1,
        }
        self.aln[i].write_record(rec, best_state)
    }

    fn handle_ordering(&mut self, best: &mut AlnBuffer, ord: Option<Ordering>) -> Result<()> {
        match ord {
            Some(Ordering::Greater) => {
                let all_before_last = best.len() - 1;
                best.drain(0..all_before_last).try_for_each(|mut b| {
                    let nr = b.get_nr();
                    b.records
                        .drain(..)
                        .try_for_each(|r| self.write_record(nr, r, None))
                })
            }
            Some(Ordering::Less) => {
                let mut last = best.pop().unwrap();
                let nr = last.get_nr();
                last.records
                    .drain(..)
                    .try_for_each(|r| self.write_record(nr, r, None))
            }
            Some(Ordering::Equal) => {
                // Neither is written. (FIXME: make configurable?)
                Ok(())
            }
            None => {
                // None of the alignments were fully unmapped or perfect matches,
                // so we need to score them to find the best.
                let first = &best.first().unwrap();
                let first_ord = first.order_mates();

                let last = &best.last().unwrap();
                let last_ord = last.order_mates();

                let first_score = stitched_fragment(&self.penalties, &first.records, first_ord)?.score()?;
                let last_score = stitched_fragment(&self.penalties, &last.records, last_ord)?.score()?;

                let mut ord = first_score.partial_cmp(&last_score);
                if ord.is_none() {
                    ord = Some(Ordering::Equal);
                }
                self.handle_ordering(best, ord)
            }
        }
    }

    fn handle_record_is_fragment_finished(
        &mut self,
        i: usize,
        rec: Record,
        best: &mut AlnBuffer,
    ) -> Result<bool> {
        if !self.is_skipped_secondary(&rec) {
            if let Some(new_readname) = self.ompare_qname(best, rec.qname()) {
                if new_readname {
                    // end of round for this alignment
                    self.aln[i].un_next(rec);
                    // FIXME: comparing more than 2 alignments?
                    if best.len() > 1 {
                        let ord = best.last().unwrap().partial_cmp(&best[0]);
                        self.handle_ordering(best, ord)?;
                    }
                    return Ok(true);
                }
                for state in best.iter_mut().rev() {
                    if state.get_nr() == i {
                        state.records.push(rec);
                        return Ok(false);
                    }
                }
            }
            best.push(FragmentState::from_record(rec, i));
        } // else skip secondary

        Ok(false)
    }

    fn handle_best(&mut self, best: &mut AlnBuffer) -> Result<()> {
        if best.len() > 1 {
            best.sort_unstable_by(|a, b| {
                b.partial_cmp(&a).unwrap_or(std::cmp::Ordering::Equal)
            });
            if best[0] > best[1] {
                best.drain(1..).try_for_each(|mut b| {
                    let nr = b.get_nr();
                    b.records
                        .drain(..)
                        .try_for_each(|r| self.write_record(nr, r, None))
                })?;
            }
        }
        let best_state = (best.len() == 1).then_some(true);

        best.drain(..).try_for_each(|mut b| {
            let nr = b.get_nr();
            b.records
                .drain(..)
                .try_for_each(|r| self.write_record(nr, r, best_state))
        })
    }

    pub fn process(&mut self) -> Result<()> {
        let mut best: AlnBuffer = smallvec![];

        let mut i = 0;
        while i != self.aln.len() {
            while let Some(rec) = self.aln[i].next_rec()? {
                if self.handle_record_is_fragment_finished(i, rec, &mut best)? {
                    i += 1;
                    if i == self.aln.len() {
                        self.handle_best(&mut best)?;
                        i = 0;
                    }
                }
            }
            i += 1;
        }
        self.handle_best(&mut best)?;
        while i > 0 {
            i -= 1;
            eprintln!(
                "Filtered from alignment {i}: {}",
                self.branch_counters[i << 1]
            );
            eprintln!(
                "Assigned to alignment {i}: {}",
                self.branch_counters[1 + (i << 1)]
            );
            eprintln!(
                "Ambiguous for alignment {i}: {}",
                self.branch_counters[16 + i]
            );
            ensure!(
                self.aln[i].next_rec()?.is_none(),
                "alignment {i} still has reads"
            );
        }
        Ok(())
    }
}
