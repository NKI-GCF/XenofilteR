use super::core::{AlnBuffer, LineByLine};
use crate::alignment::FragmentState;
use crate::alignment::MdCigFlags;
use crate::alignment::SimpleRec;
use crate::filter_algorithm::line_by_line::READ_CT;
use anyhow::{anyhow, ensure, Result};
use noodles::sam::alignment::RecordBuf;
use smallvec::{smallvec, SmallVec};
use std::cmp::{Ord, Ordering};

enum Resolution {
    Ordered(Ordering),
    Scored(f64),
}

pub(super) enum Decision {
    First,
    Last,
    Ambiguous,
    ConfDelta(u8),
    VariantRescued(u8),
}

impl<R: SimpleRec> LineByLine<R> {
    pub(crate) fn process(&mut self) -> Result<()> {
        let mut best: AlnBuffer<R> = smallvec![];

        let mut i = 0;
        while i != self.aln.len() {
            while let Some(rec) = self.aln[i].next_rec()? {
                if self.handle_record_is_fragment_finished(i, rec, &mut best)? {
                    break;
                }
            }
            let mut decision = None;
            if best.len() > 1 {
                decision = match self.resolve(&best)? {
                    Resolution::Ordered(ord) => self.handle_ordering(&mut best, ord)?,
                    Resolution::Scored(delta) => self.apply_delta(&mut best, delta)?,
                };
                assert!(!best.is_empty());
            }
            i += 1;
            if i == self.aln.len() {
                if best.is_empty() {
                    break;
                }
                self.handle_best(&mut best, decision)?;
                i = 0;
            }
        }
        while i > 0 {
            i -= 1;
            self.print_counters(i);
            ensure!(
                self.aln[i].next_rec()?.is_none(),
                "alignment {i} still has reads"
            );
        }
        Ok(())
    }

    // Only ever takes `best` by shared reference — `cmp_perfect`'s MdCigFlags
    // borrow directly from `best`'s elements, so nothing here may request
    // `&mut AlnBuffer<R>` while those are alive.
    fn resolve(&mut self, best: &AlnBuffer<R>) -> Result<Resolution> {
        let mut ord = best[0].partial_cmp(&best[best.len() - 1]);
        if ord.is_none() {
            let (mcfs1, mcfs2) = best[0].cmp_perfect(&best[best.len() - 1], &mut ord)?;
            if ord.is_none() {
                let delta = self.score_delta(best, mcfs1, mcfs2)?;
                return Ok(Resolution::Scored(delta));
            }
        }
        Ok(Resolution::Ordered(
            ord.expect("ord must be Some when not Scored"),
        ))
    }

    fn score_delta<'b>(
        &mut self,
        best: &'b AlnBuffer<R>,
        mcfs1: SmallVec<[MdCigFlags<'b>; READ_CT]>,
        mcfs2: SmallVec<[MdCigFlags<'b>; READ_CT]>,
    ) -> Result<f64> {
        let first = best.first().unwrap();
        let last = best.last().unwrap();
        let score1 = self.score_candidate(first, mcfs1, first.get_nr())?;
        let score2 = self.score_candidate(last, mcfs2, last.get_nr())?;
        Ok(score1 - score2)
    }

    fn apply_delta(&mut self, best: &mut AlnBuffer<R>, mut delta: f64) -> Result<Option<Decision>> {
        match delta {
            d if d > self.ambiguous_log_threshold => self.handle_greater_than(best)?,
            d if d < -self.ambiguous_log_threshold => {
                self.handle_less_than(best)?;
                delta = -delta;
            }
            _ => return Ok(self.add_decision_tag.then_some(Decision::Ambiguous)),
        }
        if self.add_decision_tag {
            let phred = (10.0 * delta / std::f64::consts::LN_10) as u32;
            Ok(Some(Decision::ConfDelta(phred.min(255) as u8)))
        } else {
            Ok(None)
        }
    }

    fn handle_ordering(
        &mut self,
        best: &mut AlnBuffer<R>,
        ord: Ordering,
    ) -> Result<Option<Decision>> {
        match ord {
            Ordering::Greater => self
                .handle_greater_than(best)
                .map(|()| self.add_decision_tag.then_some(Decision::First)),
            Ordering::Less => self
                .handle_less_than(best)
                .map(|()| self.add_decision_tag.then_some(Decision::Last)),
            Ordering::Equal => Ok(self.add_decision_tag.then_some(Decision::Ambiguous)),
        }
    }
    fn handle_record_is_fragment_finished(
        &mut self,
        i: usize,
        rec: R,
        best: &mut AlnBuffer<R>,
    ) -> Result<bool> {
        if !(self.is_secondary_skipped)(&rec)? {
            let name = rec
                .name()
                .ok_or_else(|| anyhow!("Record has no read name"))?;
            if let Some(new_readname) = (self.is_new_qname)(best, name.as_ref()) {
                if new_readname {
                    #[cfg(test)]
                    if self.aln.is_empty() {
                        return Ok(true);
                    }
                    // end of round for this alignment
                    self.aln[i].un_next(rec)?;
                    return Ok(true);
                }
                for state in best.iter_mut().rev() {
                    if state.get_nr() == i {
                        state.add_record(rec)?;
                        return Ok(false);
                    }
                }
            }
            best.push(FragmentState::from_record(rec, i)?);
        } // else skip secondary

        Ok(false)
    }
    fn handle_greater_than(&mut self, best: &mut AlnBuffer<R>) -> Result<()> {
        let mut last = best.pop().unwrap();
        let nr = last.get_nr();
        last.drain_records().try_for_each(|r| {
            self.write_record(nr, r.as_record_buf(self.aln[nr].header())?, Some(false))
        })
    }
    fn handle_less_than(&mut self, best: &mut AlnBuffer<R>) -> Result<()> {
        let all_before_last = best.len() - 1;
        best.drain(0..all_before_last).try_for_each(|mut b| {
            let nr = b.get_nr();
            b.drain_records().try_for_each(|r| {
                self.write_record(nr, r.as_record_buf(self.aln[nr].header())?, Some(false))
            })
        })
    }
    fn handle_best(&mut self, best: &mut AlnBuffer<R>, decision: Option<Decision>) -> Result<()> {
        let best_state = (best.len() == 1).then_some(true);

        best.drain(..).try_for_each(|mut b| {
            let nr = b.get_nr();
            b.drain_records().try_for_each(|r| {
                if best_state.is_none() && (self.is_unmapped_skipped)(&r)? {
                    self.branch_counters[24 + nr] += 1;
                    return Ok(());
                }
                let header = self
                    .aln
                    .get(nr)
                    .map(|a| a.header())
                    .ok_or_else(|| anyhow!("No alignment for index {nr}"))?;
                let mut record_buf = RecordBuf::try_from_alignment_record(header, &r)?;
                match decision {
                    Some(Decision::ConfDelta(decision)) => {
                        self.add_aux_tags(&mut record_buf, b"XF", decision)?
                    }
                    Some(Decision::VariantRescued(decision)) => {
                        self.add_aux_tags(&mut record_buf, b"XR", decision)?
                    }
                    _ => { /* no tag to add */ }
                }
                self.write_record(nr, record_buf, best_state)
            })
        })
    }
}

#[cfg(test)]
mod tests;
