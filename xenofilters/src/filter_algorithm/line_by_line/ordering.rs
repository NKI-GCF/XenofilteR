use super::core::{AlnBuffer, LineByLine};
use crate::alignment::FragmentState;
use anyhow::{Result, anyhow, ensure};
use crate::alignment::SimpleRec;
use noodles::sam::alignment::RecordBuf;
use std::cmp::{Ord, Ordering};
use smallvec::smallvec;

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
                let last = &best[best.len() - 1];

                let mut ord = best[0].partial_cmp(last);

                if ord.is_none() {
                    ord = best[0].cmp_perfect(last)?;
                }
                decision = self.handle_ordering(&mut best, ord)?;
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
    fn handle_record_is_fragment_finished(
        &mut self,
        i: usize,
        rec: R,
        best: &mut AlnBuffer<R>,
    ) -> Result<bool> {
        if !(self.is_secondary_skipped)(&rec)? {
            let name = rec.name().ok_or_else(|| anyhow!("Record has no read name"))?;
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
        last.drain_records()
            .try_for_each(|r| self.write_record(nr, r.as_record_buf(self.aln[nr].header())?, Some(false)))
    }
    fn handle_less_than(&mut self, best: &mut AlnBuffer<R>) -> Result<()> {
        let all_before_last = best.len() - 1;
        best.drain(0..all_before_last).try_for_each(|mut b| {
            let nr = b.get_nr();
            b.drain_records()
                .try_for_each(|r| self.write_record(nr, r.as_record_buf(self.aln[nr].header())?, Some(false)))
        })
    }
    fn handle_ordering(
        &mut self,
        best: &mut AlnBuffer<R>,
        ord: Option<Ordering>,
    ) -> Result<Option<Decision>> {
        match ord {
            Some(Ordering::Greater) => self
                .handle_greater_than(best)
                .map(|()| self.add_decision_tag.then_some(Decision::First)),
            Some(Ordering::Less) => self
                .handle_less_than(best)
                .map(|()| self.add_decision_tag.then_some(Decision::Last)),
            Some(Ordering::Equal) => Ok(self.add_decision_tag.then_some(Decision::Ambiguous)),
            None => {
                // None of the alignments were fully unmapped or perfect matches,
                // so we need to score them to find the best.
                assert!(best.len() >= 2, "score path requires at least two candidates");
                let (first_score, last_score) = {
                    let first = best.first().unwrap();
                    let last  = best.last().unwrap();
                    (self.score_candidate(first, first.get_nr())?,
                     self.score_candidate(last,  last.get_nr())?)
                };

                let mut delta = first_score - last_score;
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
        }
    }
    fn handle_best(
        &mut self,
        best: &mut AlnBuffer<R>,
        decision: Option<Decision>,
    ) -> Result<()> {
        let best_state = (best.len() == 1).then_some(true);

        best.drain(..).try_for_each(|mut b| {
            let nr = b.get_nr();
            b.drain_records().try_for_each(|r| {
                if best_state.is_none() && (self.is_unmapped_skipped)(&r)? {
                    self.branch_counters[24 + nr] += 1;
                    return Ok(());
                }
                let header = self.aln.get(nr).map(|a| a.header()).ok_or_else(|| anyhow!("No alignment for index {nr}"))?;
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

#[cfg(test)]
fn debug_print_best<R: SimpleRec>(best: &AlnBuffer<R>, last: &AlnBuffer<R>, ord: Option<std::cmp::Ordering>) {
    // FIXME: this does not print or test all reads in the buffer, just the first one.
    let best_rec = &best[0].get_records()[0];
    let last_rec = &last[0].get_records()[0];
    assert_eq!(best_rec.name(), last_rec.name());
    eprintln!(
        "{}: {} vs {} => {:?}",
        std::str::from_utf8(best_rec.name().as_ref().unwrap()).unwrap_or("<?>"),
        best[0].get_nr(),
        best.last().unwrap().get_nr(),
        ord
    );
}
