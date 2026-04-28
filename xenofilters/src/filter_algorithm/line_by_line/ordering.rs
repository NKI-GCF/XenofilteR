use super::core::{AlnBuffer, LineByLine};
use crate::alignment::{build_fragment, stringify_record};
use crate::fragment_state::FragmentState;
use anyhow::{Result, anyhow};
use rust_htslib::bam::record::{Aux, Record};
use std::cmp::Ordering;
use smallvec::{SmallVec, smallvec};
use crate::variant::VntPerRec;

pub(crate) enum Decision {
    First,
    Last,
    Ambiguous,
    ConfDelta(u8),
    VariantRescued(u8),
}

fn alignment_range(rec: &Record) -> Option<(i32, i64, i64)> {
    if rec.is_unmapped() {
        None
    } else {
        Some((rec.tid(), rec.pos(), rec.cigar().end_pos()))
    }
}

impl LineByLine {
    fn handle_greater_than(&mut self, best: &mut AlnBuffer) -> Result<()> {
        let mut last = best.pop().unwrap();
        let nr = last.get_nr();
        last.records
            .drain(..)
            .try_for_each(|r| self.write_record(nr, r, Some(false)))
    }
    fn handle_less_than(&mut self, best: &mut AlnBuffer) -> Result<()> {
        let all_before_last = best.len() - 1;
        best.drain(0..all_before_last).try_for_each(|mut b| {
            let nr = b.get_nr();
            b.records
                .drain(..)
                .try_for_each(|r| self.write_record(nr, r, Some(false)))
        })
    }
    fn score_candidate(
        &self,
        state: &FragmentState,
        aln_idx: usize,
    ) -> Result<f64> {
        let mut ranges: SmallVec<[(usize, i64, i64); 2]> = SmallVec::with_capacity(state.records.len());
        for rec in &state.records {
            if let Some((tid, start, end)) = alignment_range(rec) {
                ranges.push((tid as usize, start, end));
            }
        }

        let mut dvnt_per_rec = SmallVec::with_capacity(state.records.len());
        for rec in &state.records {
            if let Some((tid, start, end)) = alignment_range(rec) {
                let ranges: SmallVec<[(usize, i64, i64); 1]> = smallvec![(tid as usize, start, end)];
                let delta_vars = self.aln.get(aln_idx)
                    .and_then(|a| a.variant_store())
                    .map(|s| s.variants_overlapping_multi(&ranges))
                    .unwrap_or_default();
                dvnt_per_rec.push(delta_vars);
            } else {
                dvnt_per_rec.push(SmallVec::new());
            }
        }

        let fragment = build_fragment(
            &state.records,
            state.order_mates()
        )?;

        fragment.score(&self.penalties, dvnt_per_rec).map_err(|e| {
            anyhow!(
                "Error scoring fragment for alignment {aln_idx}: {}\n{}",
                e,
                state
                    .records
                    .iter()
                    .map(stringify_record)
                    .collect::<Vec<String>>()
                    .join("\n")
            )
        })
    }

    pub(super) fn handle_ordering(
        &mut self,
        best: &mut AlnBuffer,
        ord: Option<Ordering>,
    ) -> Result<Option<Decision>> {
        match ord {
            Some(Ordering::Greater) => self
                .handle_greater_than(best)
                .map(|_| self.add_decision_tag.then_some(Decision::First)),
            Some(Ordering::Less) => self
                .handle_less_than(best)
                .map(|_| self.add_decision_tag.then_some(Decision::Last)),
            Some(Ordering::Equal) => Ok(self.add_decision_tag.then_some(Decision::Ambiguous)),
            None => {
                // None of the alignments were fully unmapped or perfect matches,
                // so we need to score them to find the best.
                let first = &best.first().unwrap();
                let first_score = self.score_candidate(first, first.get_nr())?;

                let last = &best.last().unwrap();
                let last_score = self.score_candidate(last, last.get_nr())?;

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
    pub(super) fn handle_best(
        &mut self,
        best: &mut AlnBuffer,
        decision: Option<Decision>,
    ) -> Result<()> {
        let best_state = (best.len() == 1).then_some(true);

        best.drain(..).try_for_each(|mut b| {
            let nr = b.get_nr();
            b.records.drain(..).try_for_each(|mut r| {
                match decision {
                    Some(Decision::ConfDelta(decision)) => {
                        self.add_aux_tags(&mut r, b"XF", Aux::U8(decision))?
                    }
                    Some(Decision::VariantRescued(decision)) => {
                        self.add_aux_tags(&mut r, b"XR", Aux::U8(decision))?
                    }
                    _ => { /* no tag to add */ }
                }
                self.write_record(nr, r, best_state)
            })
        })
    }
    pub(super) fn handle_record_is_fragment_finished(
        &mut self,
        i: usize,
        rec: Record,
        best: &mut AlnBuffer,
    ) -> Result<bool> {
        if !(self.is_secondary_skipped)(&rec) {
            if let Some(new_readname) = (self.is_new_qname)(best, rec.qname()) {
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
                        state.records.push(rec);
                        return Ok(false);
                    }
                }
            }
            best.push(FragmentState::from_record(rec, i));
        } // else skip secondary

        Ok(false)
    }
}
