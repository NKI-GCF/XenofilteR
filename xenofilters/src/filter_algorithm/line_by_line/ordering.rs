use super::core::{AlnBuffer, LineByLine};
use crate::alignment::{Fragment, stringify_record};
use crate::fragment_state::FragmentState;
use anyhow::{Result, anyhow};
use noodles::bam::record::Record;
use std::cmp::Ordering;
use smallvec::SmallVec;
use noodles::sam::alignment::RecordBuf;

pub(crate) enum Decision {
    First,
    Last,
    Ambiguous,
    ConfDelta(u8),
    VariantRescued(u8),
}

impl LineByLine {
    fn handle_greater_than(&mut self, best: &mut AlnBuffer) -> Result<()> {
        let mut last = best.pop().unwrap();
        let nr = last.get_nr();
        last.records
            .drain(..)
            .try_for_each(|r| self.write_record(nr, &r, Some(false)))
    }
    fn handle_less_than(&mut self, best: &mut AlnBuffer) -> Result<()> {
        let all_before_last = best.len() - 1;
        best.drain(0..all_before_last).try_for_each(|mut b| {
            let nr = b.get_nr();
            b.records
                .drain(..)
                .try_for_each(|r| self.write_record(nr, &r, Some(false)))
        })
    }
    fn score_candidate(
        &self,
        state: &FragmentState<Record>,
        aln_idx: usize,
    ) -> Result<f64> {
        let mut dvnt_per_rec = SmallVec::with_capacity(state.records.len());
        let mut segment = SmallVec::new();
        let mut md_cig_flags = SmallVec::with_capacity(state.records.len());
        let aln = self.aln.get(aln_idx).ok_or_else(|| anyhow!("No alignment for index {aln_idx}"))?;

        for idx in state.order_mates(&self.aln) {
            let rec = &state.records[idx];
            let ops = &state.ops[idx];
            if ops.is_unmapped() {
                dvnt_per_rec.push(SmallVec::new());
            } else {
                let tid = rec.reference_sequence_id().transpose()?
                    .ok_or_else(|| anyhow!("Mapped record has no reference sequence ID"))?;
                let start = rec.alignment_start().transpose()?
                    .ok_or_else(|| anyhow!("Mapped record has no alignment start"))?
                    .get();
                let end = start + rec.cigar().len();
                let delta_vars = aln.variant_store()
                    .map(|s| s.overlapping_multi(tid, start, end))
                    .unwrap_or_default();
                dvnt_per_rec.push(delta_vars);
            }
            if !ops.is_secondary() {
                segment.push(rec);
                md_cig_flags.push(&state.ops[idx]);
            } else if ops.is_last_segment() {
                break;
            }
        }
        Fragment::new(&self.penalties, segment, md_cig_flags)?.score(dvnt_per_rec).map_err(|e| {
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
                .map(|()| self.add_decision_tag.then_some(Decision::First)),
            Some(Ordering::Less) => self
                .handle_less_than(best)
                .map(|()| self.add_decision_tag.then_some(Decision::Last)),
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
            b.records.drain(..).try_for_each(|r| {
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
                self.write_record(nr, &record_buf, best_state)
            })
        })
    }
    pub(super) fn handle_record_is_fragment_finished(
        &mut self,
        i: usize,
        rec: Record,
        best: &mut AlnBuffer,
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
}
