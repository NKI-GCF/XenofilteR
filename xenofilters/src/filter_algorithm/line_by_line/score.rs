use anyhow::{Result, anyhow};
use crate::alignment::{Fragment, stringify_record, SimpleRec,FragmentState, MdCigFlags};
use crate::variant::FragEvalVec;
use smallvec::SmallVec;
use super::core::LineByLine;

impl<R: SimpleRec> LineByLine<R> {
    pub(super) fn score_candidate<'v>(
        &'v mut self,
        state: &FragmentState<R>,
        aln_idx: usize,
    ) -> Result<f64> {
        let mut segment: SmallVec<[&R; 8]> = SmallVec::new();
        let mut md_cig_flags = SmallVec::with_capacity(state.get_records().len());
        let aln = self.aln.get(aln_idx).ok_or_else(|| anyhow!("No alignment for index {aln_idx}"))?;
        let mut dvnt: FragEvalVec<'_> = SmallVec::new();

        for idx in state.order_mates() {
            let rec = &state.get_records()[idx];
let flags = state.flags(idx).ok_or_else(|| anyhow!("No flags for record index {idx} in alignment {aln_idx}"))?;
            // TODO: add test to confirm secondary is always after non-secondary after order_mates.
            if flags.is_secondary() {
                if flags.is_last_segment() {
                    break;
                }
                continue;
            }
            if flags.is_unmapped() {
                dvnt.push(SmallVec::new());
            } else {
                let tid = rec.ref_seq_id().transpose()?
                    .ok_or_else(|| anyhow!("Mapped record has no reference sequence ID"))?;
                let start = rec.alignment_start().transpose()?
                    .ok_or_else(|| anyhow!("Mapped record has no alignment start"))?
                    .get();
                let end = start + rec.cigar().len();
                let delta_vars = if let Some(store) = aln.variant_store() {
                    store.overlapping_multi(tid, start, end)
                } else {
                    SmallVec::new()
                };
                dvnt.push(delta_vars);
            }
            segment.push(rec);
            md_cig_flags.push(MdCigFlags::try_from_record(rec, flags)?);
        }
        Fragment::new(&self.penalties, segment, md_cig_flags)?.score(&mut self.scratch, &mut dvnt).map_err(move |e| {
            anyhow!(
                "Error scoring fragment for alignment {aln_idx}: {}\n{}",
                e,
                state
                    .get_records()
                    .iter()
                    .map(stringify_record)
                    .collect::<Vec<String>>()
                    .join("\n")
            )
        })
    }
}
