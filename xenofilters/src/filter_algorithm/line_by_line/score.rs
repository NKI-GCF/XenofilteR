use super::core::LineByLine;
use crate::alignment::{stringify_record, Fragment, FragmentState, MdCigFlags, SimpleRec};
use crate::filter_algorithm::line_by_line::READ_CT;
use crate::variant::FragEvalVec;
use anyhow::{anyhow, Result};
use smallvec::SmallVec;

impl<R: SimpleRec> LineByLine<R> {
    pub(super) fn score_candidate<'r>(
        &mut self,
        state: &'r FragmentState<R>,
        mcfs: SmallVec<[MdCigFlags<'r>; READ_CT]>,
        aln_idx: usize,
    ) -> Result<f64> {
        let mut segment: SmallVec<[&R; READ_CT]> = SmallVec::new();
        let mut seg_mcfs: SmallVec<[MdCigFlags; READ_CT]> = SmallVec::new();
        let mut mcfs_opt: SmallVec<[Option<MdCigFlags>; READ_CT]> =
            mcfs.into_iter().map(Some).collect();

        let aln = self
            .aln
            .get(aln_idx)
            .ok_or_else(|| anyhow!("No alignment for index {aln_idx}"))?;
        let store = aln.variant_store();
        let has_variants = store.is_some();

        let mut dvnt: FragEvalVec<'_> = SmallVec::new();
        let mut supplementary_penalty = 0.0;

        for idx in state.order_mates() {
            let rec = &state.get_records()[idx];
            let flags = state
                .flags(idx)
                .ok_or_else(|| anyhow!("No flags for record index {idx} in alignment {aln_idx}"))?;

            // Supplementary alignments contribute BOTH:
            //   1. A chimeric-junction penalty:
            //      gap_open + chimeric_junction_bases × gap_extend
            //   2. Per-base NW scoring via the segment below
            //      (so their actual alignment quality is still reflected in the score)
            if flags.is_supplementary() {
                supplementary_penalty += self.penalties.chimeric_junction_penalty;
            }

            if flags.is_unmapped() || !has_variants {
                dvnt.push(SmallVec::new());
            } else {
                let tid = rec
                    .ref_seq_id()
                    .transpose()?
                    .ok_or_else(|| anyhow!("Mapped record has no reference sequence ID"))?;
                let start = rec
                    .alignment_start()
                    .transpose()?
                    .ok_or_else(|| anyhow!("Mapped record has no alignment start"))?
                    .get();
                let cig_len = mcfs_opt[idx]
                    .as_ref()
                    .ok_or_else(|| anyhow!("MdCigFlags missing for record index {idx}"))?
                    .get_cigar()
                    .len();
                let end = start + cig_len;
                dvnt.push(store.as_ref().unwrap().overlapping_multi(tid, start, end));
            }

            if !flags.is_secondary() {
                segment.push(rec);
                seg_mcfs.push(
                    mcfs_opt[idx]
                        .take()
                        .ok_or_else(|| anyhow!("MdCigFlags already consumed for index {idx}"))?,
                );
            } else if flags.is_last_segment() {
                break;
            }
        }

        if segment.is_empty() {
            // Fragment consists only of supplementary records (unusual).
            return Ok(supplementary_penalty);
        }

        let base_score = Fragment::new(&self.penalties, segment, seg_mcfs)?
            .score(&mut self.scratch, &mut dvnt)
            .map_err(|e| {
                anyhow!(
                    "Error scoring fragment for alignment {aln_idx}: {}\n{}",
                    e,
                    state
                        .get_records()
                        .iter()
                        .map(stringify_record)
                        .collect::<Vec<_>>()
                        .join("\n")
                )
            })?;
        Ok(base_score + supplementary_penalty)
    }
}
