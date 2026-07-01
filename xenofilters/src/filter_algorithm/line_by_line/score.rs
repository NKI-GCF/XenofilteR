use super::core::LineByLine;
use crate::alignment::{
    mate_slot, segment_id, stringify_record, Fragment, FragmentState, MdCigFlags, SimpleRec,
};
use crate::filter_algorithm::line_by_line::READ_CT;
use crate::variant::FragEvalVec;
use crate::Error; // Assuming your error enum lives here
use smallvec::SmallVec;

impl<R: SimpleRec> LineByLine<R> {
    pub(super) fn score_candidate<'r>(
        &mut self,
        state: &'r FragmentState<R>,
        mcfs: SmallVec<[MdCigFlags<'r>; READ_CT]>,
        aln_idx: usize,
        cancel_slot: [bool; 2],
    ) -> Result<f64, Error> {
        let mut segment: SmallVec<[&R; READ_CT]> = SmallVec::new();
        let mut seg_mcfs: SmallVec<[MdCigFlags; READ_CT]> = SmallVec::new();
        let mut mcfs_opt: SmallVec<[Option<MdCigFlags>; READ_CT]> =
            mcfs.into_iter().map(Some).collect();

        let aln = self
            .aln
            .get(aln_idx)
            .ok_or(Error::NoAlignmentForIndex { aln_idx })?;
        let store = aln.variant_store();
        let has_variants = store.is_some();

        let mut dvnt: FragEvalVec<'_> = SmallVec::new();
        let mut supplementary_penalty = 0.0;

        for idx in state.order_mates() {
            let rec = &state.get_records()[idx];
            let flags = state
                .flags(idx)
                .ok_or(Error::NoFlagsForRecordInAlignment { idx, aln_idx })?;

            // Skip every record (primary or supplementary) belonging to a
            // unanimously-cancelled mate slot: its contribution is provably
            // identical across all competing streams, so it is excluded from
            // the NW segment entirely rather than scored and subtracted out.
            let slot = mate_slot(segment_id(flags));
            if cancel_slot[slot] {
                mcfs_opt[idx] = None; // release the borrow; never consumed
                continue;
            }

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
                    .ok_or(Error::MappedRecordNoReferenceSequenceId)?;
                let start = rec
                    .alignment_start()
                    .transpose()?
                    .ok_or(Error::MappedRecordNoAlignmentStart)?
                    .get();
                let cig_len = mcfs_opt[idx]
                    .as_ref()
                    .ok_or(Error::MdCigFlagsMissing { idx })?
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
                        .ok_or(Error::MdCigFlagsAlreadyConsumed { idx })?,
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
                let state_str = state
                    .get_records()
                    .iter()
                    .map(stringify_record)
                    .collect::<Vec<_>>()
                    .join("\n");
                Error::FragmentScoringError {
                    aln_idx,
                    message: e.to_string(),
                    state: state_str,
                }
            })?;
        Ok(base_score + supplementary_penalty)
    }
}
