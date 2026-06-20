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
        let mut segment:  SmallVec<[&R; READ_CT]>              = SmallVec::new();
        let mut seg_mcfs: SmallVec<[MdCigFlags; READ_CT]>      = SmallVec::new();
        let mut mcfs_opt: SmallVec<[Option<MdCigFlags>; READ_CT]> =
            mcfs.into_iter().map(Some).collect();

        let aln = self.aln.get(aln_idx)
            .ok_or_else(|| anyhow!("No alignment for index {aln_idx}"))?;

        // Clone the Arc<dyn StoreTrait> (atomic increment, O(1)).
        // `None` when no variant file was supplied for this stream.
        let store = aln.variant_store();
        let has_variants = store.is_some();

        let mut dvnt: FragEvalVec<'_> = SmallVec::new();

        for idx in state.order_mates() {
            let rec   = &state.get_records()[idx];
            let flags = state.flags(idx)
                .ok_or_else(|| anyhow!("No flags for record index {idx} in alignment {aln_idx}"))?;

            if flags.is_unmapped() || !has_variants {
                // No variant lookup needed: push an empty eval vec.
                dvnt.push(SmallVec::new());
            } else {
                let tid = rec.ref_seq_id()
                    .transpose()?
                    .ok_or_else(|| anyhow!("Mapped record has no reference sequence ID"))?;
                let start = rec.alignment_start()
                    .transpose()?
                    .ok_or_else(|| anyhow!("Mapped record has no alignment start"))?
                    .get();
                let cig_len = mcfs_opt[idx]
                    .as_ref()
                    .ok_or_else(|| anyhow!("MdCigFlags missing for record index {idx}"))?
                    .get_cigar()
                    .len();
                let end = start + cig_len;
                // store is Some here (checked above via has_variants).
                dvnt.push(store.as_ref().unwrap().overlapping_multi(tid, start, end));
            }

            if !flags.is_secondary() {
                segment.push(rec);
                seg_mcfs.push(mcfs_opt[idx].take()
                    .ok_or_else(|| anyhow!("MdCigFlags already consumed for index {idx}"))?);
            } else if flags.is_last_segment() {
                break;
            }
        }

        Fragment::new(&self.penalties, segment, seg_mcfs)?
            .score(&mut self.scratch, &mut dvnt)
            .map_err(|e| {
                anyhow!(
                    "Error scoring fragment for alignment {aln_idx}: {}\n{}",
                    e,
                    state.get_records().iter()
                        .map(stringify_record)
                        .collect::<Vec<String>>()
                        .join("\n")
                )
            })
    }
}
