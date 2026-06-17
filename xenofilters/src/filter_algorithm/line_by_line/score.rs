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
        // mcfs is in raw record order; order_mates()/secondary-filtering below
        // reorders, so pull entries out by index rather than assuming alignment.
        let mut mcfs_opt: SmallVec<[Option<MdCigFlags>; READ_CT]> =
            mcfs.into_iter().map(Some).collect();

        let aln = self
            .aln
            .get(aln_idx)
            .ok_or_else(|| anyhow!("No alignment for index {aln_idx}"))?;
        let mut dvnt: FragEvalVec<'_> = SmallVec::new();

        for idx in state.order_mates() {
            let rec = &state.get_records()[idx];
            let flags = state
                .flags(idx)
                .ok_or_else(|| anyhow!("No flags for record index {idx} in alignment {aln_idx}"))?;
            if flags.is_unmapped() {
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
                let delta_vars = if let Some(store) = aln.variant_store() {
                    store.overlapping_multi(tid, start, end)
                } else {
                    SmallVec::new()
                };
                dvnt.push(delta_vars);
            }
            if !flags.is_secondary() {
                segment.push(rec);
                seg_mcfs.push(mcfs_opt[idx].take().ok_or_else(|| {
                    anyhow!("MdCigFlags already consumed for record index {idx}")
                })?);
            } else if flags.is_last_segment() {
                break;
            }
        }

        Fragment::new(&self.penalties, segment, seg_mcfs)?
            .score(&mut self.scratch, &mut dvnt)
            .map_err(move |e| {
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
