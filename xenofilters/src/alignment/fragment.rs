mod context;
mod record;

use crate::alignment::{
    align_alt_to_read, weighted_ref_score, BaseOp, MdCigFlags, ScoreOpIter, VariantWindow,
};
use crate::filter_algorithm::line_by_line::{Scratch, READ_CT};
use crate::penalty::{Penalty, MAX_Q};
use crate::variant::{Eval, FragEvalVec, VNT_CT};
use crate::Error;
use noodles::sam::alignment::Record;
use smallvec::SmallVec;

use context::{VariantCtx, WindowCtx};
pub(crate) use record::SimpleRec;

pub(crate) struct Fragment<'r, R> {
    pen: &'r Penalty,
    seg: SmallVec<[&'r R; READ_CT]>,
    md_cig_flags: SmallVec<[MdCigFlags<'r>; READ_CT]>,
    seg_start: SmallVec<[usize; READ_CT]>,
    seg_i: usize,
    refpos: usize,
    nt_i: usize,
}

/// Weighted Interval Scheduling over variant `Eval` entries with positive delta.
///
/// Sorts by `Eval::end()` then runs the classic O(n log n) WIS recurrence:
///   dp[i] = max(dp[i-1], delta_i + dp[latest non-overlapping predecessor])
///
/// Returns the maximum total rescue delta achievable with a non-overlapping
/// subset of variants.  Only entries with `delta() > 0.0` are considered;
/// entries with zero or negative delta cannot contribute to rescue.
///
/// `p_variant > 0.5` is a structural precondition for any entry to have
/// a positive delta; callers need not filter further.
pub fn wis_max_rescue_delta<'v>(
    dvnt: &mut FragEvalVec<'v>,
    dp: &mut SmallVec<[f64; READ_CT]>,
) -> f64 {
    let mut variants: SmallVec<[&Eval; VNT_CT]> =
        dvnt.iter().flatten().filter(|v| v.delta() > 0.0).collect();

    if variants.is_empty() {
        return 0.0;
    }

    // Sort by end coordinate for Weighted Interval Scheduling
    variants.sort_unstable_by_key(|v| v.end());

    dp.clear();
    dp.resize(variants.len(), 0.0);
    dp[0] = variants[0].delta();

    for i in 1..variants.len() {
        let current_weight = variants[i].delta();
        let current_start = variants[i].start();

        let mut incl_weight = current_weight;

        let latest = variants[..i].partition_point(|v| v.end() <= current_start);
        if latest > 0 {
            incl_weight += dp[latest - 1];
        }

        dp[i] = dp[i - 1].max(incl_weight);
    }
    *dp.last().unwrap_or(&0.0)
}

impl<'r, R: SimpleRec> Fragment<'r, R> {
    pub(crate) fn new(
        pen: &'r Penalty,
        seg: SmallVec<[&'r R; READ_CT]>,
        md_cig_flags: SmallVec<[MdCigFlags<'r>; READ_CT]>,
    ) -> Result<Self, Error> {
        let seg_start: SmallVec<[usize; READ_CT]> = seg
            .iter()
            .map(|r| {
                r.alignment_start()
                    .transpose()
                    .map(|o| o.map(|p| p.get().saturating_sub(1)).unwrap_or(0))
            })
            .collect::<Result<_, _>>()?;
        let refpos = seg_start[0];
        Ok(Self {
            pen,
            seg,
            md_cig_flags,
            seg_start,
            seg_i: 0,
            refpos,
            nt_i: 0,
        })
    }
    pub(crate) fn score<'v>(
        &mut self,
        scratch: &mut Scratch,
        dvnt: &mut FragEvalVec<'v>,
    ) -> Result<f64, Error> {
        // Ensure per-score mutable scanning state is reset so a Fragment
        // instance can be reused for multiple score() calls without leaking
        // seg_i/refpos/nt_i from previous runs (exposed by unit tests).
        self.seg_i = 0;
        self.refpos = self.seg_start[0];
        self.nt_i = 0;

        // Variants fully covered mid-scan get moved here (see
        // evaluate_variants_in_window) so they aren't re-processed/double-counted
        // by later windows in the same segment, but still reach wis_max_rescue_delta.
        let mut finished: FragEvalVec<'v> = (0..dvnt.len()).map(|_| SmallVec::new()).collect();

        let mut score = self.score_segment(scratch, dvnt, &mut finished, 0)?;

        for i in 1..self.seg.len() {
            self.seg_i = i;
            self.refpos = self.seg_start[i];
            self.nt_i = 0;
            score += self.score_segment(scratch, dvnt, &mut finished, i)?;
        }

        for (d, f) in dvnt.iter_mut().zip(finished.iter_mut()) {
            d.extend(f.drain(..));
        }

        let variant_delta = wis_max_rescue_delta(dvnt, &mut scratch.dp);
        scratch.last_variant_delta = variant_delta;
        Ok(score + variant_delta)
    }

    /// Score all aligned bases in segment `seg_i` and evaluate any overlapping
    /// variants from `dvnt[dvnt_i]`.
    ///
    /// Per-base scores accumulate into the return value.  Variants whose
    /// reference span falls fully within a scored window are moved to `finished`
    /// to prevent double-counting across adjacent windows; `wis_max_rescue_delta`
    /// merges them back before weighted-interval scheduling.
    fn score_segment<'v>(
        &mut self,
        scratch: &mut Scratch,
        dvnt: &mut FragEvalVec<'v>,
        finished: &mut FragEvalVec<'v>,
        i: usize,
    ) -> Result<f64, Error> {
        let mut score = 0.0;

        for prior_seg_i in 0..self.seg_i {
            if self.md_cig_flags[prior_seg_i].is_last_segment()
                != self.md_cig_flags[self.seg_i].is_last_segment()
            {
                continue;
            }
            let seg_ref_start = self.seg_start[prior_seg_i];
            let seg_ref_end = seg_ref_start + self.seg[prior_seg_i].cigar().len();
            let ctx = WindowCtx::new(i, prior_seg_i, seg_ref_start, seg_ref_end, 0.0);
            self.evaluate_variants_in_window(scratch, dvnt, finished, ctx)?;
        }

        let mut iter = ScoreOpIter::new(&self.md_cig_flags[self.seg_i]).peekable();

        while let Some(op_res) = iter.next() {
            let op = op_res?;
            let ref_start = self.refpos;
            let mut ref_score = 0.0;
            match op {
                BaseOp::Match => {
                    ref_score += self.pen.log_likelihood_match[self.q(self.seg_i, self.nt_i)?];
                    self.nt_i += 1;
                    self.refpos += 1;
                    while matches!(iter.peek(), Some(Ok(BaseOp::Match))) {
                        iter.next();
                        ref_score +=
                            self.pen.log_likelihood_match[self.q(self.seg_i, self.nt_i)?];
                        self.nt_i += 1;
                        self.refpos += 1;
                    }
                }
                BaseOp::Clip(len) => {
                    for _ in 0..len {
                        ref_score +=
                            self.pen.log_likelihood_mismatch[self.q(self.seg_i, self.nt_i)?];
                        self.nt_i += 1;
                    }
                    // refpos does NOT advance - soft clips consume read bases, not reference bases
                }
                BaseOp::Mis => {
                    ref_score += self.pen.log_likelihood_mismatch[self.q(self.seg_i, self.nt_i)?];
                    self.nt_i += 1;
                    self.refpos += 1;
                    while matches!(iter.peek(), Some(Ok(BaseOp::Mis))) {
                        iter.next();
                        ref_score +=
                            self.pen.log_likelihood_mismatch[self.q(self.seg_i, self.nt_i)?];
                        self.nt_i += 1;
                        self.refpos += 1;
                    }
                }
                BaseOp::Del(len) => {
                    self.refpos += len;
                    ref_score += self.pen.gap_open + (len as f64) * self.pen.gap_extend;
                }
                BaseOp::Ins(len) => {
                    self.nt_i += len;
                    ref_score += self.pen.gap_open + (len as f64) * self.pen.gap_extend;
                }
                BaseOp::RefSkip(len) => {
                    self.refpos += len;
                }
            }
            let ctx = WindowCtx::new(i, self.seg_i, ref_start, self.refpos, ref_score);
            self.evaluate_variants_in_window(scratch, dvnt, finished, ctx)?;
            score += ref_score;
        }

        for next_seg_i in (self.seg_i + 1)..self.seg.len() {
            if self.md_cig_flags[next_seg_i].is_last_segment()
                != self.md_cig_flags[self.seg_i].is_last_segment()
            {
                break;
            }
            let seg_ref_start = self.seg_start[next_seg_i];
            let seg_ref_end = seg_ref_start + self.seg[next_seg_i].sequence().len();
            let ctx = WindowCtx::new(i, next_seg_i, seg_ref_start, seg_ref_end, 0.0);
            self.evaluate_variants_in_window(scratch, dvnt, finished, ctx)?;
        }

        Ok(score)
    }

    fn evaluate_variants_in_window<'v>(
        &self,
        scratch: &mut Scratch,
        dvnt: &mut FragEvalVec<'v>,
        finished: &mut FragEvalVec<'v>,
        ctx: WindowCtx,
    ) -> Result<(), Error> {
        let mut i = 0;
        while i < dvnt[ctx.dvnt_i].len() && dvnt[ctx.dvnt_i][i].start() < ctx.ref_end {
            if let Some((weighted_ref_score, alt_score)) =
                self.score_variant_against_segment(scratch, dvnt, ctx.to_variant_ctx(i))?
            {
                dvnt[ctx.dvnt_i][i].update(weighted_ref_score, alt_score);
                // only the ref span needs to be within the window.
                let fully_processed = dvnt[ctx.dvnt_i][i].ref_end() <= ctx.ref_end;
                if fully_processed {
                    let done = dvnt[ctx.dvnt_i].remove(i);
                    finished[ctx.dvnt_i].push(done);
                    continue;
                }
            } else {
                // If a full scoring was not possible via score_variant_against_segment
                // (score_variant_against_segment returned None), add the expected
                // p-weighted reference contribution for the overlap portion (if any).
                // Do not use raw ctx.ref_score here because that is an unweighted
                // per-window score and would mix semantics.
                let v_eval = &dvnt[ctx.dvnt_i][i];
                let vnt = v_eval.vnt();

                // Keep the existing deferral rule for multi-base deletions: if the
                // variant ref span is longer than 1 and the ref_end is not yet within
                // this window, skip (defer) so we don't score partial deletions here.
                if vnt.ref_allele().len() > 1 && v_eval.ref_end() > ctx.ref_end {
                    // Defer to a later window that covers the full deletion: no update.
                } else if let Some(window) = VariantWindow::compute(
                    ctx.ref_start,
                    ctx.ref_end,
                    v_eval.start(),
                    v_eval.ref_end(),
                ) {
                    let seg_ref_start = self.seg_start[ctx.seg_i];
                    let p_variant = vnt.p_variant();
                    let wref =
                        weighted_ref_score(window, seg_ref_start, p_variant, self.pen, |nt_i| {
                            self.q(ctx.seg_i, nt_i)
                        })?;
                    dvnt[ctx.dvnt_i][i].update(wref, 0.0);
                }
            }
            i += 1;
        }
        Ok(())
    }

    fn requires_revcmp(&self, seg_i: usize) -> bool {
        if seg_i == self.seg_i {
            false
        } else {
            let current_seg_ori = self.md_cig_flags[self.seg_i].is_reverse_complemented();
            let other_seg_ori = self.md_cig_flags[seg_i].is_reverse_complemented();
            current_seg_ori != other_seg_ori
        }
    }
    /// Score a variant's alt allele against read bases from a specific segment.
    /// `ref_start` and `ref_end` define the reference window to score within this segment.
    fn score_variant_against_segment<'v>(
        &self,
        scratch: &mut Scratch,
        dvnt: &FragEvalVec<'v>,
        ctx: VariantCtx,
    ) -> Result<Option<(f64, f64)>, Error> {
        let vnt_eval = &dvnt[ctx.dvnt_i][ctx.dvnt_j];
        let vnt = vnt_eval.vnt();

        // Multi-base deletions: defer scoring until the full ref span is in window.
        // Partial overlap would pass incorrect ref_len to align_alt_to_read.
        // SNVs (ref_len=1) and insertions (ref_len=1, alt_len>1) are always safe.
        if vnt.ref_allele().len() > 1 && vnt_eval.ref_end() > ctx.ref_end {
            return Ok(None); // defer to a later window that covers the full deletion
        }
        let window = match VariantWindow::compute(
            ctx.ref_start,
            ctx.ref_end,
            vnt_eval.start(),
            vnt_eval.ref_end(),
        ) {
            Some(w) => w,
            None => return Ok(None),
        };

        let alt = vnt.alt_allele();
        let p_variant = vnt.p_variant();
        let seg_ref_start = self.seg_start[ctx.seg_i];

        let weighted_ref =
            weighted_ref_score(window, seg_ref_start, p_variant, self.pen, |nt_i| {
                self.q(ctx.seg_i, nt_i)
            })?;

        let read_offset = window.read_offset(seg_ref_start);
        let revcmp = self.requires_revcmp(ctx.seg_i);
        let seq = self.seg[ctx.seg_i].sequence();
        let seq_len = seq.len();

        let alt_score = align_alt_to_read(
            alt,
            read_offset,
            window.ref_len,
            p_variant,
            self.pen,
            |fwd_nt_i| {
                let (read_base, q_nt_i) = if revcmp {
                    let ri = seq_len - 1 - fwd_nt_i;
                    (complement(seq.get(ri).unwrap_or(b'N')), ri)
                } else {
                    (seq.get(fwd_nt_i).unwrap_or(b'N'), fwd_nt_i)
                };
                Ok((read_base, self.q(ctx.seg_i, q_nt_i)?))
            },
            scratch,
        )?;

        #[cfg(test)]
        eprintln!(
            "[score_variant_against_segment] vnt=[{},{}) window={window:?} read_offset={read_offset} weighted_ref={weighted_ref} alt_score={alt_score}",
            vnt_eval.start(),
            vnt_eval.ref_end()
        );

        Ok(Some((weighted_ref, alt_score)))
    }
    #[cfg(test)]
    pub(crate) fn seg_len(&self) -> usize {
        self.seg.len()
    }
}

impl<'r, R> Fragment<'r, R>
where
    R: Record + SimpleRec,
{
    pub(crate) fn q(&self, seg_i: usize, nt_i: usize) -> Result<usize, Error> {
        self.seg[seg_i]
            .quality_at(nt_i)
            .map(|q| (q as usize).min(MAX_Q - 1))
            .ok_or(Error::QualityScoreOutOfBounds { nt_i, seg_i })
    }
}

fn complement(b: u8) -> u8 {
    match b {
        b'A' => b'T',
        b'C' => b'G',
        b'G' => b'C',
        b'T' => b'A',
        _ => b'N',
    }
}

#[cfg(test)]
pub(crate) mod tests;
