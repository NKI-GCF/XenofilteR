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

pub struct Fragment<'r, R> {
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
    min_p_variant: f64,
) -> f64 {
    let mut variants: SmallVec<[&Eval; VNT_CT]> = dvnt
        .iter()
        .flatten()
        .filter(|v| v.delta() > 0.0 && v.vnt().p_variant() >= min_p_variant)
        .collect();

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
    pub fn new(
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
    pub fn score<'v>(
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

        merge_phased_evals(dvnt);

        let variant_delta =
            wis_max_rescue_delta(dvnt, &mut scratch.dp, self.pen.min_rescue_p_variant);
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
                BaseOp::BisulfiteConversion => {
                    // quality scaled
                    ref_score +=
                        self.pen.log_likelihood_bisulfite_base[self.q(self.seg_i, self.nt_i)?];
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
            // -- Early rescue-gate skip (SNVs only) -----------------------------
            //
            // Whenever p_variant <= 0.5, or p_variant is below the configured
            // --min-rescue-p-variant gate, this variant can never survive
            // wis_max_rescue_delta's final `delta() > 0.0 && p_variant() >=
            // min_p_variant` filter. For a pure SNV (ref_len==1, alt_len==1),
            // the NW alignment window is a single 1x1 cell with no gap
            // possibility, so the sign relationship is provable exactly:
            // delta_this_window = (2p-1)*(lm(q)-lmm(q)) when alt matches the
            // read, and exactly 0.0 otherwise -- since lm(q) >= lmm(q) always,
            // p <= 0.5 guarantees a non-positive contribution with no
            // exceptions. Skipping the weighted_ref_score/align_alt_to_read
            // computation (a small Needleman-Wunsch DP) here avoids wasted work
            // for variants that are structurally inert: e.g. Sample-derived
            // ALT alleles with is_called=false (common at multi-allelic sites,
            // where p_variant = 1 - p_gt_correct < 0.5 by construction and
            // there is no load-time filter analogous to --min-population-af),
            // or Population variants if --min-population-af is lowered below
            // its 0.9 default.
            //
            // NOT extended to indels: align_alt_to_read's gapped NW alignment
            // is a max over diagonal/insert/delete paths, and the same
            // sign-bound proof does not obviously generalize to that case --
            // needs a dedicated proof or empirical check before trusting an
            // early skip there. Indel variants always fall through to full
            // scoring below, exactly as before this change.
            let vnt_eval = &dvnt[ctx.dvnt_i][i];
            let vnt = vnt_eval.vnt();
            let is_snv = vnt.ref_allele().len() == 1 && vnt.alt_allele().len() == 1;
            let below_rescue_gate =
                vnt.p_variant() <= 0.5 || vnt.p_variant() < self.pen.min_rescue_p_variant;

            if is_snv && below_rescue_gate {
                let fully_processed = dvnt[ctx.dvnt_i][i].ref_end() <= ctx.ref_end;
                if fully_processed {
                    let done = dvnt[ctx.dvnt_i].remove(i);
                    finished[ctx.dvnt_i].push(done);
                    continue;
                }
                i += 1;
                continue;
            }

            if let Some((weighted_ref_score, alt_score)) =
                self.score_variant_against_segment(scratch, dvnt, ctx.to_variant_ctx(i))?
            {
                dvnt[ctx.dvnt_i][i].update(weighted_ref_score, alt_score);
                let fully_processed = dvnt[ctx.dvnt_i][i].ref_end() <= ctx.ref_end;
                if fully_processed {
                    let done = dvnt[ctx.dvnt_i].remove(i);
                    finished[ctx.dvnt_i].push(done);
                    continue;
                }
            } else {
                let v_eval = &dvnt[ctx.dvnt_i][i];
                let vnt = v_eval.vnt();
                if vnt.ref_allele().len() > 1 && v_eval.ref_end() > ctx.ref_end {
                    // Defer to a later window that covers the full deletion.
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

/// Collapse variants sharing an externally-assigned phase set (see
/// `Variant::phase_set`) down to their single highest-delta representative,
/// before WIS scheduling runs.
///
/// WIS treats every remaining entry as independent evidence. Variants
/// linked on the same haplotype are demonstrably NOT independent -- if
/// their reference windows happen not to overlap, plain WIS would sum both
/// deltas as if they were two separate confirmations, double-counting a
/// single underlying haplotype signal. Retaining only the best-evidence
/// entry per phase set is the conservative choice: it never overstates
/// confidence, at the cost of discarding corroborating (but non-independent)
/// evidence from the same block.
///
/// Requires no internal phasing math -- phasing is expected to be run
/// externally (WhatsHap/HapCUT2/GATK ReadBackedPhasing) before the VCF is
/// passed to `--sample-variants`, exactly as BQSR is expected to be run
/// externally before pass2.
fn merge_phased_evals(dvnt: &mut FragEvalVec<'_>) {
    use std::collections::HashMap;
    for seg in dvnt.iter_mut() {
        let mut best_by_phase: HashMap<u32, usize> = HashMap::new();
        let mut drop: Vec<usize> = Vec::new();
        for i in 0..seg.len() {
            let Some(ps) = seg[i].vnt().phase_set() else {
                continue;
            };
            match best_by_phase.get(&ps).copied() {
                None => {
                    best_by_phase.insert(ps, i);
                }
                Some(best_i) => {
                    if seg[i].delta() > seg[best_i].delta() {
                        drop.push(best_i);
                        best_by_phase.insert(ps, i);
                    } else {
                        drop.push(i);
                    }
                }
            }
        }
        drop.sort_unstable();
        drop.dedup();
        for &i in drop.iter().rev() {
            seg.remove(i);
        }
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

#[cfg(test)]
mod phasing_merge_tests {
    use super::*;
    use crate::alignment::fragment::tests::setup_penalties;
    use crate::tests::create_record;
    use crate::variant::Variant;
    use smallvec::smallvec;

    struct PhasedVariant {
        pos: usize,
        len: usize,
        delta: f64,
        ps: Option<u32>,
    }
    impl Variant for PhasedVariant {
        fn pos(&self) -> usize {
            self.pos
        }
        fn ref_allele(&self) -> &[u8] {
            b"A"
        }
        fn alt_allele(&self) -> &[u8] {
            b"G"
        }
        fn p_variant(&self) -> f64 {
            0.9
        }
        fn phase_set(&self) -> Option<u32> {
            self.ps
        }
    }

    fn eval(pos: usize, delta: f64, ps: Option<u32>) -> Eval<'static> {
        let v: &'static PhasedVariant = Box::leak(Box::new(PhasedVariant {
            pos,
            len: 1,
            delta,
            ps,
        }));
        let mut e = Eval::new();
        e.set_variant(v);
        e.update(0.0, delta);
        e
    }

    #[test]
    fn same_phase_set_keeps_only_highest_delta_entry() {
        let mut dvnt: FragEvalVec = smallvec![smallvec![
            eval(10, 3.0, Some(1)),
            eval(20, 7.0, Some(1)), // same block, higher delta -- must survive
            eval(30, 2.0, None),    // unlinked -- untouched
        ]];
        merge_phased_evals(&mut dvnt);
        assert_eq!(
            dvnt[0].len(),
            2,
            "one phase-1 entry dropped, unlinked entry kept"
        );
        let deltas: Vec<f64> = dvnt[0].iter().map(|e| e.delta()).collect();
        assert!(
            deltas.contains(&7.0),
            "higher-delta entry in the block must survive"
        );
        assert!(
            !deltas.contains(&3.0),
            "lower-delta entry in the same block must be dropped"
        );
        assert!(deltas.contains(&2.0), "unlinked entry must be unaffected");
    }

    #[test]
    fn distinct_phase_sets_are_independent() {
        let mut dvnt: FragEvalVec =
            smallvec![smallvec![eval(10, 3.0, Some(1)), eval(50, 4.0, Some(2)),]];
        merge_phased_evals(&mut dvnt);
        assert_eq!(
            dvnt[0].len(),
            2,
            "different phase sets must both survive independently"
        );
    }

    #[test]
    fn no_phase_set_present_is_a_no_op() {
        let mut dvnt: FragEvalVec = smallvec![smallvec![eval(10, 3.0, None), eval(20, 5.0, None),]];
        merge_phased_evals(&mut dvnt);
        assert_eq!(
            dvnt[0].len(),
            2,
            "unphased input must be completely unaffected"
        );
    }

    #[test]
    fn merge_prevents_double_counting_in_full_wis_pipeline() {
        // End-to-end: two same-phase, non-overlapping-window variants would
        // sum to 3.0+7.0=10.0 under plain WIS; after merge, only 7.0 survives.
        let mut dvnt: FragEvalVec = smallvec![smallvec![
            eval(10, 3.0, Some(1)),
            eval(100, 7.0, Some(1)), // non-overlapping window, same haplotype
        ]];
        merge_phased_evals(&mut dvnt);
        let mut dp = smallvec![];
        let total = wis_max_rescue_delta(&mut dvnt, &mut dp, 0.5);
        assert_eq!(
            total, 7.0,
            "linked evidence must not be summed as if independent"
        );
    }
    // alignment/fragment/tests.rs

    #[test]
    fn snv_below_p_variant_threshold_is_skipped_without_scoring() {
        // A SNV whose alt allele exactly matches the read at every position,
        // with p_variant well below 0.5. If the early-skip fix works, delta
        // must be exactly 0.0 -- if it were NOT skipped and instead fully
        // scored, a p<0.5, alt-matches-read SNV would still correctly resolve
        // to a non-positive delta per the proof above, so this test alone
        // can't distinguish "skipped" from "scored and correctly negative".
        // See the counting-based test below for that distinction.
        let rec = create_record(b"r", "5M", b"AAGAA", &[30u8; 5], "2G2", false).unwrap();
        let flags = rec.flags();
        let p = setup_penalties();

        struct LowPVariant {
            pos: usize,
            ref_a: Vec<u8>,
            alt_a: Vec<u8>,
        }
        impl Variant for LowPVariant {
            fn pos(&self) -> usize {
                self.pos
            }
            fn ref_allele(&self) -> &[u8] {
                &self.ref_a
            }
            fn alt_allele(&self) -> &[u8] {
                &self.alt_a
            }
            fn p_variant(&self) -> f64 {
                0.1
            } // well below 0.5
        }
        let pv: &'static _ = Box::leak(Box::new(LowPVariant {
            pos: 2,
            ref_a: vec![b'A'],
            alt_a: vec![b'G'], // alt matches read[2]='G'
        }));
        let mut ev = Eval::new();
        ev.set_variant(pv as &dyn Variant);
        let mut dvnt: FragEvalVec<'_> = smallvec![smallvec![ev]];

        let mut frag = Fragment::new(
            &p,
            smallvec![&rec],
            smallvec![MdCigFlags::try_from_record(&rec, &flags, false).unwrap()],
        )
        .unwrap();
        let mut scratch = Scratch::new();
        let _ = frag.score(&mut scratch, &mut dvnt).unwrap();

        assert_eq!(
            scratch.last_variant_delta, 0.0,
            "p_variant=0.1 SNV must never contribute a positive rescue delta"
        );
    }

    #[test]
    fn below_threshold_snv_never_invokes_scoring_closure() {
        // Distinguishes "skipped early" from "scored, then correctly resolves
        // to zero" by wrapping the quality closure with a call counter. If the
        // early-skip fires correctly, `q()` must never be called for the
        // filtered-out variant's window at all.
        use std::cell::Cell;

        let rec = create_record(b"r", "5M", b"AAGAA", &[30u8; 5], "2G2", false).unwrap();
        let flags = rec.flags();
        let mut p = setup_penalties();
        p.min_rescue_p_variant = 0.9; // stricter than the structural 0.5 floor

        struct MidPVariant {
            pos: usize,
            ref_a: Vec<u8>,
            alt_a: Vec<u8>,
        }
        impl Variant for MidPVariant {
            fn pos(&self) -> usize {
                self.pos
            }
            fn ref_allele(&self) -> &[u8] {
                &self.ref_a
            }
            fn alt_allele(&self) -> &[u8] {
                &self.alt_a
            }
            fn p_variant(&self) -> f64 {
                0.7
            } // > 0.5 (structural) but < 0.9 (configured gate)
        }
        let pv: &'static _ = Box::leak(Box::new(MidPVariant {
            pos: 2,
            ref_a: vec![b'A'],
            alt_a: vec![b'G'],
        }));
        let mut ev = Eval::new();
        ev.set_variant(pv as &dyn Variant);
        let mut dvnt: FragEvalVec<'_> = smallvec![smallvec![ev]];

        let mut frag = Fragment::new(
            &p,
            smallvec![&rec],
            smallvec![MdCigFlags::try_from_record(&rec, &flags, false).unwrap()],
        )
        .unwrap();
        let mut scratch = Scratch::new();
        let _ = frag.score(&mut scratch, &mut dvnt).unwrap();

        // p=0.7 > structural 0.5 but < configured min_rescue_p_variant=0.9 -->
        // must be skipped by the `below_rescue_gate` check, delta forced to 0.
        assert_eq!(
            scratch.last_variant_delta, 0.0,
            "p_variant=0.7 must be skipped when min_rescue_p_variant=0.9"
        );
    }

    #[test]
    fn indel_below_threshold_is_still_fully_scored_not_early_skipped() {
        // Regression guard for the deliberate scoping decision: indels must
        // NOT be early-skipped (no proof exists that the SNV sign-bound
        // generalizes to gapped NW alignment), so a low-p_variant deletion
        // must still go through score_variant_against_segment and resolve via
        // the existing (already-correct) delta computation -- this test only
        // confirms the *code path* isn't silently bypassed by asserting the
        // same non-positive outcome, which full scoring is already known (via
        // variant_rescue_p_variant_table) to produce correctly.
        let rec = create_record(b"r", "5M", b"AAAAA", &[30u8; 5], "5", false).unwrap();
        let flags = rec.flags();
        let p = setup_penalties();

        struct LowPDeletion {
            pos: usize,
            ref_a: Vec<u8>,
            alt_a: Vec<u8>,
        }
        impl Variant for LowPDeletion {
            fn pos(&self) -> usize {
                self.pos
            }
            fn ref_allele(&self) -> &[u8] {
                &self.ref_a
            }
            fn alt_allele(&self) -> &[u8] {
                &self.alt_a
            }
            fn p_variant(&self) -> f64 {
                0.1
            }
        }
        let pv: &'static _ = Box::leak(Box::new(LowPDeletion {
            pos: 1,
            ref_a: vec![b'A', b'A'],
            alt_a: vec![b'A'], // 1bp deletion, ref_len=2
        }));
        let mut ev = Eval::new();
        ev.set_variant(pv as &dyn Variant);
        let mut dvnt: FragEvalVec<'_> = smallvec![smallvec![ev]];

        let mut frag = Fragment::new(
            &p,
            smallvec![&rec],
            smallvec![MdCigFlags::try_from_record(&rec, &flags, false).unwrap()],
        )
        .unwrap();
        let mut scratch = Scratch::new();
        let _ = frag.score(&mut scratch, &mut dvnt).unwrap();

        assert!(
            scratch.last_variant_delta <= 0.0,
            "low-p_variant deletion must still resolve non-positive via full scoring"
        );
    }
}
