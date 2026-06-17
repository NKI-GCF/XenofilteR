use crate::alignment::MdCigFlags;
use crate::alignment::{AlignmentError, BaseOp, ScoreOpIter};
use crate::filter_algorithm::line_by_line::{Scratch, READ_CT};
use crate::penalty::{Penalty, MAX_Q};
use crate::variant::{Eval, FragEvalVec, VNT_CT};
use anyhow::{anyhow, Result};
use noodles::sam::alignment::Record;
use noodles::sam::alignment::RecordBuf;
use noodles::sam::Header;
use smallvec::SmallVec;

pub(crate) struct Fragment<'r, R> {
    pen: &'r Penalty,
    seg: SmallVec<[&'r R; READ_CT]>,
    md_cig_flags: SmallVec<[MdCigFlags<'r>; READ_CT]>,
    seg_start: SmallVec<[usize; READ_CT]>,
    seg_i: usize,
    refpos: usize,
    nt_i: usize,
}

pub(super) fn maximize_delta<'v>(
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
    ) -> Result<Self, AlignmentError> {
        let seg_start: SmallVec<[usize; READ_CT]> = seg
            .iter()
            .map(|r| {
                r.alignment_start()
                    .transpose()
                    .map(|o| o.map(|p| p.get()).unwrap_or(0))
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
    ) -> Result<f64, AlignmentError> {
        let mut score = self.score_with_variant(scratch, dvnt, 0)?;

        for i in 1..self.seg.len() {
            self.seg_i = i;
            self.refpos = self.seg_start[i];
            self.nt_i = 0;
            score += self.score_with_variant(scratch, dvnt, i)?;
        }
        Ok(score + maximize_delta(dvnt, &mut scratch.dp))
    }

    fn score_with_variant<'v>(
        &mut self,
        scratch: &mut Scratch,
        dvnt: &mut FragEvalVec<'v>,
        i: usize,
    ) -> Result<f64> {
        let mut score = 0.0;

        // Score variant bases that reach into segments before the current one.
        for prior_seg_i in 0..self.seg_i {
            if self.md_cig_flags[prior_seg_i].is_last_segment()
                != self.md_cig_flags[self.seg_i].is_last_segment()
            {
                continue; // skip read 1 read(s) when processing read 2
            }
            let seg_ref_start = self.seg_start[prior_seg_i];
            let seg_ref_end = seg_ref_start + self.seg[prior_seg_i].cigar().len();

            // prior (or following) segment bases contribute to alt scoring only, no incurrence.
            self.score_variants_in_window(
                scratch,
                dvnt,
                i,
                prior_seg_i,
                seg_ref_start,
                seg_ref_end,
                0.0,
            )?;
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
                BaseOp::Relocate { penalty_score, pos } => {
                    self.refpos = pos;
                    ref_score += penalty_score;
                }
                BaseOp::RefSkip(len) => {
                    self.refpos += len;
                }
            }

            self.score_variants_in_window(
                scratch,
                dvnt,
                i,
                self.seg_i,
                ref_start,
                self.refpos,
                ref_score,
            )?;
            score += ref_score;
        }

        for next_seg_i in (self.seg_i + 1)..self.seg.len() {
            if self.md_cig_flags[next_seg_i].is_last_segment()
                != self.md_cig_flags[self.seg_i].is_last_segment()
            {
                break; // skip segments from read 2 when processing read 1
            }
            let seg_ref_start = self.seg_start[next_seg_i];
            let seg_ref_end = seg_ref_start + self.seg[next_seg_i].sequence().len();

            self.score_variants_in_window(
                scratch,
                dvnt,
                i,
                next_seg_i,
                seg_ref_start,
                seg_ref_end,
                0.0,
            )?;
        }

        Ok(score)
    }
    fn score_variants_in_window<'v>(
        &self,
        scratch: &mut Scratch,
        dvnt: &mut FragEvalVec<'v>,
        dvnt_i: usize,
        seg_i: usize,
        start: usize,
        end: usize,
        ref_score: f64, // unweighted, for non-variant positions
    ) -> Result<()> {
        let mut i = 0;
        while i < dvnt[dvnt_i].len() && dvnt[dvnt_i][i].start() < end {
            if let Some((weighted_ref_score, alt_score)) =
                self.score_variant_in_seg(scratch, dvnt, dvnt_i, i, seg_i, start, end)?
            {
                // Use weighted_ref_score instead of ref_score for variant positions
                dvnt[dvnt_i][i].update(weighted_ref_score, alt_score);
                let fully_processed =
                    dvnt[dvnt_i][i].ref_end() <= end && dvnt[dvnt_i][i].alt_end() <= end;
                if fully_processed {
                    dvnt[dvnt_i].remove(i);
                    continue;
                }
            } else {
                dvnt[dvnt_i][i].update(ref_score, 0.0);
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
    fn score_variant_in_seg<'v>(
        &self,
        scratch: &mut Scratch,
        dvnt: &FragEvalVec<'v>,
        dvnt_i: usize,
        dvnt_j: usize,
        seg_i: usize,
        ref_start: usize,
        ref_end: usize,
    ) -> Result<Option<(f64, f64)>> {
        // (weighted_ref_score, alt_score)
        let vnt_eval = &dvnt[dvnt_i][dvnt_j];
        let vnt_start = vnt_eval.start();

        // Clamp the ref window to the variant's ref span
        let eff_ref_start = ref_start.max(vnt_start);
        let eff_ref_end = ref_end.min(vnt_eval.ref_end());
        let ref_len = eff_ref_end - eff_ref_start;

        // Derive the alt slice from the ref offset (valid for MNPs; see caveat below)
        let ref_consumed = eff_ref_start - vnt_eval.start();
        let vnt = vnt_eval.vnt();
        let alt = vnt.alt_allele();
        let alt_slice = &alt[ref_consumed.min(alt.len())..(ref_consumed + ref_len).min(alt.len())];
        if alt_slice.is_empty() && ref_len == 0 {
            return Ok(None);
        }

        let p_variant = vnt.p_variant();
        let seg_ref_start = self.seg_start[seg_i];
        let nt_offset = eff_ref_start - seg_ref_start;

        // Weighted ref score over the overlapping bases
        let mut weighted_ref_score = 0.0;
        for j in 0..(eff_ref_end - eff_ref_start) {
            let q = self.q(seg_i, nt_offset + j)?;
            let lm = self.pen.log_likelihood_match[q];
            let lmm = self.pen.log_likelihood_mismatch[q];
            // ref path: read matches ref, so weight by (1 - p_variant)
            weighted_ref_score += (1.0 - p_variant) * lm + p_variant * lmm;
        }

        let n = alt.len();
        let m = n.max(ref_len);

        scratch.resize_nw(m + 1);

        // Initialise first row: gaps in the read (deletions relative to alt)
        scratch.prev[0].reinit(self.pen.gap_open, self.pen.gap_extend, 0);
        for (j, p) in scratch.prev.iter_mut().enumerate().take(m + 1).skip(1) {
            p.reinit(self.pen.gap_open, self.pen.gap_extend, -(j as i32));
        }

        let nt_i_base = eff_ref_start - seg_ref_start;
        let revcmp = self.requires_revcmp(seg_i);
        let seq = self.seg[seg_i].sequence();
        let seq_len = seq.len();

        for i in 1..=n {
            let alt_base = alt[i - 1];
            scratch.curr[0].reinit(self.pen.gap_open, self.pen.gap_extend, i as i32);

            for j in 1..=m {
                let fwd_nt_i = nt_i_base + j - 1;
                let (read_base, q_nt_i) = if revcmp {
                    let ri = seq_len - 1 - fwd_nt_i;
                    let b = seq.get(ri).unwrap_or(b'N');
                    (complement(b), ri)
                } else {
                    (seq.get(fwd_nt_i).unwrap_or(b'N'), fwd_nt_i)
                };
                let q = self.q(seg_i, q_nt_i)?;

                let lm = self.pen.log_likelihood_match[q];
                let lmm = self.pen.log_likelihood_mismatch[q];

                // Weight match/mismatch score by variant probability,
                // mirroring score_alt_match / score_ref_match
                let per_base_score = if alt_base == read_base {
                    p_variant * lm + (1.0 - p_variant) * lmm
                } else {
                    (1.0 - p_variant) * lm + p_variant * lmm
                };

                scratch.curr[j].m = per_base_score
                    + scratch.prev[j - 1]
                        .m
                        .max(scratch.prev[j - 1].i)
                        .max(scratch.prev[j - 1].d);

                scratch.curr[j].i =
                    (scratch.curr[j - 1].m + self.pen.gap_open + self.pen.gap_extend)
                        .max(scratch.curr[j - 1].i + self.pen.gap_extend);

                scratch.curr[j].d = (scratch.prev[j].m + self.pen.gap_open + self.pen.gap_extend)
                    .max(scratch.prev[j].d + self.pen.gap_extend);
            }

            scratch.swap_nw();
        }

        let alt_score = scratch.prev[m]
            .m
            .max(scratch.prev[m].i)
            .max(scratch.prev[m].d);

        Ok(Some((weighted_ref_score, alt_score)))
    }
}

pub(crate) trait SimpleRec: Record + PartialEq {
    fn quality_at(&self, i: usize) -> Option<u8>;
    fn ref_seq_id(&self) -> Option<Result<usize, std::io::Error>>;
    fn as_record_buf(&self, header: &Header) -> Result<RecordBuf, std::io::Error>;
}

impl SimpleRec for noodles::bam::Record {
    fn quality_at(&self, i: usize) -> Option<u8> {
        self.quality_scores().as_ref().get(i).copied()
    }
    fn ref_seq_id(&self) -> Option<Result<usize, std::io::Error>> {
        self.reference_sequence_id()
    }
    fn as_record_buf(&self, header: &Header) -> Result<RecordBuf, std::io::Error> {
        RecordBuf::try_from_alignment_record(header, self)
    }
}

impl SimpleRec for noodles::sam::alignment::RecordBuf {
    fn quality_at(&self, i: usize) -> Option<u8> {
        self.quality_scores().as_ref().get(i).copied()
    }
    fn ref_seq_id(&self) -> Option<Result<usize, std::io::Error>> {
        self.reference_sequence_id().map(Ok)
    }
    fn as_record_buf(&self, _header: &Header) -> Result<RecordBuf, std::io::Error> {
        Ok(self.clone())
    }
}

impl<'r, R> Fragment<'r, R>
where
    R: Record + SimpleRec,
{
    fn q(&self, seg_i: usize, nt_i: usize) -> Result<usize> {
        self.seg[seg_i]
            .quality_at(nt_i)
            .map(|q| (q as usize).min(MAX_Q - 1))
            .ok_or_else(|| anyhow!("Quality score index {nt_i} out of bounds for segment {seg_i}"))
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
