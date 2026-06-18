use crate::alignment::{MdCigFlags, VariantWindow, weighted_ref_score, align_alt_to_read};
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
        // Variants fully covered mid-scan get moved here (see
        // score_variants_in_window) so they aren't re-processed/double-counted
        // by later windows in the same segment, but still reach maximize_delta.
        let mut finished: FragEvalVec<'v> = (0..dvnt.len()).map(|_| SmallVec::new()).collect();

        let mut score = self.score_with_variant(scratch, dvnt, &mut finished, 0)?;

        for i in 1..self.seg.len() {
            self.seg_i = i;
            self.refpos = self.seg_start[i];
            self.nt_i = 0;
            score += self.score_with_variant(scratch, dvnt, &mut finished, i)?;
        }

        for (d, f) in dvnt.iter_mut().zip(finished.iter_mut()) {
            d.extend(f.drain(..));
        }

        Ok(score + maximize_delta(dvnt, &mut scratch.dp))
    }

    fn score_with_variant<'v>(
        &mut self,
        scratch: &mut Scratch,
        dvnt: &mut FragEvalVec<'v>,
        finished: &mut FragEvalVec<'v>,
        i: usize,
    ) -> Result<f64> {
        let mut score = 0.0;

        for prior_seg_i in 0..self.seg_i {
            if self.md_cig_flags[prior_seg_i].is_last_segment()
                != self.md_cig_flags[self.seg_i].is_last_segment()
            {
                continue;
            }
            let seg_ref_start = self.seg_start[prior_seg_i];
            let seg_ref_end = seg_ref_start + self.seg[prior_seg_i].cigar().len();
            self.score_variants_in_window(
                scratch, dvnt, finished, i, prior_seg_i, seg_ref_start, seg_ref_end, 0.0,
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
                scratch, dvnt, finished, i, self.seg_i, ref_start, self.refpos, ref_score,
            )?;
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
            self.score_variants_in_window(
                scratch, dvnt, finished, i, next_seg_i, seg_ref_start, seg_ref_end, 0.0,
            )?;
        }

        Ok(score)
    }

    fn score_variants_in_window<'v>(
        &self,
        scratch: &mut Scratch,
        dvnt: &mut FragEvalVec<'v>,
        finished: &mut FragEvalVec<'v>,
        dvnt_i: usize,
        seg_i: usize,
        start: usize,
        end: usize,
        ref_score: f64,
    ) -> Result<()> {
        let mut i = 0;
        while i < dvnt[dvnt_i].len() && dvnt[dvnt_i][i].start() < end {
            if let Some((weighted_ref_score, alt_score)) =
                self.score_variant_in_seg(scratch, dvnt, dvnt_i, i, seg_i, start, end)?
            {
                dvnt[dvnt_i][i].update(weighted_ref_score, alt_score);
                let fully_processed =
                    dvnt[dvnt_i][i].ref_end() <= end && dvnt[dvnt_i][i].alt_end() <= end;
                if fully_processed {
                    let done = dvnt[dvnt_i].remove(i);
                    finished[dvnt_i].push(done);
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
        let vnt_eval = &dvnt[dvnt_i][dvnt_j];
        let window = match VariantWindow::compute(ref_start, ref_end, vnt_eval.start(), vnt_eval.ref_end()) {
            Some(w) => w,
            None => return Ok(None),
        };

        let vnt = vnt_eval.vnt();
        let alt = vnt.alt_allele();
        let p_variant = vnt.p_variant();
        let seg_ref_start = self.seg_start[seg_i];

        let weighted_ref = weighted_ref_score(window, seg_ref_start, p_variant, self.pen, |nt_i| {
            self.q(seg_i, nt_i)
        })?;

        let read_offset = window.read_offset(seg_ref_start);
        let revcmp = self.requires_revcmp(seg_i);
        let seq = self.seg[seg_i].sequence();
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
                Ok((read_base, self.q(seg_i, q_nt_i)?))
            },
            scratch,
        )?;

        #[cfg(test)]
        eprintln!(
            "[score_variant_in_seg] vnt=[{},{}) window={window:?} read_offset={read_offset} weighted_ref={weighted_ref} alt_score={alt_score}",
            vnt_eval.start(), vnt_eval.ref_end()
        );

        Ok(Some((weighted_ref, alt_score)))
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
