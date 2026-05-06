use crate::{MAX_Q, Penalty};
use anyhow::{Result, anyhow};
use rust_htslib::bam::record::Record;
use smallvec::SmallVec;
use crate::variant::{DeltaPerRec, VariantEval, VntPerRec};
use crate::alignment::{AlignmentError, UnifiedOpIterator, UnifiedOp};

#[derive(Clone, Copy)]
struct Cell {
    m: f64, // Match/Mismatch
    i: f64, // Insertion (gap in Alt)
    d: f64, // Deletion (gap in Read)
}

impl Cell {
    fn reinit(&mut self, gap_open: f64, gap_extend: f64, i: i32) {
        self.m = f64::NEG_INFINITY;
        self.i = f64::NEG_INFINITY;
        self.d = f64::NEG_INFINITY;
        match i {
            0 => self.m = 0.0,
            i if i < 0 => self.i = gap_open + (i.abs() as f64) * gap_extend,
            i => self.d = gap_open + (i as f64) * gap_extend,
        }
    }
}

impl Default for Cell {
    fn default() -> Self {
        Self { m: -f64::INFINITY, i: -f64::INFINITY, d: -f64::INFINITY }
    }
}

pub(crate) struct Fragment<'r> {
    pen: &'r Penalty,
    seg: SmallVec<[&'r Record; 2]>,
    seg_i: usize,
    refpos: i64,
    nt_i: usize,
    prev: Vec<Cell>,
    curr: Vec<Cell>,
}

impl<'r> Fragment<'r> {
    pub(crate) fn new(pen: &'r Penalty, seg: SmallVec<[&'r Record; 2]>) -> Self {
        let refpos = seg[0].pos();
        Self {
            pen,
            seg,
            seg_i: 0,
            refpos,
            nt_i: 0,
            prev: Vec::new(),
            curr: Vec::new(),
        }
    }
    pub(crate) fn score<'v>(&mut self, mut vnt_per_rec: VntPerRec<'v>) -> Result<f64, AlignmentError>
    {
        let mut score = self.score_with_variant(&mut vnt_per_rec[0])?;

        for i in 1..self.seg.len() {
            self.seg_i = i;
            self.refpos = self.seg[i].pos();
            self.nt_i = 0;
            score += self.score_with_variant(&mut vnt_per_rec[i])?;
        }
        Ok(score + self.maximize_delta(vnt_per_rec))
    }
    fn maximize_delta<'v>(&self, delta: VntPerRec<'v>) -> f64 {
        let mut variants: SmallVec<[&VariantEval; 4]> = delta
            .iter()
            .flatten()
            .filter(|v| v.delta() > 0.0)
            .collect();

        if variants.is_empty() {
            return 0.0;
        }

        // Sort by end coordinate for Weighted Interval Scheduling
        variants.sort_unstable_by_key(|v| v.end());

        let mut dp = vec![0.0; variants.len()];
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

    fn score_with_variant(&mut self, vnt: &mut DeltaPerRec<'_>) -> Result<f64> {
        let mut score = 0.0;

        // Score variant bases that reach into segments before the current one.
        for prior_seg_i in 0..self.seg_i {
            if self.seg[prior_seg_i].is_last_in_template() != self.seg[self.seg_i].is_last_in_template() {
                continue; // skip read 1 read(s) when processing read 2
            }
            let seg_ref_start = self.seg[prior_seg_i].pos();
            let seg_ref_end = seg_ref_start + self.seg[prior_seg_i].cigar().len() as i64;

            // prior (or following) segment bases contribute to alt scoring only, no incurrence.
            self.score_variants_in_window(vnt, prior_seg_i, seg_ref_start, seg_ref_end, 0.0);
        }

        let op_iter = UnifiedOpIterator::new(self.seg[self.seg_i])
            .map_err(|e| anyhow!("Failed to create UnifiedOpIterator: {}", e))?;

        for op_res in op_iter {
            let op = op_res?;
            let ref_start = self.refpos;
            let mut ref_score = 0.0;
            match op {
                UnifiedOp::Match(len) => {
                    for _ in 0..len {
                        ref_score += self.pen.log_likelihood_match[self.q(self.seg_i, self.nt_i)];
                        self.nt_i += 1;
                        self.refpos += 1;
                    }
                },
                UnifiedOp::Mis(len) => {
                    for _ in 0..len {
                        ref_score += self.pen.log_likelihood_mismatch[self.q(self.seg_i, self.nt_i)];
                        self.nt_i += 1;
                        self.refpos += 1;
                    }
                }
                UnifiedOp::Del(len) => {
                    self.refpos += len as i64;
                    ref_score += self.pen.gap_open + (len as f64) * self.pen.gap_extend;
                }
                UnifiedOp::Ins(len) => {
                    self.nt_i += len as usize;
                    ref_score += self.pen.gap_open + (len as f64) * self.pen.gap_extend;
                },
                UnifiedOp::Relocate { penalty_score, pos } => {
                    self.refpos = pos;
                    ref_score += penalty_score;
                }
                UnifiedOp::RefSkip(len) => {
                    self.refpos += len as i64;
                }
            }

            self.score_variants_in_window(vnt, self.seg_i, ref_start, self.refpos, ref_score);
            score += ref_score;
        }

        for next_seg_i in (self.seg_i + 1)..self.seg.len() {
            if self.seg[next_seg_i].is_last_in_template() != self.seg[self.seg_i].is_last_in_template() {
                break; // skip segments from read 2 when processing read 1
            }
            let seg_ref_start = self.seg[next_seg_i].pos();
            let seg_ref_end = seg_ref_start + self.seg[next_seg_i].seq().len() as i64;

            self.score_variants_in_window(vnt, next_seg_i, seg_ref_start, seg_ref_end, 0.0);
        }

        Ok(score)
    }
    fn score_variants_in_window(
        &mut self,
        vnt: &mut DeltaPerRec<'_>,
        seg_i: usize,
        start: i64,
        end: i64,
        ref_score: f64,  // unweighted, for non-variant positions
    ) {
        let mut i = 0;
        while i < vnt.len() && vnt[i].start() < end {
            if let Some((weighted_ref_score, alt_score)) =
                self.score_variant_in_seg(&vnt[i], seg_i, start, end)
            {
                // Use weighted_ref_score instead of ref_score for variant positions
                vnt[i].update(weighted_ref_score, alt_score);
                let fully_processed = vnt[i].ref_end() <= end && vnt[i].alt_end() <= end;
                if fully_processed {
                    vnt.remove(i);
                    continue;
                }
            } else {
                vnt[i].update(ref_score, 0.0);
            }
            i += 1;
        }
    }
    /// Score a variant's alt allele against read bases from a specific segment.
    /// `ref_start` and `ref_end` define the reference window to score within this segment.
    fn score_variant_in_seg(
        &mut self,
        vnt_eval: &VariantEval,
        seg_i: usize,
        ref_start: i64,
        ref_end: i64,
    ) -> Option<(f64, f64)> {  // (weighted_ref_score, alt_score)
        let vnt_start = vnt_eval.start();

        // Clamp the ref window to the variant's ref span
        let eff_ref_start = ref_start.max(vnt_start);
        let eff_ref_end = ref_end.min(vnt_eval.ref_end());
        let ref_len = (eff_ref_end - eff_ref_start) as usize;

        // Derive the alt slice from the ref offset (valid for MNPs; see caveat below)
        let ref_consumed = (eff_ref_start - vnt_eval.start()) as usize;
        let vnt = vnt_eval.vnt();
        let alt = vnt.alt_allele();
        let alt_slice = &alt[ref_consumed.min(alt.len())..(ref_consumed + ref_len).min(alt.len())];
        if alt_slice.is_empty() && ref_len == 0 { return None; }

        let p_variant = vnt.p_variant();

        // Weighted ref score over the overlapping bases
        let mut weighted_ref_score = 0.0;
        for j in 0..(eff_ref_end - eff_ref_start) as usize {
            let nt_i = (eff_ref_start - self.seg[seg_i].pos()) as usize + j;
            let q = self.q(seg_i, nt_i);
            let lm  = self.pen.log_likelihood_match[q];
            let lmm = self.pen.log_likelihood_mismatch[q];
            // ref path: read matches ref, so weight by (1 - p_variant)
            weighted_ref_score += (1.0 - p_variant) * lm + p_variant * lmm;
        }

        let n = alt.len();
        let m = n.max(ref_len);

        self.prev.resize(m + 1, Cell::default());
        self.curr.resize(m + 1, Cell::default());

        // Initialise first row: gaps in the read (deletions relative to alt)
        self.prev[0].reinit(self.pen.gap_open, self.pen.gap_extend, 0);
        for j in 1..=m {
            self.prev[j].reinit(self.pen.gap_open, self.pen.gap_extend, -(j as i32));
        }

        let seg_ref_start = self.seg[seg_i].pos();
        let nt_i_base = (ref_start - seg_ref_start) as usize;

        for i in 1..=n {
            let alt_base = alt[i - 1];
            self.curr[0].reinit(self.pen.gap_open, self.pen.gap_extend, i as i32);

            for j in 1..=m {
                let mut nt_i = nt_i_base + j - 1;
                let read_base = if seg_i != self.seg_i && self.seg[seg_i].is_reverse() != self.seg[self.seg_i].is_reverse() {
                    let encoded = self.seg[seg_i].seq().encoded;
                    nt_i = encoded.len() - 1 - nt_i;
                    match encoded[nt_i] {
                        b'A' => b'T',
                        b'C' => b'G',
                        b'G' => b'C',
                        b'T' => b'A',
                        _ => b'N',
                    }
                } else {
                    self.seg[seg_i].seq().encoded[nt_i]
                };
                let q = self.q(seg_i, nt_i);

                let lm = self.pen.log_likelihood_match[q];
                let lmm = self.pen.log_likelihood_mismatch[q];

                // Weight match/mismatch score by variant probability,
                // mirroring score_alt_match / score_ref_match
                let per_base_score = if alt_base == read_base {
                    p_variant * lm + (1.0 - p_variant) * lmm
                } else {
                    (1.0 - p_variant) * lm + p_variant * lmm
                };

                self.curr[j].m = per_base_score
                    + self.prev[j - 1].m
                        .max(self.prev[j - 1].i)
                        .max(self.prev[j - 1].d);

                self.curr[j].i = (self.curr[j - 1].m + self.pen.gap_open + self.pen.gap_extend)
                    .max(self.curr[j - 1].i + self.pen.gap_extend);

                self.curr[j].d = (self.prev[j].m + self.pen.gap_open + self.pen.gap_extend)
                    .max(self.prev[j].d + self.pen.gap_extend);
            }

            std::mem::swap(&mut self.prev, &mut self.curr);
        }

        let alt_score = self.prev[m].m.max(self.prev[m].i).max(self.prev[m].d);

        Some((weighted_ref_score, alt_score))
    }
    fn q(&self, seg_i: usize, nt_i: usize) -> usize {
        (self.seg[seg_i].qual()[nt_i] as usize).min(MAX_Q - 1)
    }
}

#[cfg(test)]
pub(crate) mod tests;
