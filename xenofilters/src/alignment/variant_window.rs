//! Geometry + scoring helpers for evaluating a candidate variant against a
//! window of aligned read. Split out of `fragment.rs` so each piece can be
//! unit-tested without a full Fragment/Scratch/BAM-record harness.

use crate::filter_algorithm::line_by_line::Scratch;
use crate::penalty::Penalty;
use crate::Error;

/// Geometry of how a `[ref_start, ref_end)` scoring window overlaps a
/// candidate variant's reference span `[vnt_start, vnt_ref_end)`.
#[derive(Debug, Clone, Copy, PartialEq)]
pub(crate) struct VariantWindow {
    pub(crate) eff_ref_start: usize,
    pub(crate) eff_ref_end: usize,
    pub(crate) ref_len: usize,
    /// Offset into the variant's own ref/alt allele where this window begins.
    pub(crate) ref_consumed: usize,
}

impl VariantWindow {
    /// `None` if the window and the variant's ref span don't actually overlap.
    pub(crate) fn compute(
        ref_start: usize,
        ref_end: usize,
        vnt_start: usize,
        vnt_ref_end: usize,
    ) -> Option<Self> {
        let eff_ref_start = ref_start.max(vnt_start);
        let eff_ref_end = ref_end.min(vnt_ref_end);
        if eff_ref_start >= eff_ref_end {
            return None;
        }
        Some(VariantWindow {
            eff_ref_start,
            eff_ref_end,
            ref_len: eff_ref_end - eff_ref_start,
            ref_consumed: eff_ref_start - vnt_start,
        })
    }

    /// 0-based read offset corresponding to `eff_ref_start` (assumes no
    /// indels between the segment's own start and this window — same
    /// assumption the original code made).
    pub(crate) fn read_offset(&self, seg_ref_start: usize) -> usize {
        self.eff_ref_start - seg_ref_start
    }
}

/// Expected score for this window if we assume "maybe still reference,
/// maybe the variant", weighted by `p_variant`.
pub(crate) fn weighted_ref_score(
    window: VariantWindow,
    seg_ref_start: usize,
    p_variant: f64,
    pen: &Penalty,
    mut quality_at: impl FnMut(usize) -> Result<usize, Error>,
) -> Result<f64, Error> {
    let read_offset = window.read_offset(seg_ref_start);
    let mut score = 0.0;
    for j in 0..window.ref_len {
        let q = quality_at(read_offset + j)?;
        score += (1.0 - p_variant) * pen.log_likelihood_match[q]
            + p_variant * pen.log_likelihood_mismatch[q];
    }
    #[cfg(test)]
    eprintln!(
        "[weighted_ref_score] window={window:?} seg_ref_start={seg_ref_start} \
         read_offset={read_offset} p_variant={p_variant} -> {score}"
    );
    Ok(score)
}

/// Gapped alignment of `alt` against the read starting at 0-based
/// `read_offset`, weighted by `p_variant`. `read_base_and_quality(j)` must
/// return the read base and quality-array index for read offset `j`
/// (handles revcomp internally on the caller side).
pub(crate) fn align_alt_to_read(
    alt: &[u8],
    read_offset: usize,
    ref_len: usize,
    p_variant: f64,
    pen: &Penalty,
    read_base_and_quality: impl Fn(usize) -> Result<(u8, usize), Error>,
    scratch: &mut Scratch,
) -> Result<f64, Error> {
    let n = alt.len();
    let m = n.max(ref_len);

    scratch.resize_nw(m + 1);
    scratch.prev[0].reinit(pen.gap_open, pen.gap_extend, 0);
    for (j, cell) in scratch.prev.iter_mut().enumerate().take(m + 1).skip(1) {
        cell.reinit(pen.gap_open, pen.gap_extend, -(j as i32));
    }

    for i in 1..=n {
        let alt_base = alt[i - 1];
        scratch.curr[0].reinit(pen.gap_open, pen.gap_extend, i as i32);

        for j in 1..=m {
            let (read_base, q) = read_base_and_quality(read_offset + j - 1)?;
            let lm = pen.log_likelihood_match[q];
            let lmm = pen.log_likelihood_mismatch[q];
            let per_base_score = if alt_base == read_base {
                p_variant * lm + (1.0 - p_variant) * lmm
            } else {
                (1.0 - p_variant) * lm + p_variant * lmm
            };

            #[cfg(test)]
            eprintln!(
                "[align_alt_to_read] i={i} j={j} read_idx={} alt={} read={} match={} per_base={per_base_score}",
                read_offset + j - 1,
                alt_base as char,
                read_base as char,
                alt_base == read_base
            );

            scratch.curr[j].m = per_base_score
                + scratch.prev[j - 1]
                    .m
                    .max(scratch.prev[j - 1].i)
                    .max(scratch.prev[j - 1].d);
            scratch.curr[j].i = (scratch.curr[j - 1].m + pen.gap_open + pen.gap_extend)
                .max(scratch.curr[j - 1].i + pen.gap_extend);
            scratch.curr[j].d = (scratch.prev[j].m + pen.gap_open + pen.gap_extend)
                .max(scratch.prev[j].d + pen.gap_extend);
        }
        scratch.swap_nw();
    }

    let alt_score = scratch.prev[m]
        .m
        .max(scratch.prev[m].i)
        .max(scratch.prev[m].d);
    #[cfg(test)]
    eprintln!("[align_alt_to_read] alt_score={alt_score}");
    Ok(alt_score)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::penalty::MAX_Q;

    fn pen() -> Penalty {
        Penalty {
            gap_open: -2.0,
            gap_extend: -0.5,
            log_likelihood_match: [0.0; MAX_Q],
            log_likelihood_mismatch: [-1.0; MAX_Q],
            chimeric_junction_penalty: -3.0,
            log_likelihood_bisulfite: 0.0,
        }
    }

    // --- VariantWindow ---

    #[test]
    fn test_window_basic_overlap() {
        let w = VariantWindow::compute(1, 6, 3, 4).unwrap();
        assert_eq!(w.eff_ref_start, 3);
        assert_eq!(w.eff_ref_end, 4);
        assert_eq!(w.ref_len, 1);
        assert_eq!(w.ref_consumed, 0);
    }

    #[test]
    fn test_window_read_offset_for_snv_in_5m_segment() {
        // This is exactly the snp_alt_support_gives_positive_delta scenario:
        // segment starts at ref pos 1, variant at ref pos 3 -> read index 2.
        let w = VariantWindow::compute(1, 6, 3, 4).unwrap();
        assert_eq!(w.read_offset(1), 2);
    }

    #[test]
    fn test_window_no_overlap() {
        // no ovelap before
        assert!(VariantWindow::compute(1, 3, 5, 6).is_none());
        // touching boundary is not overlap
        assert!(VariantWindow::compute(1, 3, 3, 6).is_none());
        // no overlap after
    }

    #[test]
    fn test_window_partial_offset_into_multibase_variant() {
        let w = VariantWindow::compute(5, 10, 3, 8).unwrap(); // vnt span [3,8)
        assert_eq!(w.eff_ref_start, 5);
        assert_eq!(w.ref_len, 3);
        assert_eq!(w.ref_consumed, 2);
    }

    // --- weighted_ref_score ---

    #[test]
    fn test_weighted_ref_score_certain_variant_is_full_mismatch_penalty() -> Result<(), Error> {
        let window = VariantWindow::compute(1, 6, 3, 4).unwrap();
        let score = weighted_ref_score(window, 1, 1.0, &pen(), |_| Ok(30))?;
        assert!((score - (-1.0)).abs() < 1e-12);
        Ok(())
    }

    #[test]
    fn test_weighted_ref_score_certain_reference_is_zero() -> Result<(), Error> {
        let window = VariantWindow::compute(1, 6, 3, 4).unwrap();
        let score = weighted_ref_score(window, 1, 0.0, &pen(), |_| Ok(30))?;
        assert!(score.abs() < 1e-12);
        Ok(())
    }

    #[test]
    fn test_weighted_ref_score_queries_correct_read_index() -> Result<(), Error> {
        let window = VariantWindow::compute(1, 6, 3, 4).unwrap();
        let mut seen = Vec::new();
        weighted_ref_score(window, 1, 0.5, &pen(), |nt_i| {
            seen.push(nt_i);
            Ok(30)
        })?;
        assert_eq!(seen, vec![2]); // must be read index 2, not 0
        Ok(())
    }

    // --- align_alt_to_read ---

    #[test]
    fn test_align_certain_alt_match_scores_zero() -> Result<(), Error> {
        let mut scratch = Scratch::new();
        let read = b"AAGAA";
        let score = align_alt_to_read(
            b"G",
            2,
            1,
            1.0,
            &pen(),
            |idx| Ok((read[idx], 30)),
            &mut scratch,
        )?;
        assert!(score.abs() < 1e-9);
        Ok(())
    }

    #[test]
    fn test_align_wrong_offset_reproduces_the_original_bug() -> Result<(), Error> {
        // Same read/allele, but offset 0 instead of 2 — this is exactly what
        // the pre-fix `nt_i_base = ref_start - seg_ref_start` produced.
        let mut scratch = Scratch::new();
        let read = b"AAGAA";
        let score = align_alt_to_read(
            b"G",
            0,
            1,
            1.0,
            &pen(),
            |idx| Ok((read[idx], 30)),
            &mut scratch,
        )?;
        assert!((score - (-1.0)).abs() < 1e-9); // read[0]='A' vs alt 'G' -> mismatch
        Ok(())
    }
}
