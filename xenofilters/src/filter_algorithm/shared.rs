//! Scoring helpers shared by `collated` and `hash_lookup` backends.
//!
//! Both backends need to:
//! 1. Take a pair of `FragmentState<R>` (one per stream).
//! 2. Resolve which is better (or ambiguous).
//! 3. Emit records to the appropriate output writers.
//!
//! This mirrors `LineByLine::resolve` / `apply_delta` / `handle_best` but is
//! decoupled from the streaming merge loop.

use crate::alignment::{FragmentState, MdCigFlags, SimpleRec};
use crate::aln_stream::AlignmentStream;
use crate::filter_algorithm::line_by_line::{Scratch, READ_CT};
use crate::penalty::Penalty;
use anyhow::{anyhow, Result};
use noodles::sam::alignment::record_buf::RecordBuf;
use noodles::sam::alignment::RecordBuf as RecordBufAlias;
use smallvec::SmallVec;
use std::cmp::Ordering;

// Re-use the Decision type from ordering so CLI tag logic is consistent.
pub(crate) use crate::filter_algorithm::line_by_line::ordering::Decision;

/// Minimal context the shared scorer needs from the owning backend struct.
pub(crate) struct ScorerCtx<'a, R: SimpleRec> {
    pub(crate) aln: &'a mut SmallVec<[Box<dyn AlignmentStream<R>>; 2]>,
    pub(crate) penalties: &'a Penalty,
    pub(crate) scratch: &'a mut Scratch,
    pub(crate) branch_counters: &'a mut [u64; 32],
    pub(crate) add_decision_tag: bool,
    pub(crate) ambiguous_log_threshold: f64,
}

enum Resolution {
    Ordered(Ordering),
    Scored(f64),
}

/// Score a single `FragmentState` against its alignment stream, returning the
/// log-likelihood sum (same as `LineByLine::score_candidate`).
fn score_candidate<R: SimpleRec>(
    ctx: &mut ScorerCtx<R>,
    state: &FragmentState<R>,
    mcfs: SmallVec<[MdCigFlags<'_>; READ_CT]>,
    aln_idx: usize,
) -> Result<f64> {
    use crate::alignment::{Fragment, stringify_record};
    use crate::variant::FragEvalVec;

    let mut segment: SmallVec<[&R; READ_CT]> = SmallVec::new();
    let mut seg_mcfs: SmallVec<[MdCigFlags; READ_CT]> = SmallVec::new();
    let mut mcfs_opt: SmallVec<[Option<MdCigFlags>; READ_CT]> =
        mcfs.into_iter().map(Some).collect();

    let aln = ctx
        .aln
        .get(aln_idx)
        .ok_or_else(|| anyhow!("No alignment for index {aln_idx}"))?;
    let mut dvnt: FragEvalVec<'_> = SmallVec::new();

    for idx in state.order_mates() {
        let rec = &state.get_records()[idx];
        let flags = state
            .flags(idx)
            .ok_or_else(|| anyhow!("No flags for record {idx}"))?;
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
                .ok_or_else(|| anyhow!("MdCigFlags missing for record {idx}"))?
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
            seg_mcfs.push(
                mcfs_opt[idx]
                    .take()
                    .ok_or_else(|| anyhow!("MdCigFlags already consumed for record {idx}"))?,
            );
        } else if flags.is_last_segment() {
            break;
        }
    }

    Fragment::new(ctx.penalties, segment, seg_mcfs)?
        .score(ctx.scratch, &mut dvnt)
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
        })
}

fn resolve<R: SimpleRec>(
    ctx: &mut ScorerCtx<R>,
    a: &FragmentState<R>,
    b: &FragmentState<R>,
) -> Result<Resolution> {
    let mut ord = a.partial_cmp(b);
    if ord.is_none() {
        let (mcfs1, mcfs2) = a.cmp_perfect(b, &mut ord)?;
        if ord.is_none() {
            let score1 = score_candidate(ctx, a, mcfs1, a.get_nr())?;
            let score2 = score_candidate(ctx, b, mcfs2, b.get_nr())?;
            return Ok(Resolution::Scored(score1 - score2));
        }
    }
    Ok(Resolution::Ordered(ord.expect("must be Some")))
}

/// Write `frag`'s records to filtered output (lost this comparison).
fn emit_filtered<R: SimpleRec>(
    ctx: &mut ScorerCtx<R>,
    frag: &mut FragmentState<R>,
) -> Result<()> {
    let nr = frag.get_nr();
    for r in frag.drain_records() {
        let header = ctx.aln[nr].header();
        let rec = r.as_record_buf(header)?;
        ctx.branch_counters[nr << 1] += 1;
        if let Some(a) = ctx.aln.get_mut(nr) {
            a.write_record(rec, Some(false))?;
        }
    }
    Ok(())
}

/// Write `frag`'s records to best / ambiguous output.
fn emit_best<R: SimpleRec>(
    ctx: &mut ScorerCtx<R>,
    frag: &mut FragmentState<R>,
    decision: Option<&Decision>,
    is_ambiguous: bool,
) -> Result<()> {
    let nr = frag.get_nr();
    let best_state = if is_ambiguous { None } else { Some(true) };
    for r in frag.drain_records() {
        let header = ctx.aln[nr].header();
        let mut rec: RecordBuf = r.as_record_buf(header)?;
        if !is_ambiguous {
            match decision {
                Some(Decision::ConfDelta(v)) => {
                    let tag = noodles::sam::alignment::record::data::field::Tag::new(b'X', b'F');
                    rec.data_mut()
                        .insert(tag, noodles::sam::alignment::record_buf::data::field::Value::from(*v));
                }
                Some(Decision::VariantRescued(v)) => {
                    let tag = noodles::sam::alignment::record::data::field::Tag::new(b'X', b'R');
                    rec.data_mut()
                        .insert(tag, noodles::sam::alignment::record_buf::data::field::Value::from(*v));
                }
                _ => {}
            }
            ctx.branch_counters[1 + (nr << 1)] += 1;
        } else {
            ctx.branch_counters[16 + nr] += 1;
        }
        if let Some(a) = ctx.aln.get_mut(nr) {
            a.write_record(rec, best_state)?;
        }
    }
    Ok(())
}

/// Score a matched pair and emit to the correct outputs.
/// `a` is always stream-0 fragment, `b` is always stream-1 fragment.
pub(crate) fn score_and_emit<R: SimpleRec>(
    ctx: &mut ScorerCtx<R>,
    mut a: FragmentState<R>,
    mut b: FragmentState<R>,
) -> Result<()> {
    match resolve(ctx, &a, &b)? {
        Resolution::Ordered(Ordering::Greater) => {
            emit_filtered(ctx, &mut b)?;
            let dec = ctx.add_decision_tag.then_some(Decision::First);
            emit_best(ctx, &mut a, dec.as_ref(), false)?;
        }
        Resolution::Ordered(Ordering::Less) => {
            emit_filtered(ctx, &mut a)?;
            let dec = ctx.add_decision_tag.then_some(Decision::Last);
            emit_best(ctx, &mut b, dec.as_ref(), false)?;
        }
        Resolution::Ordered(Ordering::Equal) => {
            let dec = ctx.add_decision_tag.then_some(Decision::Ambiguous);
            emit_best(ctx, &mut a, dec.as_ref(), true)?;
            emit_best(ctx, &mut b, dec.as_ref(), true)?;
        }
        Resolution::Scored(delta) => {
            let (decision, winner, loser) = if delta > ctx.ambiguous_log_threshold {
                let phred = if ctx.add_decision_tag {
                    let p = (10.0 * delta / std::f64::consts::LN_10) as u32;
                    Some(Decision::ConfDelta(p.min(255) as u8))
                } else {
                    None
                };
                (phred, &mut a as *mut _, &mut b as *mut _)
            } else if delta < -ctx.ambiguous_log_threshold {
                let phred = if ctx.add_decision_tag {
                    let p = (10.0 * (-delta) / std::f64::consts::LN_10) as u32;
                    Some(Decision::ConfDelta(p.min(255) as u8))
                } else {
                    None
                };
                (phred, &mut b as *mut _, &mut a as *mut _)
            } else {
                // Ambiguous
                let dec = ctx.add_decision_tag.then_some(Decision::Ambiguous);
                emit_best(ctx, &mut a, dec.as_ref(), true)?;
                emit_best(ctx, &mut b, dec.as_ref(), true)?;
                return Ok(());
            };
            // SAFETY: winner and loser are disjoint borrows of a and b.
            let (winner, loser) = unsafe { (&mut *winner, &mut *loser) };
            emit_filtered(ctx, loser)?;
            emit_best(ctx, winner, decision.as_ref(), false)?;
        }
    }
    Ok(())
}

/// Emit a fragment for which no match was found in the other stream.
/// Treated as: this stream wins (written to best output), other stream
/// has nothing to filter. This is a safety net for malformed input.
pub(crate) fn emit_unmatched<R: SimpleRec>(
    ctx: &mut ScorerCtx<R>,
    mut frag: FragmentState<R>,
) -> Result<()> {
    emit_best(ctx, &mut frag, None, false)
}
