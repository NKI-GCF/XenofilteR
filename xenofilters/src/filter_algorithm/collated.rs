//! [`CollatedMatcher`] — fragment-matching for individually-collated BAM streams
//! whose inter-stream name order may differ.
//!
//! Each input stream must have all records for a given read name contiguous
//! (collated / query-name-grouped), but the two streams need not present
//! fragments in the same order. A `HashMap` on each side buffers unmatched
//! fragments until their counterpart arrives.
//!
//! BED ambiguous regions and diagnostic VCF positions are queried via tabix
//! for early-assignment decisions.
//!
//! Memory usage: O(name-order skew). Output order: not guaranteed.
//! Single-threaded only.

pub(crate) mod reader;
#[cfg(test)]
pub(crate) mod tests;

use crate::Error;
use crate::alignment::{FragmentState, SimpleRec, PreAssessResult, pre_assess_alignments, MateKind, mate_kind::MateClassifiable };
use crate::aln_stream::AlignmentStream;
use crate::config::{CommonArgs, StripReadSuffix};
use crate::filter_algorithm::line_by_line::{Scratch, ordering::Decision};
use crate::penalty::Penalty;
use crate::region::tabix_query::{TabixBed, TabixVcf};
use noodles::sam::alignment::record::data::field::Tag;
use noodles::sam::alignment::record_buf::{RecordBuf, data::field::Value};
use reader::{CollatedReader, canonical_name};
use smallvec::SmallVec;
use std::collections::HashMap;
use ahash::RandomState;
use crate::config::CollatedArgs;
use crate::config::args::resolve_threshold;

pub(crate) struct CollatedMatcher<R: SimpleRec> {
    a: CollatedReader<R>,
    b: CollatedReader<R>,
    waiting_a: HashMap<Box<[u8]>, FragmentState<R>, RandomState>,
    waiting_b: HashMap<Box<[u8]>, FragmentState<R>, RandomState>,
    penalties: Penalty,
    scratch: Scratch,
    pub(crate) routing_counters: SmallVec<[u64; 8]>,
    add_decision_tag: bool,
    ambiguous_log_threshold: f64,
    strip: StripReadSuffix,
    bed: [Option<TabixBed>; 2],
    vcf: [Option<TabixVcf>; 2],
}

impl<R: SimpleRec> CollatedMatcher<R> {
    pub(crate) fn new_from_collated(
        args: &CollatedArgs,
        aln:  SmallVec<[Box<dyn AlignmentStream<RecordBuf>>; 2]>,
        bed:  [Option<TabixBed>; 2],
        vcf:  [Option<TabixVcf>; 2],
    ) -> Result<CollatedMatcher<RecordBuf>, Error> {
        let ambiguous_log_threshold = resolve_threshold(
            args.common.scoring.ambiguous_threshold, false,
        );
        let mut it = aln.into_iter();
        let a0 = it.next().unwrap();
        let a1 = it.next().unwrap();
        let bisulfite = args.common.scoring.bisulfite; 
        Ok(CollatedMatcher {
            a: CollatedReader::new(a0, bisulfite, 0),
            b: CollatedReader::new(a1, bisulfite, 1),
            waiting_a: HashMap::with_hasher(RandomState::new()),
            waiting_b: HashMap::with_hasher(RandomState::new()),
            penalties: args.common.scoring.to_penalty(),
            scratch: Scratch::new(),
            routing_counters: SmallVec::from_elem(0, 2 * 4), // 2 streams × 4 counters per stream
            add_decision_tag: args.common.io.add_decision_tag,
            ambiguous_log_threshold,
            bed,
            vcf,
            strip: args.common.io.strip_read_suffix,
        })
    }

    // CONCURRENCY STUB — CollatedMatcher parallel worker pool
    //
    // `score_pair` / `nw_score_fragment` are embarrassingly parallel once a pair
    // has been extracted from `waiting_a` / `waiting_b`.  A crossbeam bounded
    // channel can dispatch pairs to N workers:
    //
    //   let (work_tx, work_rx) = bounded::<(FragmentState<R>, FragmentState<R>)>(N*2);
    //   Workers call score_pair and send ScoredFragment back to the IO thread.
    //   Writers remain on the IO thread (no Mutex needed).
    //
    // Output order is NOT guaranteed (acceptable for Collated).
    // N-STREAM: scales to N waiting maps; memory is O(name-order skew × streams).
    pub(crate) fn process(&mut self, config: &CommonArgs) -> Result<(), Error> {
        loop {
            let fa = self.a.next_fragment(config.io.strip_read_suffix)?;
            let fb = self.b.next_fragment(config.io.strip_read_suffix)?;
            match (fa, fb) {
                (None, None) => break,
                (Some(fa), Some(fb)) => {
                    self.handle_fragment(fa)?;
                    self.handle_fragment(fb)?;
                }
                (Some(fa), None) => self.handle_fragment(fa)?,
                (None, Some(fb)) => self.handle_fragment(fb)?,
            }
        }
        // Safety net: valid input leaves both maps empty.
        let drain_a: Vec<_> = self.waiting_a.drain().map(|(_, v)| v).collect();
        let drain_b: Vec<_> = self.waiting_b.drain().map(|(_, v)| v).collect();
        for frag in drain_a {
            self.write_fragment(frag, None, Some(true))?;
        }
        for frag in drain_b {
            self.write_fragment(frag, None, Some(true))?;
        }
        config.print_routing_counters(&self.routing_counters, "Collated");
        Ok(())
    }

    fn handle_fragment(&mut self, frag: FragmentState<R>) -> Result<(), Error> {
        let nr = frag.get_nr();
        let key = canonical_name(frag.first_qname(), self.strip);
        if nr == 0 {
            if let Some(other) = self.waiting_b.remove(key.as_ref()) {
                self.score_pair(frag, other)?;
            } else {
                self.waiting_a.insert(key, frag);
            }
        } else if let Some(other) = self.waiting_a.remove(key.as_ref()) {
            self.score_pair(other, frag)?;
        } else {
            self.waiting_b.insert(key, frag);
        }
        Ok(())
    }
    /// Returns true if any primary alignment in `frag` overlaps a BED ambiguous
    /// region or diagnostic VCF position for `aln_idx`. Forces full NW scoring
    /// even when `cmp_perfect` would give a quick answer.
    fn fragment_overlaps_region(
        &mut self,
        frag: &FragmentState<R>,
        aln_idx: usize,
    ) -> Result<bool, Error> {
        let bed = match &self.bed[aln_idx] { Some(b) => b, None => return Ok(false) };
        for (rec, flags) in frag.get_records().iter().zip(frag.flags.iter()) {
            if flags.is_secondary() || flags.is_unmapped() { continue; }
            let is_rev  = flags.is_reverse_complemented();
            let ref_id  = rec.ref_seq_id().transpose()?.ok_or(Error::NoRefSeqId)?;
            let start   = rec.alignment_start().transpose()?.ok_or(Error::NoAlignmentStart)?
                             .get().saturating_sub(1);
            let end     = start + rec.cigar().as_ref().len();
            // `any_overlap` now strand-aware via BED column 6.
            if bed.any_overlap(ref_id, start, end, is_rev)? { return Ok(true); }
        }
        Ok(false)
    }

    fn score_pair(&mut self, a: FragmentState<R>, b: FragmentState<R>) -> Result<(), Error> {
        // Tier 1: unmapped fast-path — before BED/VCF I/O.
        let mut ord = a.partial_cmp(&b);
        if let Some(o) = ord {
            return self.apply_ordered(a, b, o);
        }

        let a_needs_scoring = self.fragment_overlaps_region(&a, 0)?;
        let b_needs_scoring = self.fragment_overlaps_region(&b, 1)?;

        let (mcfs1, mcfs2) = a.cmp_perfect(&b, &mut ord)?;

        // Tier 2: perfect-match fast-path (no region forces scoring).
        if let Some(o) = ord
            && !a_needs_scoring
            && !b_needs_scoring
        {
            drop(mcfs1);
            drop(mcfs2);
            return self.apply_ordered(a, b, o);
        }

        // Tier 2.5: unified pre-assessment — single CIGAR+MD walk per record.
        // Guard: only when no BED/VCF region forces full scoring (diagnostic variants
        // must be scored via NW to properly account for variant rescue).
        if !a_needs_scoring
            && !b_needs_scoring
            && let PreAssessResult::EarlyDecision(pa_ord) = pre_assess_alignments(&mcfs1, &mcfs2)
        {
            drop(mcfs1);
            drop(mcfs2);
            return self.apply_ordered(a, b, pa_ord);
        }

        // 2-way mate cancellation: a mate slot cancels only when BOTH streams
        // agree on the same non-Other kind. Region overlap on either stream
        // forces full scoring of every mate (variant rescue must see it).
        let cancel_slot = if a_needs_scoring || b_needs_scoring {
            [false, false]
        } else {
            let mk_a = a.mate_kinds();
            let mk_b = b.mate_kinds();
            [
                matches!((mk_a[0], mk_b[0]), (Some(x), Some(y)) if x == y && x != MateKind::Other),
                matches!((mk_a[1], mk_b[1]), (Some(x), Some(y)) if x == y && x != MateKind::Other),
            ]
        };

        // Tier 3: full NW scoring.
        let s1 = a.score_state_nw(mcfs1, self.a.variant_store().as_deref(), &self.penalties, &mut self.scratch, cancel_slot, self.a.positive_regions())?;
        let s2 = b.score_state_nw(mcfs2, self.b.variant_store().as_deref(), &self.penalties, &mut self.scratch, cancel_slot, self.b.positive_regions())?;
        let delta = s1 - s2;

        if delta > self.ambiguous_log_threshold {
            let dec = self.phred_decision(delta);
            self.emit_discarded(b)?;
            self.write_fragment(a, dec.as_ref(), Some(true))?;
        } else if delta < -self.ambiguous_log_threshold {
            let dec = self.phred_decision(-delta);
            self.emit_discarded(a)?;
            self.write_fragment(b, dec.as_ref(), Some(true))?;
        } else {
            let dec = self.add_decision_tag.then_some(Decision::Ambiguous);
            self.write_fragment(a, dec.as_ref(), None)?;
            self.write_fragment(b, dec.as_ref(), None)?;
        }
        Ok(())
    }

    fn apply_ordered(
        &mut self,
        a: FragmentState<R>,
        b: FragmentState<R>,
        ord: std::cmp::Ordering,
    ) -> Result<(), Error> {
        use std::cmp::Ordering::{Equal, Greater, Less};
        match ord {
            Greater => {
                self.emit_discarded(b)?;
                let dec = self.add_decision_tag.then_some(Decision::First);
                self.write_fragment(a, dec.as_ref(), Some(true))?;
            }
            Less => {
                self.emit_discarded(a)?;
                let dec = self.add_decision_tag.then_some(Decision::Last);
                self.write_fragment(b, dec.as_ref(), Some(true))?;
            }
            Equal => {
                let dec = self.add_decision_tag.then_some(Decision::Ambiguous);
                self.write_fragment(a, dec.as_ref(), None)?;
                self.write_fragment(b, dec.as_ref(), None)?;
            }
        }
        Ok(())
    }

    fn phred_decision(&self, abs_delta: f64) -> Option<Decision> {
        self.add_decision_tag.then(|| {
            let p = (10.0 * abs_delta / std::f64::consts::LN_10) as u32;
            Decision::PhredConfidence(p.min(255) as u8)
        })
    }

    fn emit_discarded(&mut self, mut frag: FragmentState<R>) -> Result<(), Error> {
        let nr = frag.get_nr();
        let base = nr * 4;
        for r in frag.drain_records() {
            let stream = if nr == 0 { &mut self.a } else { &mut self.b };
            let rec = r.as_record_buf(stream.header())?;
            self.routing_counters[base] += 1;
            stream.write_record(rec, Some(false))?;
        }
        Ok(())
    }

    fn write_fragment(
        &mut self,
        mut frag: FragmentState<R>,
        decision: Option<&Decision>,
        best_state: Option<bool>,
    ) -> Result<(), Error> {
        let nr = frag.get_nr();
        let base = nr * 4;
        let is_ambiguous = best_state.is_none();
        for r in frag.drain_records() {
            let stream = if nr == 0 { &mut self.a } else { &mut self.b };
            let mut rec: RecordBuf = r.as_record_buf(stream.header())?;
            if !is_ambiguous {
                match decision {
                    Some(Decision::PhredConfidence(v)) => {
                        rec.data_mut().insert(Tag::new(b'X', b'F'), Value::from(*v));
                    }
                    Some(Decision::VariantRescued(v)) => {
                        rec.data_mut().insert(Tag::new(b'X', b'R'), Value::from(*v));
                    }
                    _ => {}
                }
                self.routing_counters[base + 1] += 1;
            } else {
                self.routing_counters[base + 2] += 1;
            }
            stream.write_record(rec, best_state)?;
        }
        Ok(())
    }
}
