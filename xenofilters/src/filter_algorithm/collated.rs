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

use crate::{Error, print_routing_counters};
use crate::alignment::{Fragment, FragmentState, MdCigFlags, SimpleRec, stringify_record,PreAssessResult, pre_assess_alignments, mate_slot, segment_id, score_state_nw, compute_cancel_slots};
use crate::aln_stream::AlignmentStream;
use crate::config::{Config, StripReadSuffix};
use crate::filter_algorithm::line_by_line::{READ_CT, Scratch, ordering::Decision};
use crate::penalty::Penalty;
use crate::region::tabix_query::{TabixBed, TabixVcf};
use crate::variant::FragEvalVec;
use noodles::sam::alignment::record::{cigar::Cigar, data::field::Tag};
use noodles::sam::alignment::record_buf::{RecordBuf, data::field::Value};
use reader::{CollatedReader, canonical_name};
use smallvec::SmallVec;
use std::collections::HashMap;
use crate::alignment::{MateClassifiable, MateKind};

pub(crate) struct CollatedMatcher<R: SimpleRec> {
    a: CollatedReader<R>,
    b: CollatedReader<R>,
    waiting_a: HashMap<Box<[u8]>, FragmentState<R>>,
    waiting_b: HashMap<Box<[u8]>, FragmentState<R>>,
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
    pub(crate) fn new(
        config: Config,
        mut aln: SmallVec<[Box<dyn AlignmentStream<R>>; 2]>,
        bed: [Option<TabixBed>; 2],
        vcf: [Option<TabixVcf>; 2],
    ) -> Result<Self, Error> {
        assert_eq!(
            aln.len(),
            2,
            "CollatedMatcher requires exactly 2 alignment streams"
        );
        let ambiguous_log_threshold = match config.ambiguous_threshold {
            0 => 0.0,
            t => (t as f64) * std::f64::consts::LN_10 / 10.0,
        };
        let penalties = config.to_penalties();
        let strip = config.strip_read_suffix;
        let add_decision_tag = config.add_decision_tag;

        let aln_len = aln.len();
        for (i, a) in aln.iter_mut().enumerate() {
            a.init_writers(&config, i)?;
        }

        let mut iter = aln.into_iter();
        let stream0 = iter.next().unwrap();
        let stream1 = iter.next().unwrap();

        Ok(Self {
            a: CollatedReader::new(stream0, strip, 0),
            b: CollatedReader::new(stream1, strip, 1),
            waiting_a: HashMap::new(),
            waiting_b: HashMap::new(),
            penalties,
            scratch: Scratch::new(),
            routing_counters: SmallVec::from_elem(0, aln_len * 4),
            add_decision_tag,
            ambiguous_log_threshold,
            strip,
            bed,
            vcf,
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
    pub(crate) fn process(&mut self) -> Result<(), Error> {
        loop {
            let fa = self.a.next_fragment()?;
            let fb = self.b.next_fragment()?;
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
        print_routing_counters(&self.routing_counters, "collated");
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
        if self.bed[aln_idx].is_none() && self.vcf[aln_idx].is_none() {
            return Ok(false);
        }
        let header = if aln_idx == 0 {
            self.a.header()
        } else {
            self.b.header()
        };
        // Build ref_id → name map once per call (headers are small).
        let ref_names: Vec<String> = header
            .reference_sequences()
            .iter()
            .map(|(name, _)| name.to_string())
            .collect();

        for rec in frag.get_records() {
            let flags = rec.flags()?;
            if flags.is_secondary() || flags.is_unmapped() {
                continue;
            }
            let ref_id = match rec.ref_seq_id().transpose()? {
                Some(id) => id,
                None => continue,
            };
            let chrom = match ref_names.get(ref_id) {
                Some(n) => n.as_str(),
                None => continue,
            };
            let start = match rec.alignment_start().transpose()? {
                Some(p) => p.get().saturating_sub(1), // 0-based for BED
                None => continue,
            };
            let end = start + rec.cigar().len();

            if let Some(bed) = &mut self.bed[aln_idx]
                && bed.overlaps(chrom, start, end)?
            {
                return Ok(true);
            }
            if let Some(vcf) = &mut self.vcf[aln_idx]
                && vcf.overlaps(chrom, start, end)?
            {
                return Ok(true);
            }
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
        let cancel_slot = compute_cancel_slots(&a, &b);

        // Tier 3: full NW scoring.
        let s1 = score_state_nw(&a, mcfs1, self.a.variant_store().as_deref(), &self.penalties, &mut self.scratch, cancel_slot)?;
        let s2 = score_state_nw(&b, mcfs2, self.b.variant_store().as_deref(), &self.penalties, &mut self.scratch, cancel_slot)?;
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
