//! [`CollatedMatcher`] — fragment-matching for individually-collated BAM streams
//! whose inter-stream name order may differ.
//!
//! Each input stream must have all records for a given read name contiguous,
//! but the two streams need not present fragments in the same order.
//! A `HashMap` on each side buffers unmatched fragments.
//!
//! BED ambiguous regions and diagnostic VCF positions are queried via tabix
//! for early-assignment decisions (random access needed because the BAM is
//! name-ordered, not position-ordered).
//!
//! Memory usage: O(name-order skew). Output order: not guaranteed.
//! Single-threaded only.

pub(crate) mod reader;
#[cfg(test)]
pub(crate) mod tests;

use crate::alignment::{Fragment, FragmentState, MdCigFlags, SimpleRec, stringify_record};
use crate::aln_stream::AlignmentStream;
use crate::config::{Config, StripReadSuffix};
use crate::filter_algorithm::line_by_line::{Scratch, READ_CT, ordering::Decision};
use crate::penalty::Penalty;
use crate::region::tabix_query::{TabixBed, TabixVcf};
use crate::variant::FragEvalVec;
use anyhow::{anyhow, Result};
use noodles::sam::alignment::record::data::field::Tag;
use noodles::sam::alignment::record_buf::data::field::Value;
use noodles::sam::alignment::record_buf::RecordBuf;
use reader::{CollatedReader, canonical_name};
use smallvec::SmallVec;
use std::collections::HashMap;

pub(crate) struct CollatedMatcher<R: SimpleRec> {
    a: CollatedReader<R>,
    b: CollatedReader<R>,
    waiting_a: HashMap<Box<[u8]>, FragmentState<R>>,
    waiting_b: HashMap<Box<[u8]>, FragmentState<R>>,
    penalties: Penalty,
    scratch: Scratch,
    pub(crate) branch_counters: [u64; 32],
    add_decision_tag: bool,
    ambiguous_log_threshold: f64,
    strip: StripReadSuffix,
    /// Optional tabix-indexed BED for early-assignment checks, per stream.
    bed: [Option<TabixBed>; 2],
    /// Optional tabix-indexed VCF for early-assignment checks, per stream.
    vcf: [Option<TabixVcf>; 2],
}

impl<R: SimpleRec> CollatedMatcher<R> {
    pub(crate) fn new(
        config: Config,
        mut aln: SmallVec<[Box<dyn AlignmentStream<R>>; 2]>,
        bed: [Option<TabixBed>; 2],
        vcf: [Option<TabixVcf>; 2],
    ) -> Result<Self> {
        assert_eq!(aln.len(), 2, "CollatedMatcher requires exactly 2 alignment streams");
        let ambiguous_log_threshold = match config.ambiguous_threshold {
            0 => 0.0,
            t => (t as f64) * std::f64::consts::LN_10 / 10.0,
        };
        let penalties = config.to_penalties();
        let strip = config.strip_read_suffix;
        let add_decision_tag = config.add_decision_tag;

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
            branch_counters: [0u64; 32],
            add_decision_tag,
            ambiguous_log_threshold,
            strip,
            bed,
            vcf,
        })
    }

    /// Convenience constructor with no BED/VCF (tests and namesorted-only paths).
    pub(crate) fn new_no_regions(
        config: Config,
        aln: SmallVec<[Box<dyn AlignmentStream<R>>; 2]>,
    ) -> Result<Self> {
        Self::new(config, aln, [None, None], [None, None])
    }

    pub(crate) fn process(&mut self) -> Result<()> {
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
        for frag in drain_a { self.emit_records_owned(frag, None, Some(true))?; }
        for frag in drain_b { self.emit_records_owned(frag, None, Some(true))?; }
        self.print_counters();
        Ok(())
    }

    fn handle_fragment(&mut self, frag: FragmentState<R>) -> Result<()> {
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

    fn score_pair(&mut self, a: FragmentState<R>, b: FragmentState<R>) -> Result<()> {
        let mut ord = a.partial_cmp(&b);

        // By nesting this block, mcfs1 and mcfs2 are dropped before `res` is matched,
        // avoiding borrow-of-moved-value errors on a and b.
        enum Res { Ordered(std::cmp::Ordering), Scored(f64) }
        let res = if ord.is_none() {
            let (mcfs1, mcfs2) = a.cmp_perfect(&b, &mut ord)?;
            if ord.is_none() {
                let s1 = self.score_one(&a, mcfs1, 0)?;
                let s2 = self.score_one(&b, mcfs2, 1)?;
                Res::Scored(s1 - s2)
            } else {
                Res::Ordered(ord.unwrap())
            }
        } else {
            Res::Ordered(ord.unwrap())
        };

        use std::cmp::Ordering::*;
        match res {
            Res::Ordered(Greater) => {
                self.emit_filtered(b)?;
                let dec = self.add_decision_tag.then_some(Decision::First);
                self.emit_records_owned(a, dec.as_ref(), Some(true))?;
            }
            Res::Ordered(Less) => {
                self.emit_filtered(a)?;
                let dec = self.add_decision_tag.then_some(Decision::Last);
                self.emit_records_owned(b, dec.as_ref(), Some(true))?;
            }
            Res::Ordered(Equal) => {
                let dec = self.add_decision_tag.then_some(Decision::Ambiguous);
                self.emit_records_owned(a, dec.as_ref(), None)?;
                self.emit_records_owned(b, dec.as_ref(), None)?;
            }
            Res::Scored(delta) if delta > self.ambiguous_log_threshold => {
                let dec = self.phred_decision(delta);
                self.emit_filtered(b)?;
                self.emit_records_owned(a, dec.as_ref(), Some(true))?;
            }
            Res::Scored(delta) if delta < -self.ambiguous_log_threshold => {
                let dec = self.phred_decision(-delta);
                self.emit_filtered(a)?;
                self.emit_records_owned(b, dec.as_ref(), Some(true))?;
            }
            Res::Scored(_) => {
                let dec = self.add_decision_tag.then_some(Decision::Ambiguous);
                self.emit_records_owned(a, dec.as_ref(), None)?;
                self.emit_records_owned(b, dec.as_ref(), None)?;
            }
        }
        Ok(())
    }

    fn phred_decision(&self, abs_delta: f64) -> Option<Decision> {
        self.add_decision_tag.then(|| {
            let p = (10.0 * abs_delta / std::f64::consts::LN_10) as u32;
            Decision::ConfDelta(p.min(255) as u8)
        })
    }

    fn score_one(
        &mut self,
        state: &FragmentState<R>,
        mcfs: SmallVec<[MdCigFlags<'_>; READ_CT]>,
        aln_idx: usize,
    ) -> Result<f64> {
        let mut segment: SmallVec<[&R; READ_CT]> = SmallVec::new();
        let mut seg_mcfs: SmallVec<[MdCigFlags; READ_CT]> = SmallVec::new();
        let mut mcfs_opt: SmallVec<[Option<MdCigFlags>; READ_CT]> =
            mcfs.into_iter().map(Some).collect();
        let mut dvnt: FragEvalVec<'_> = SmallVec::new();

        for idx in state.order_mates() {
            let rec = &state.get_records()[idx];
            let flags = state.flags(idx)
                .ok_or_else(|| anyhow!("No flags for record {idx}"))?;
            if flags.is_unmapped() {
                dvnt.push(SmallVec::new());
            } else {
                let tid = rec.ref_seq_id().transpose()?
                    .ok_or_else(|| anyhow!("No reference sequence ID"))?;
                let start = rec.alignment_start().transpose()?
                    .ok_or_else(|| anyhow!("No alignment start"))?.get();
                let cig_len = mcfs_opt[idx].as_ref()
                    .ok_or_else(|| anyhow!("MdCigFlags missing for {idx}"))?.get_cigar().len();
                let end = start + cig_len;
                let delta_vars = {
                    let stream = if aln_idx == 0 { &self.a } else { &self.b };
                    match stream.variant_store() {
                        Some(store) => store.overlapping_multi(tid, start, end),
                        None => SmallVec::new(),
                    }
                };
                dvnt.push(delta_vars);
            }
            if !flags.is_secondary() {
                segment.push(rec);
                seg_mcfs.push(mcfs_opt[idx].take()
                    .ok_or_else(|| anyhow!("MdCigFlags consumed for {idx}"))?);
            } else if flags.is_last_segment() {
                break;
            }
        }

        Fragment::new(&self.penalties, segment, seg_mcfs)?
            .score(&mut self.scratch, &mut dvnt)
            .map_err(|e| anyhow!(
                "Error scoring fragment for alignment {aln_idx}: {}\n{}",
                e,
                state.get_records().iter().map(stringify_record).collect::<Vec<_>>().join("\n")
            ))
    }

    fn emit_filtered(&mut self, mut frag: FragmentState<R>) -> Result<()> {
        let nr = frag.get_nr();
        for r in frag.drain_records() {
            let stream = if nr == 0 { &mut self.a } else { &mut self.b };
            let rec = r.as_record_buf(stream.header())?;
            self.branch_counters[nr << 1] += 1;
            stream.write_record(rec, Some(false))?;
        }
        Ok(())
    }

    fn emit_records_owned(
        &mut self,
        mut frag: FragmentState<R>,
        decision: Option<&Decision>,
        best_state: Option<bool>,
    ) -> Result<()> {
        let nr = frag.get_nr();
        let is_ambiguous = best_state.is_none();
        for r in frag.drain_records() {
            let stream = if nr == 0 { &mut self.a } else { &mut self.b };
            let mut rec: RecordBuf = r.as_record_buf(stream.header())?;
            if !is_ambiguous {
                match decision {
                    Some(Decision::ConfDelta(v)) => {
                        rec.data_mut().insert(Tag::new(b'X', b'F'), Value::from(*v));
                    }
                    Some(Decision::VariantRescued(v)) => {
                        rec.data_mut().insert(Tag::new(b'X', b'R'), Value::from(*v));
                    }
                    _ => {}
                }
                self.branch_counters[1 + (nr << 1)] += 1;
            } else {
                self.branch_counters[16 + nr] += 1;
            }
            stream.write_record(rec, best_state)?;
        }
        Ok(())
    }

    pub(crate) fn print_counters(&self) {
        for i in 0..2 {
            eprintln!("collated[filter:{}]: {}", i, self.branch_counters[i << 1]);
            eprintln!("collated[out:{}]: {}", i, self.branch_counters[1 + (i << 1)]);
            eprintln!("collated[ambig:{}]: {}", i, self.branch_counters[16 + i]);
        }
    }
}
