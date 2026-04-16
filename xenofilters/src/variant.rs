mod population_variant;
mod sample_variant;
mod variant_store;

pub(super) use population_variant::{PopulationVariant, parse_population_record};
pub(super) use sample_variant::{SampleVariant, parse_sample_record};
pub(super) use variant_store::{VariantStore, VariantStoreTrait};
use crate::Penalty;
use anyhow::{Result, anyhow};
use rust_htslib::bcf::{self, Read};
use std::path::Path;
use crate::alignment::{Segment,AlignmentError};
use crate::alignment::{UnifiedOp, UnifiedOpIterator};
use crate::MAX_Q;
use smallvec::SmallVec;

type SliceInfo<'a> = (bool, &'a [u8], &'a [u8]); // (is_revcmp, bases, quals)
type SlicesForVariant<'a> = SmallVec<[SliceInfo<'a>; 1]>;

fn reverse_complement_encoded_slice(bases: &[u8]) -> Vec<u8> {
    bases.iter().rev().map(|&b| match b {
        b'A' => b'T',
        b'T' => b'A',
        b'C' => b'G',
        b'G' => b'C',
        _ => b'N',
    }).collect()
}

/// Trait for any object that can be scored against an alignment.
pub trait Variant {
    /// The 1-based reference position of the variant.
    fn pos(&self) -> i64;

    /// The reference allele
    fn ref_allele(&self) -> &[u8];

    /// The alternate allele
    fn alt_allele(&self) -> &[u8];

    fn get_read_window_len(&self) -> usize {
        // may be longer than the variant itself due to indels.
        self.alt_allele().len().max(self.ref_allele().len())
    }
    fn get_size(&self) -> usize {
        let n = self.alt_allele().len();
        let m = self.get_read_window_len();
        (n + 1) * (m + 1)
    }

    /// Checks if a read chunk matches this variant's ALT allele.
    /// Default implementation handles simple string equality.
    fn matches_alt(&self, read_bases: &[u8]) -> bool {
        self.alt_allele() == read_bases
    }
    fn matches_ref(&self, read_bases: &[u8]) -> bool {
        self.ref_allele() == read_bases
    }

    /// Provides an adjusted score for a read chunk that **matches the ALT allele**.
    /// This is called when the read *disagrees* with the reference but *agrees* with the variant.
    fn score_alt_match(&self, penalties: &Penalty, quals: &[u8]) -> f64;

    /// Provides an adjusted score for a read chunk that **matches the REF allele**.
    /// This is called when a variant is present, but the read *agrees* with the reference.
    fn score_ref_match(&self, penalties: &Penalty, quals: &[u8]) -> f64;

    fn len(&self) -> usize { self.ref_allele().len().max(self.alt_allele().len()) }

    /// End position (1-based, exclusive) — allows easy overlap check
    fn end(&self) -> i64 {
        self.pos() + self.ref_allele().len() as i64
    }

    /// Does this variant overlap [read_start, read_end) ?
    fn overlaps(&self, read_start: i64, read_end: i64) -> bool {
        let v_start = self.pos();
        let v_end   = self.end();
        v_start < read_end && v_end > read_start
    }
    fn extract_context(
        &self,
        segments: &[Segment],
        penalty: &Penalty
    ) -> Result<(Vec<u8>, Vec<u8>, f64), AlignmentError> {
        let mut bases = Vec::new();
        let mut quals = Vec::new();
        let mut incurred = 0.0;

        let v_start = self.pos();
        let v_end = self.end(); // exclusive

        for seg in segments {
            // Optimization: Skip segment entirely if it doesn't overlap variant bounds
            if seg.ref_end <= v_start || seg.ref_start >= v_end {
                continue;
            }

            let op_iter = UnifiedOpIterator::new(seg.rec)
                .map_err(|e| anyhow!("Failed to create UnifiedOpIterator: {}", e))?;

            let mut ref_pos = seg.ref_start;
            let mut offset = 0;

            for op_res in op_iter {
                let op = op_res?;

                match op {
                    UnifiedOp::Match(len) | UnifiedOp::Mis(len) => {
                        for i in 0..len as usize {
                            let curr_pos = ref_pos + i as i64;
                            if curr_pos >= v_start && curr_pos < v_end {
                                let q = seg.rec.qual()[offset + i];
                                let b = seg.rec.seq().encoded[offset + i];
                                bases.push(b);
                                quals.push(q);

                                // Reverse exactly what was added in Phase 1
                                if let UnifiedOp::Mis(_) = op {
                                    incurred += penalty.log_likelihood_mismatch[(q as usize).min(MAX_Q - 1)];
                                } else {
                                    incurred += penalty.log_likelihood_match[(q as usize).min(MAX_Q - 1)];
                                }
                            }
                        }
                        ref_pos += len as i64;
                        offset += len as usize;
                    }
                    UnifiedOp::Ins(len) => {
                        // Insertions consume 0 ref bases, so they sit exactly at ref_pos
                        if ref_pos >= v_start && ref_pos < v_end {
                            for i in 0..len as usize {
                                let q = seg.rec.qual()[offset + i];
                                let b = seg.rec.seq().encoded[offset + i];
                                bases.push(b);
                                quals.push(q);
                            }
                            incurred += penalty.gap_open + (len as f64 - 1.0) * penalty.gap_extend;
                        }
                        offset += len as usize;
                    }
                    UnifiedOp::Del(len) => {
                        let end_pos = ref_pos + len as i64;
                        // If deletion overlaps the variant window
                        if ref_pos < v_end && end_pos > v_start {
                            incurred += penalty.gap_open + (len as f64 - 1.0) * penalty.gap_extend;
                        }
                        ref_pos += len as i64;
                    }
                    UnifiedOp::RefSkip(len) => ref_pos += len as i64,
                    UnifiedOp::Relocate { pos, .. } => ref_pos = pos,
                }
            }
        }

        Ok((bases, quals, incurred))
    }

}

/// Generic VCF reader that accepts a parser function.
pub fn vcf_reader<V: Variant>(
    f: &Path,
    parser: impl Fn(&mut bcf::Record) -> Result<Vec<V>>,
) -> Result<VariantStore<V>> {
    let mut bcf_reader = bcf::Reader::from_path(f)
        .map_err(|e| anyhow!("Failed to open VCF/BCF {}: {}", f.display(), e))?;

    let mut per_chr: Vec<Vec<V>> = Vec::new();
    let mut max_variant_len: i64 = 1; // at least 1
    let mut is_sorted = true;

    for record_result in bcf_reader.records() {
        let mut record = record_result?;
        let rid = record.rid().ok_or_else(|| anyhow!("Record missing RID"))? as usize;
        let variants = parser(&mut record)?;

        while per_chr.len() <= rid {
            per_chr.push(Vec::new());
        }
        let mut last_pos = None;

        for v in variants {
            let pos = v.pos();
            if is_sorted {
                if let Some(last) = last_pos {
                    if pos < last {
                        eprintln!("Variants in {} are not sorted by position.", f.display());
                        is_sorted = false;
                    }
                }
                last_pos = Some(pos);
            }
            let span = v.end() - pos;
            if span > max_variant_len {
                max_variant_len = span;
            }
            per_chr[rid].push(v);
        }
    }
    if !is_sorted {
        // Sort each chromosome once, for binary search.
        for chr in &mut per_chr {
            chr.sort_by_key(|v| v.pos());
        }
    }

    Ok(VariantStore {
        per_chr,
        max_variant_len,
    })
}
