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
    fn matches_seq(&self, slices: SlicesForVariant, allele: &[u8]) -> bool {
        let mut offset = 0;

        for (is_revcmp, bases, _) in slices {
            let base_len = bases.len().min(allele.len().saturating_sub(offset));
            if is_revcmp {
                for b in bases[..base_len].iter().rev() {
                    match (b, allele[offset]) {
                        (b'A', b'T') | (b'T', b'A') | (b'C', b'G') | (b'G', b'C') => offset += 1,
                        _ => return false,
                    }
                }
            } else if bases[..base_len] != allele[offset..offset + base_len] {
                return false;
            }
            offset += base_len;
        }
        true
    }
    fn matches_alt(&self, slices: SlicesForVariant) -> bool {
        self.matches_seq(slices, self.alt_allele())

    }
    fn matches_ref(&self, slices: SlicesForVariant) -> bool {
        self.matches_seq(slices, self.ref_allele())
    }

    fn score_alt(&self, slices: SlicesForVariant, penalties: &Penalty, quals: &[u8]) -> f64 {
        let mut score = 0.0;
        let mut offset = 0;
        let alt = self.alt_allele();

        for (is_revcmp, bases, base_quals) in slices {
            let base_len = bases.len().min(alt.len().saturating_sub(offset));
            for i in 0..base_len {
                let q = (base_quals[i] as usize).min(MAX_Q - 1);
                score += if is_revcmp {
                    match (bases[i], alt[offset + i]) {
                        (b'A', b'T') | (b'T', b'A') | (b'C', b'G') | (b'G', b'C') => {
                            penalties.log_likelihood_match[q]
                        },
                        _ => penalties.log_likelihood_mismatch[q],
                    }
                } else if bases[i] == alt[offset + i] {
                    penalties.log_likelihood_match[q]
                } else {
                    penalties.log_likelihood_mismatch[q]
                };
            }
            offset += base_len;
        }
        score
    }

    /// Provides an adjusted score for a read chunk that **matches the ALT allele**.
    /// This is called when the read *disagrees* with the reference but *agrees* with the variant.
    fn score_alt_match(&self, penalties: &Penalty, quals: &[u8]) -> f64;

    /// Provides an adjusted score for a read chunk that **matches the REF allele**.
    /// This is called when a variant is present, but the read *agrees* with the reference.
    fn score_ref_match(&self, penalties: &Penalty, quals: &[u8]) -> f64;

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
