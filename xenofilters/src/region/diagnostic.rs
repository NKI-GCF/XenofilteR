//! In-memory store of species-diagnostic variant positions loaded from a VCF.
//!
//! A variant in stream N's diagnostic VCF means: a read aligned to stream N's
//! reference that overlaps this position carries evidence for species N.
//! Reads overlapping any diagnostic position must go through full scoring —
//! they cannot be early-assigned even if the alignment is otherwise perfect.
//!
//! The diagnostic allele is derived from MD + CIGAR in pass 1 (no sequence
//! needed). Confirmation that the read carries the diagnostic allele vs a
//! third allele is left to the full scorer.
//!
//! Store layout mirrors [`crate::variant::store::Store`]:
//! per-chromosome sorted `Vec` with `partition_point` binary search.

use anyhow::{anyhow, Result};
use noodles::bcf;
use noodles::vcf;
use std::collections::HashMap;
use std::path::Path;

/// A single diagnostic variant position.
#[derive(Debug, Clone)]
pub(crate) struct DiagnosticSite {
    /// 1-based reference position (matching VCF/BAM convention).
    pub(crate) pos: usize,
    /// Length of the reference allele (for overlap span).
    pub(crate) ref_len: usize,
    /// The diagnostic allele (alt). A read carrying this base at `pos`
    /// relative to the reference is evidence for this stream's species.
    pub(crate) alt: Vec<u8>,
}

impl DiagnosticSite {
    fn end(&self) -> usize {
        self.pos + self.ref_len
    }
}

/// Per-chromosome sorted list of diagnostic variant sites.
#[derive(Debug, Default)]
pub(crate) struct DiagnosticVariants {
    per_ref: HashMap<usize, Vec<DiagnosticSite>>,
    max_ref_len: usize,
}

impl DiagnosticVariants {
    /// Load from a BCF/VCF file. Reference sequence names are resolved to
    /// integer IDs using `name_to_id`.
    pub(crate) fn from_vcf(
        path: &Path,
        name_to_id: &HashMap<String, usize>,
    ) -> Result<Self> {
        let mut reader = bcf::io::reader::Builder::default()
            .build_from_path(path)
            .map_err(|e| anyhow!("Cannot open VCF/BCF {}: {}", path.display(), e))?;
        let header = reader.read_header()?;

        let mut per_ref: HashMap<usize, Vec<DiagnosticSite>> = HashMap::new();
        let mut max_ref_len = 1usize;

        for result in reader.records() {
            let record = result?;
            let chrom_idx = record.reference_sequence_id()?;

            // Resolve BCF reference sequence index to BAM header ref_id.
            let chrom_name = header
                .contigs()
                .get_index(chrom_idx)
                .map(|(name, _)| name.to_string())
                .ok_or_else(|| anyhow!("BCF contig index {chrom_idx} not in header"))?;
            let ref_id = match name_to_id.get(&chrom_name) {
                Some(&id) => id,
                None => continue,
            };

            let pos = record
                .variant_start()
                .transpose()?
                .map(|p| p.get())
                .unwrap_or(0);

            let ref_allele = record.reference_bases();
            let ref_bytes = ref_allele.as_ref();
            let ref_len = ref_bytes.len();

            if ref_len > max_ref_len {
                max_ref_len = ref_len;
            }

            // Each ALT allele becomes a separate diagnostic site.
            let alts = record.alternate_bases();
            for alt_bytes in alts.as_ref().split(|&b| b == b',') {
                per_ref.entry(ref_id).or_default().push(DiagnosticSite {
                    pos,
                    ref_len,
                    alt: alt_bytes.to_vec(),
                });
            }
        }

        for sites in per_ref.values_mut() {
            sites.sort_unstable_by_key(|s| s.pos);
        }

        Ok(Self { per_ref, max_ref_len })
    }

    /// Returns `true` if any diagnostic site overlaps `[read_start, read_end)`.
    /// `read_start` and `read_end` are 1-based, matching BAM alignment_start.
    pub(crate) fn overlaps(&self, ref_id: usize, read_start: usize, read_end: usize) -> bool {
        let Some(sites) = self.per_ref.get(&ref_id) else {
            return false;
        };
        // First site that could overlap: pos >= read_start - max_ref_len.
        let lower = read_start.saturating_sub(self.max_ref_len);
        let first = sites.partition_point(|s| s.pos < lower);
        sites[first..].iter().any(|s| s.pos < read_end && s.end() > read_start)
    }

    pub(crate) fn is_empty(&self) -> bool {
        self.per_ref.is_empty()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn site(pos: usize, ref_len: usize) -> DiagnosticSite {
        DiagnosticSite { pos, ref_len, alt: b"A".to_vec() }
    }

    fn store(sites: &[(usize, usize, usize)]) -> DiagnosticVariants {
        let mut per_ref: HashMap<usize, Vec<DiagnosticSite>> = HashMap::new();
        let mut max_ref_len = 1;
        for &(rid, pos, ref_len) in sites {
            if ref_len > max_ref_len { max_ref_len = ref_len; }
            per_ref.entry(rid).or_default().push(site(pos, ref_len));
        }
        for v in per_ref.values_mut() {
            v.sort_unstable_by_key(|s| s.pos);
        }
        DiagnosticVariants { per_ref, max_ref_len }
    }

    #[test]
    fn test_snv_overlap() {
        let s = store(&[(0, 100, 1)]);
        assert!(s.overlaps(0, 100, 110));
    }

    #[test]
    fn test_snv_no_overlap_before() {
        let s = store(&[(0, 100, 1)]);
        assert!(!s.overlaps(0, 50, 100));
    }

    #[test]
    fn test_snv_no_overlap_after() {
        let s = store(&[(0, 100, 1)]);
        assert!(!s.overlaps(0, 101, 200));
    }

    #[test]
    fn test_deletion_overlap_span() {
        // Deletion at pos 100, ref_len 5 → spans [100, 105)
        let s = store(&[(0, 100, 5)]);
        assert!(s.overlaps(0, 103, 110));
        assert!(!s.overlaps(0, 105, 200));
    }

    #[test]
    fn test_wrong_ref() {
        let s = store(&[(0, 100, 1)]);
        assert!(!s.overlaps(1, 100, 110));
    }
}
