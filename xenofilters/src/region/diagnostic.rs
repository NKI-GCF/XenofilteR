//! In-memory store of species-diagnostic variant positions loaded from a VCF.
//!
//! A variant in stream N's diagnostic VCF means: a read aligned to stream N's
//! reference overlapping this position carries evidence for species N.
//! Reads overlapping any diagnostic position must go through full scoring.

use crate::Error;
use noodles::bcf;
use std::collections::HashMap;
use std::path::Path;

#[derive(Debug, Clone)]
pub(crate) struct DiagnosticSite {
    pub(crate) pos: usize,
    pub(crate) ref_len: usize,
}

impl DiagnosticSite {
    fn end(&self) -> usize {
        self.pos + self.ref_len
    }
}

#[derive(Debug, Default)]
pub(crate) struct DiagnosticVariants {
    per_ref: HashMap<usize, Vec<DiagnosticSite>>,
    max_ref_len: usize,
}

impl DiagnosticVariants {
    pub(crate) fn from_vcf(
        path: &Path,
        name_to_id: &HashMap<String, usize>,
    ) -> Result<Self, Error> {
        let mut reader = bcf::io::reader::Builder::default()
            .build_from_path(path)
            .map_err(|e| Error::CannotOpenVcfBcf {
                path: path.to_path_buf(),
                source: e,
            })?;
        let header = reader.read_header()?;

        let mut per_ref: HashMap<usize, Vec<DiagnosticSite>> = HashMap::new();
        let mut max_ref_len = 1usize;

        for result in reader.records() {
            let record = result?;
            let chrom_idx = record.reference_sequence_id()?;

            let chrom_name = header
                .contigs()
                .get_index(chrom_idx)
                .map(|(name, _)| name.to_string())
                .ok_or(Error::BcfContigMissing { chrom_idx })?;
            let ref_id = match name_to_id.get(&chrom_name) {
                Some(&id) => id,
                None => continue,
            };

            let pos = record
                .variant_start()
                .transpose()?
                .map(|p| p.get())
                .unwrap_or(0);

            // Bind to a variable to extend the lifetime of the temporary.
            let ref_bases = record.reference_bases();
            let ref_bytes: &[u8] = ref_bases.as_ref();
            let ref_len = ref_bytes.len().max(1);
            if ref_len > max_ref_len {
                max_ref_len = ref_len;
            }

            per_ref
                .entry(ref_id)
                .or_default()
                .push(DiagnosticSite { pos, ref_len });
        }

        for sites in per_ref.values_mut() {
            sites.sort_unstable_by_key(|s| s.pos);
        }

        Ok(Self {
            per_ref,
            max_ref_len,
        })
    }

    pub(crate) fn overlaps(&self, ref_id: usize, read_start: usize, read_end: usize) -> bool {
        let Some(sites) = self.per_ref.get(&ref_id) else {
            return false;
        };
        let lower = read_start.saturating_sub(self.max_ref_len);
        let first = sites.partition_point(|s| s.pos < lower);
        sites[first..]
            .iter()
            .any(|s| s.pos < read_end && s.end() > read_start)
    }

    pub(crate) fn is_empty(&self) -> bool {
        self.per_ref.is_empty()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn site(pos: usize, ref_len: usize) -> DiagnosticSite {
        DiagnosticSite { pos, ref_len }
    }

    fn store(sites: &[(usize, usize, usize)]) -> DiagnosticVariants {
        let mut per_ref: HashMap<usize, Vec<DiagnosticSite>> = HashMap::new();
        let mut max_ref_len = 1;
        for &(rid, pos, ref_len) in sites {
            if ref_len > max_ref_len {
                max_ref_len = ref_len;
            }
            per_ref.entry(rid).or_default().push(site(pos, ref_len));
        }
        for v in per_ref.values_mut() {
            v.sort_unstable_by_key(|s| s.pos);
        }
        DiagnosticVariants {
            per_ref,
            max_ref_len,
        }
    }

    #[test]
    fn test_snv_overlap() {
        assert!(store(&[(0, 100, 1)]).overlaps(0, 100, 110));
    }

    #[test]
    fn test_snv_no_overlap_before() {
        assert!(!store(&[(0, 100, 1)]).overlaps(0, 50, 100));
    }

    #[test]
    fn test_snv_no_overlap_after() {
        assert!(!store(&[(0, 100, 1)]).overlaps(0, 101, 200));
    }

    #[test]
    fn test_deletion_span() {
        let s = store(&[(0, 100, 5)]);
        assert!(s.overlaps(0, 103, 110));
        assert!(!s.overlaps(0, 105, 200));
    }

    #[test]
    fn test_wrong_ref() {
        assert!(!store(&[(0, 100, 1)]).overlaps(1, 100, 110));
    }
}
