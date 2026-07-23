//! In-memory store of species-diagnostic variant positions loaded from a VCF.
//!
//! A variant in stream N's diagnostic VCF means: a read aligned to stream N's
//! reference overlapping this position carries evidence for species N.
//! Reads overlapping any diagnostic position must go through full scoring.

use crate::Error;
use crate::variant::store::load_lappers;
use noodles::bcf;
use rust_lapper::{Interval, Lapper};
use std::collections::HashMap;
use std::path::Path;

#[derive(Debug, Clone, Eq, PartialEq)]
pub(crate) struct DiagnosticSite {
    pub(crate) pos: usize,
    pub(crate) ref_len: usize,
}

#[derive(Debug, Default)]
pub(crate) struct SegregateVariants {
    pub(crate) per_chr: HashMap<usize, Lapper<usize, DiagnosticSite>>,
}

impl SegregateVariants {
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

        let per_chr = load_lappers(
            reader.records().map(|r| r.map_err(Error::from)),
            name_to_id,
            |record| {
                record
                    .reference_sequence_id()
                    .ok()
                    .and_then(|idx| header.contigs().get_index(idx).map(|(n, _)| n.to_string()))
                    .unwrap_or_default()
            },
            |record, _ref_id| {
                let pos = record
                    .variant_start()
                    .transpose()?
                    .map(|p| p.get())
                    .unwrap_or(0);
                let ref_len = record.reference_bases().as_ref().len().max(1);
                let site = DiagnosticSite { pos, ref_len };
                Ok(Some(Interval {
                    start: pos,
                    stop: pos + ref_len,
                    val: site,
                }))
            },
        )?;

        Ok(Self { per_chr })
    }

    pub(crate) fn overlaps(&self, ref_id: usize, start: usize, end: usize) -> bool {
        self.per_chr
            .get(&ref_id)
            .is_some_and(|lapper| lapper.find(start, end).next().is_some())
    }

    pub(crate) fn is_empty(&self) -> bool {
        self.per_chr.values().all(|lapper| lapper.is_empty())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn site(pos: usize, ref_len: usize) -> DiagnosticSite {
        DiagnosticSite { pos, ref_len }
    }
    fn store(sites: &[(usize, usize, usize)]) -> SegregateVariants {
        let mut raw_intervals: HashMap<usize, Vec<Interval<usize, DiagnosticSite>>> =
            HashMap::new();
        for &(rid, pos, ref_len) in sites {
            raw_intervals.entry(rid).or_default().push(Interval {
                start: pos,
                stop: pos + ref_len,
                val: DiagnosticSite { pos, ref_len },
            });
        }
        let per_chr = raw_intervals
            .into_iter()
            .map(|(rid, ivs)| (rid, Lapper::new(ivs)))
            .collect();
        SegregateVariants { per_chr }
    }

    #[test]
    fn test_snv_overlap_cases() {
        let s = store(&[(0, 100, 1)]);

        // Format: (q_ref, q_start, q_end, expected_overlap)
        let cases = [
            (0, 100, 110, true),  // overlaps
            (0, 50, 100, false),  // no overlap before
            (0, 101, 200, false), // no overlap after
            (1, 100, 110, false), // wrong ref
        ];

        for (q_ref, q_start, q_end, expected) in cases {
            assert_eq!(s.overlaps(q_ref, q_start, q_end), expected);
        }
    }

    #[test]
    fn test_deletion_span() {
        let s = store(&[(0, 100, 5)]);
        assert!(s.overlaps(0, 103, 110));
        assert!(!s.overlaps(0, 105, 200));
    }
}
