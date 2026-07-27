//! In-memory store of ambiguous genomic regions loaded from a BED file.
//!
//! A read whose alignment overlaps any region here cannot be early-assigned.
//! Regions are stored per reference sequence (keyed by ref-id), sorted by
//! start position, and queried via binary search.

use crate::Error;
use crate::variant::store::load_lappers;
use noodles::bed;
use rust_lapper::{Interval, Lapper};
use std::collections::HashMap;
use std::fs::File;
use std::io::BufReader;
use std::path::Path;

/// A half-open interval `[start, end)` on a single reference sequence.
/// Coordinates are 0-based, matching BED convention.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) struct Region {
    pub(crate) start: usize,
    pub(crate) end: usize,
}

/// Per-chromosome sorted list of ambiguous regions.
#[derive(Debug, Default)]
pub(crate) struct AmbiguousRegions {
    per_chr: HashMap<usize, Lapper<usize, Region>>,
}

impl AmbiguousRegions {
    /// Load from a BED file. Reference sequence names are resolved to integer
    /// IDs via `name_to_id` (derived from the BAM header).
    pub(crate) fn from_bed(
        path: &Path,
        name_to_id: &HashMap<String, usize>,
    ) -> Result<Self, Error> {
        let file = File::open(path).map_err(|e| Error::CannotOpenBedFile {
            path: path.to_path_buf(),
            source: e,
        })?;
        let mut reader = bed::io::Reader::<3, _>::new(BufReader::new(file));
        let mut record = bed::Record::<3>::default();

        let records = std::iter::from_fn(|| match reader.read_record(&mut record) {
            Ok(0) => None,
            Ok(_) => Some(Ok(record.clone())),
            Err(e) => Some(Err(Error::BedParseError {
                path: path.to_path_buf(),
                source: e,
            })),
        });

        let per_chr = load_lappers(
            records,
            name_to_id,
            |rec: &bed::Record<3>| {
                String::try_from(rec.reference_sequence_name()).unwrap_or_default()
            },
            |rec, _ref_id| {
                let start = usize::from(
                    rec.feature_start()
                        .map_err(|e| Error::BedStartError(e.to_string()))?,
                );
                let end = usize::from(
                    rec.feature_end()
                        .ok_or(Error::BedRecordMissingEnd)?
                        .map_err(|e| Error::BedEndError(e.to_string()))?,
                );

                Ok(Some(Interval {
                    start,
                    stop: end,
                    val: Region { start, end },
                }))
            },
        )?;

        Ok(Self { per_chr })
    }

    /// Returns `true` if `[read_start, read_end)` overlaps any ambiguous region.
    pub(crate) fn overlaps(&self, ref_id: usize, start: usize, end: usize) -> bool {
        self.per_chr
            .get(&ref_id)
            .is_some_and(|lapper| lapper.find(start, end).next().is_some())
    }

    pub(crate) fn is_empty(&self) -> bool {
        self.per_chr.values().all(|lapper| lapper.is_empty())
    }

    #[cfg(test)]
    pub(crate) fn from_test(per_ref: HashMap<usize, Vec<Region>>) -> Self {
        let per_chr = per_ref
            .into_iter()
            .map(|(rid, regions)| {
                let ivs = regions
                    .into_iter()
                    .map(|r| Interval {
                        start: r.start,
                        stop: r.end,
                        val: r,
                    })
                    .collect();
                (rid, Lapper::new(ivs))
            })
            .collect();
        Self { per_chr }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn store(regions: &[(usize, usize, usize)]) -> AmbiguousRegions {
        let mut per_ref: HashMap<usize, Vec<Region>> = HashMap::new();
        for &(rid, start, end) in regions {
            per_ref.entry(rid).or_default().push(Region { start, end });
        }
        AmbiguousRegions::from_test(per_ref)
    }

    #[test]
    fn test_overlap_cases() {
        let s = store(&[(0, 100, 200)]);
        crate::tests::common::run_overlap_cases(
            &[(0, 120, 150, true), (0, 50, 150, true), (0, 150, 250, true),
              (0, 50, 100, false), (0, 200, 300, false), (1, 100, 200, false)],
            |r, s2, e| s.overlaps(r, s2, e),
        );
    }

    #[test]
    fn test_multiple_regions() {
        let s = store(&[(0, 10, 20), (0, 100, 200), (0, 300, 400)]);
        assert!(s.overlaps(0, 150, 250));
        assert!(!s.overlaps(0, 200, 300));
        assert!(s.overlaps(0, 15, 25));
    }
}
