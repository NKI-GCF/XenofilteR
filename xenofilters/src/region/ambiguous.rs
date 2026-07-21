//! In-memory store of ambiguous genomic regions loaded from a BED file.
//!
//! A read whose alignment overlaps any region here cannot be early-assigned.
//! Regions are stored per reference sequence (keyed by ref-id), sorted by
//! start position, and queried via binary search.

use crate::region::interval_store::{load_into_store, Interval, IntervalStore};
use crate::Error;
use noodles::bed;
use std::collections::HashMap;
use std::fs::File;
use std::io::BufReader;
use std::path::Path;

/// A half-open interval `[start, end)` on a single reference sequence.
/// Coordinates are 0-based, matching BED convention.
#[derive(Debug, Clone, Copy)]
pub(crate) struct Region {
    pub(crate) start: usize,
    pub(crate) end: usize,
}

impl Interval for Region {
    fn start(&self) -> usize {
        self.start
    }
    fn end(&self) -> usize {
        self.end
    }
}

/// Per-chromosome sorted list of ambiguous regions.
#[derive(Debug, Default)]
pub(crate) struct AmbiguousRegions {
    store: IntervalStore<Region>,
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

        // bed::io::Reader has no Iterator impl in noodles; adapt with from_fn.
        let records = std::iter::from_fn(|| match reader.read_record(&mut record) {
            Ok(0) => None,
            Ok(_) => Some(Ok(record.clone())),
            Err(e) => Some(Err(Error::BedParseError {
                path: path.to_path_buf(),
                source: e,
            })),
        });

        let store = load_into_store(
            records,
            name_to_id,
            |rec: &bed::Record<3>| {
                String::try_from(rec.reference_sequence_name()).unwrap_or_default()
            },
            |rec, _ref_id| {
                let start = rec
                    .feature_start()
                    .map_err(|e| Error::BedStartError(e.to_string()))?;
                let end = rec
                    .feature_end()
                    .ok_or(Error::BedRecordMissingEnd)?
                    .map_err(|e| Error::BedEndError(e.to_string()))?;
                Ok(Some(Region {
                    start: usize::from(start),
                    end: usize::from(end),
                }))
            },
        )?;
        Ok(Self { store })
    }

    /// Returns `true` if `[read_start, read_end)` overlaps any ambiguous region.
    pub(crate) fn overlaps(&self, ref_id: usize, start: usize, end: usize) -> bool {
        self.store.overlaps(ref_id, start, end)
    }

    pub(crate) fn is_empty(&self) -> bool {
        self.store.is_empty()
    }

    #[cfg(test)]
    pub(crate) fn from_test(per_ref: HashMap<usize, Vec<Region>>) -> Self {
        let mut store = IntervalStore::new();
        for (rid, regions) in per_ref {
            for r in regions {
                store.insert(rid, r);
            }
        }
        store.sort();
        Self { store }
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

        let cases = [
            (0, 120, 150, true),  // contained
            (0, 50, 150, true),   // straddles left
            (0, 150, 250, true),  // straddles right
            (0, 50, 100, false),  // touching left edge
            (0, 200, 300, false), // touching right edge
            (1, 100, 200, false), // wrong ref_id
        ];

        for (q_ref, q_start, q_end, expected) in cases {
            assert_eq!(s.overlaps(q_ref, q_start, q_end), expected);
        }
    }

    #[test]
    fn test_multiple_regions() {
        let s = store(&[(0, 10, 20), (0, 100, 200), (0, 300, 400)]);
        assert!(s.overlaps(0, 150, 250));
        assert!(!s.overlaps(0, 200, 300));
        assert!(s.overlaps(0, 15, 25));
    }
}
