//! In-memory store of ambiguous genomic regions loaded from a BED file.
//!
//! A read whose alignment overlaps any region here cannot be early-assigned.
//! Regions are stored per reference sequence (keyed by ref-id), sorted by
//! start position, and queried via binary search.

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

/// Per-chromosome sorted list of ambiguous regions.
#[derive(Debug, Default)]
pub(crate) struct AmbiguousRegions {
    per_ref: HashMap<usize, Vec<Region>>,
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
        let mut per_ref: HashMap<usize, Vec<Region>> = HashMap::new();

        let mut record = bed::Record::<3>::default();
        loop {
            let n = reader
                .read_record(&mut record)
                .map_err(|e| Error::BedParseError {
                    path: path.to_path_buf(),
                    source: e,
                })?;
            if n == 0 {
                break;
            }
            let chrom = record.reference_sequence_name();
            let chrom_str: String =
                String::try_from(chrom).map_err(|e| Error::BedInvalidChromName(e.to_string()))?;
            let id = match name_to_id.get(chrom_str.as_str()) {
                Some(&id) => id,
                None => continue,
            };
            let start = record
                .feature_start()
                .map_err(|e| Error::BedStartError(e.to_string()))?;
            let end = record
                .feature_end()
                .ok_or(Error::BedRecordMissingEnd)?
                .map_err(|e| Error::BedEndError(e.to_string()))?;

            per_ref.entry(id).or_default().push(Region {
                start: usize::from(start),
                end: usize::from(end),
            });
        }

        for regions in per_ref.values_mut() {
            regions.sort_unstable_by_key(|r| r.start);
        }

        Ok(Self { per_ref })
    }

    /// Returns `true` if `[read_start, read_end)` overlaps any ambiguous region.
    pub(crate) fn overlaps(&self, ref_id: usize, read_start: usize, read_end: usize) -> bool {
        let Some(regions) = self.per_ref.get(&ref_id) else {
            return false;
        };
        let first = regions.partition_point(|r| r.end <= read_start);
        regions[first..].iter().any(|r| r.start < read_end)
    }

    pub(crate) fn is_empty(&self) -> bool {
        self.per_ref.is_empty()
    }

    #[cfg(test)]
    pub(crate) fn from_test(per_ref: HashMap<usize, Vec<Region>>) -> Self {
        Self { per_ref }
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
        for v in per_ref.values_mut() {
            v.sort_unstable_by_key(|r| r.start);
        }
        AmbiguousRegions { per_ref }
    }

    #[test]
    fn test_no_overlap_touching_left_edge() {
        assert!(!store(&[(0, 100, 200)]).overlaps(0, 50, 100));
    }

    #[test]
    fn test_no_overlap_touching_right_edge() {
        assert!(!store(&[(0, 100, 200)]).overlaps(0, 200, 300));
    }

    #[test]
    fn test_overlap_contained() {
        assert!(store(&[(0, 100, 200)]).overlaps(0, 120, 150));
    }

    #[test]
    fn test_overlap_straddles_left() {
        assert!(store(&[(0, 100, 200)]).overlaps(0, 50, 150));
    }

    #[test]
    fn test_overlap_straddles_right() {
        assert!(store(&[(0, 100, 200)]).overlaps(0, 150, 250));
    }

    #[test]
    fn test_wrong_ref_id() {
        assert!(!store(&[(0, 100, 200)]).overlaps(1, 100, 200));
    }

    #[test]
    fn test_multiple_regions() {
        let s = store(&[(0, 10, 20), (0, 100, 200), (0, 300, 400)]);
        assert!(s.overlaps(0, 150, 250));
        assert!(!s.overlaps(0, 200, 300));
        assert!(s.overlaps(0, 15, 25));
    }
}
