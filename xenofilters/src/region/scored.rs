//! Scored genomic regions: positive-score BED and strand-aware ambiguous BED.
//!
//! BED format: chrom, start(0-based), end, name, score(0-1000), strand(+/-/.)
//! Columns 4-6 are optional. Missing strand = any. Missing score = 1000.

use crate::Error;
use std::path::Path;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum Strand {
    Plus,
    Minus,
    Any,
}

impl Strand {
    fn from_byte(b: u8) -> Self {
        match b {
            b'+' => Strand::Plus,
            b'-' => Strand::Minus,
            _ => Strand::Any,
        }
    }
    /// True when the BED strand is compatible with `read_is_reverse`.
    pub(crate) fn matches(&self, read_is_reverse: bool) -> bool {
        match self {
            Strand::Any => true,
            Strand::Plus => !read_is_reverse,
            Strand::Minus => read_is_reverse,
        }
    }
}

/// How the BED score column is converted to a log-likelihood delta.
/// Variants are intentionally simple; add via PR after empirical evaluation.
#[derive(Clone, Copy, Debug)]
pub(crate) enum ScoreFn {
    /// `bed_score / 1000 × weight`  — linear scale (default weight = 1.0)
    Linear(f64),
    /// `ln(bed_score + 1) / ln(1001) × weight`
    Log(f64),
    /// Constant `weight` regardless of bed_score
    Constant(f64),
    /// `(overlap_bases / region_len) × (bed_score / 1000) × weight`
    OverlapFraction(f64),
}

impl Default for ScoreFn {
    fn default() -> Self {
        ScoreFn::Linear(1.0)
    }
}

impl ScoreFn {
    /// Compute the score delta added to the stream's NW score.
    pub(crate) fn apply(&self, bed_score: f64, overlap_bases: usize, region_len: usize) -> f64 {
        match self {
            ScoreFn::Linear(w) => (bed_score / 1000.0) * w,
            ScoreFn::Log(w) => (bed_score + 1.0).ln() / (1001.0_f64).ln() * w,
            ScoreFn::Constant(w) => *w,
            ScoreFn::OverlapFraction(w) => {
                let frac = if region_len == 0 {
                    0.0
                } else {
                    overlap_bases as f64 / region_len as f64
                };
                frac * (bed_score / 1000.0) * w
            }
        }
    }
}

impl std::str::FromStr for ScoreFn {
    type Err = String;
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        // Format: "fn[:weight]", e.g. "linear", "log:2.0", "constant:5.0"
        let (name, weight) = match s.split_once(':') {
            Some((n, w)) => (n, w.parse::<f64>().map_err(|e| e.to_string())?),
            None => (s, 1.0f64),
        };
        match name {
            "linear" => Ok(ScoreFn::Linear(weight)),
            "log" => Ok(ScoreFn::Log(weight)),
            "constant" => Ok(ScoreFn::Constant(weight)),
            "overlap_fraction" => Ok(ScoreFn::OverlapFraction(weight)),
            other => Err(format!("unknown score fn '{other}'")),
        }
    }
}

#[derive(Debug, Clone)]
pub(crate) struct ScoredRegion {
    pub(crate) ref_id: usize, // 0-based reference sequence index
    pub(crate) start: usize,  // 0-based
    pub(crate) end: usize,    // exclusive
    pub(crate) score: f64,    // BED score column (0-1000)
    pub(crate) strand: Strand,
}

impl ScoredRegion {
    pub(crate) fn len(&self) -> usize {
        self.end.saturating_sub(self.start)
    }

    pub(crate) fn overlap_bases(&self, read_start: usize, read_end: usize) -> usize {
        let s = self.start.max(read_start);
        let e = self.end.min(read_end);
        e.saturating_sub(s)
    }
}

/// In-memory collection of strand-aware, scored BED regions.
/// Sorted per ref_id by start for binary-search lookups.
pub(crate) struct ScoredRegions {
    /// per_ref[ref_id] = sorted Vec of ScoredRegion
    per_ref: Vec<Vec<ScoredRegion>>,
}

impl ScoredRegions {
    pub(crate) fn empty() -> Self {
        Self {
            per_ref: Vec::new(),
        }
    }

    pub(crate) fn from_bed(
        path: &Path,
        name_to_id: &std::collections::HashMap<String, usize>,
    ) -> Result<Self, Error> {
        use std::io::{BufRead, BufReader};
        let f = std::fs::File::open(path)?;
        let mut per_ref: Vec<Vec<ScoredRegion>> = Vec::new();
        for line in BufReader::new(f).lines() {
            let line = line?;
            let line = line.trim();
            if line.is_empty() || line.starts_with('#') || line.starts_with("track") {
                continue;
            }
            let cols: Vec<&str> = line.split('\t').collect();
            if cols.len() < 3 {
                continue;
            }
            let chrom = cols[0];
            let start: usize = cols[1].parse()?;
            let end: usize = cols[2].parse()?;
            let score = cols
                .get(4)
                .and_then(|s| s.parse::<f64>().ok())
                .unwrap_or(1000.0);
            let strand = cols
                .get(5)
                .and_then(|s| s.as_bytes().first().copied())
                .map(Strand::from_byte)
                .unwrap_or(Strand::Any);
            let ref_id = match name_to_id.get(chrom) {
                Some(&id) => id,
                None => {
                    tracing::debug!("BED chrom '{chrom}' not in reference; skipped");
                    continue;
                }
            };
            if per_ref.len() <= ref_id {
                per_ref.resize_with(ref_id + 1, Vec::new);
            }
            per_ref[ref_id].push(ScoredRegion {
                ref_id,
                start,
                end,
                score,
                strand,
            });
        }
        for v in &mut per_ref {
            v.sort_unstable_by_key(|r| r.start);
        }
        Ok(Self { per_ref })
    }
    /// Strand-agnostic overlap check — used by call sites that don't have
    /// read orientation available (e.g. HashLookup's raw MappedRecord path,
    /// which does not currently track is_reverse per record in this context).
    /// Equivalent to `any_overlap(ref_id, start, end, /* any strand */)`.
    pub(crate) fn overlaps(&self, ref_id: usize, start: usize, end: usize) -> bool {
        let regions = match self.per_ref.get(ref_id) {
            Some(v) => v.as_slice(),
            None => return false,
        };
        let lo = regions.partition_point(|r| r.end <= start);
        regions[lo..]
            .iter()
            .take_while(|r| r.start < end)
            .next()
            .is_some()
    }

    pub(crate) fn is_empty(&self) -> bool {
        self.per_ref.is_empty()
    }

    /// Returns all regions overlapping [read_start, read_end) on ref_id
    /// whose strand matches read_is_reverse.
    pub(crate) fn overlapping(
        &self,
        ref_id: usize,
        read_start: usize,
        read_end: usize,
        read_is_reverse: bool,
    ) -> impl Iterator<Item = &ScoredRegion> {
        let regions = match self.per_ref.get(ref_id) {
            Some(v) => v.as_slice(),
            None => &[],
        };
        // Binary search for first region whose end > read_start.
        let lo = regions.partition_point(|r| r.end <= read_start);
        regions[lo..]
            .iter()
            .take_while(move |r| r.start < read_end)
            .filter(move |r| r.strand.matches(read_is_reverse))
    }

    /// True if any strand-compatible region overlaps [read_start, read_end).
    /// Used by ambiguous-region forcing.
    pub(crate) fn any_overlap(
        &self,
        ref_id: usize,
        read_start: usize,
        read_end: usize,
        read_is_reverse: bool,
    ) -> bool {
        self.overlapping(ref_id, read_start, read_end, read_is_reverse)
            .next()
            .is_some()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn strand_aware_overlap_table() {
        struct Row {
            read_rev: bool,
            strand: Strand,
            want: bool,
        }
        let region = ScoredRegion {
            ref_id: 0,
            start: 100,
            end: 200,
            score: 500.0,
            strand: Strand::Plus,
        };
        let cases = [
            Row {
                read_rev: false,
                strand: Strand::Plus,
                want: true,
            },
            Row {
                read_rev: true,
                strand: Strand::Plus,
                want: false,
            },
            Row {
                read_rev: false,
                strand: Strand::Minus,
                want: false,
            },
            Row {
                read_rev: true,
                strand: Strand::Minus,
                want: true,
            },
            Row {
                read_rev: false,
                strand: Strand::Any,
                want: true,
            },
            Row {
                read_rev: true,
                strand: Strand::Any,
                want: true,
            },
        ];
        for c in &cases {
            let r = ScoredRegion {
                strand: c.strand,
                ..region.clone()
            };
            assert_eq!(r.strand.matches(c.read_rev), c.want);
        }
    }

    #[test]
    fn score_fn_table() {
        struct Row {
            fn_: ScoreFn,
            bed: f64,
            ob: usize,
            rlen: usize,
            lo: f64,
            hi: f64,
        }
        let cases = [
            Row {
                fn_: ScoreFn::Linear(1.0),
                bed: 1000.,
                ob: 50,
                rlen: 100,
                lo: 0.999,
                hi: 1.001,
            },
            Row {
                fn_: ScoreFn::Constant(5.0),
                bed: 0.,
                ob: 0,
                rlen: 1,
                lo: 4.999,
                hi: 5.001,
            },
            Row {
                fn_: ScoreFn::Log(1.0),
                bed: 1000.,
                ob: 0,
                rlen: 1,
                lo: 0.999,
                hi: 1.001,
            },
            Row {
                fn_: ScoreFn::OverlapFraction(10.0),
                bed: 1000.,
                ob: 50,
                rlen: 100,
                lo: 4.999,
                hi: 5.001,
            },
        ];
        for c in &cases {
            let got = c.fn_.apply(c.bed, c.ob, c.rlen);
            assert!(
                got >= c.lo && got <= c.hi,
                "{got} not in [{}, {}]",
                c.lo,
                c.hi
            );
        }
    }
}
