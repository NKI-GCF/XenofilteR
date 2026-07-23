//! Scored genomic regions: positive-score BED and strand-aware ambiguous BED.
//!
//! BED format: chrom, start(0-based), end, name, score(0-1000), strand(+/-/.)
//! Columns 4-6 are optional. Missing strand = any. Missing score = 1000.

use crate::Error;
use crate::variant::store::load_lappers;
use rust_lapper::{Interval, Lapper};
use std::collections::HashMap;
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
            b'+' => Self::Plus,
            b'-' => Self::Minus,
            _ => Self::Any,
        }
    }
    /// True when the BED strand is compatible with `read_is_reverse`.
    pub(crate) fn matches(&self, read_is_reverse: bool) -> bool {
        match self {
            Self::Any => true,
            Self::Plus => !read_is_reverse,
            Self::Minus => read_is_reverse,
        }
    }
}

/// How the BED score column is converted to a log-likelihood delta.
/// Variants are intentionally simple; add via PR after empirical evaluation.
#[derive(Clone, Copy, Debug)]
pub(crate) enum ScoreFn {
    /// `bed_score / 1000 * weight`  -- linear scale (default weight = 1.0)
    Linear(f64),
    /// `ln(bed_score + 1) / ln(1001) * weight`
    Log(f64),
    /// Constant `weight` regardless of bed_score
    Constant(f64),
    /// `(overlap_bases / region_len) * (bed_score / 1000) * weight`
    OverlapFraction(f64),
}

impl ScoreFn {
    /// Compute the score delta added to the stream's NW score.
    pub(crate) fn apply(&self, bed_score: f64, overlap_bases: usize, region_len: usize) -> f64 {
        match self {
            Self::Linear(w) => (bed_score / 1000.0) * w,
            Self::Log(w) => (bed_score + 1.0).ln() / 1001.0_f64.ln() * w,
            Self::Constant(w) => *w,
            Self::OverlapFraction(w) => {
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

impl Default for ScoreFn {
    fn default() -> Self {
        Self::Linear(1.0)
    }
}

impl std::str::FromStr for ScoreFn {
    type Err = String;
    fn from_str(s: &str) -> Result<Self, String> {
        // Format: "fn[:weight]", e.g. "linear", "log:2.0", "constant:5.0"
        let (name, weight) = match s.split_once(':') {
            Some((n, w)) => (n, w.parse::<f64>().map_err(|e| e.to_string())?),
            None => (s, 1.0),
        };
        match name {
            "linear" => Ok(Self::Linear(weight)),
            "log" => Ok(Self::Log(weight)),
            "constant" => Ok(Self::Constant(weight)),
            "overlap_fraction" => Ok(Self::OverlapFraction(weight)),
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

impl PartialEq for ScoredRegion {
    fn eq(&self, other: &Self) -> bool {
        self.ref_id == other.ref_id
            && self.start == other.start
            && self.end == other.end
            && self.strand == other.strand
    }
}

impl Eq for ScoredRegion {}

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
    per_chr: HashMap<usize, Lapper<usize, ScoredRegion>>,
}

impl ScoredRegions {
    pub(crate) fn from_bed(
        path: &Path,
        name_to_id: &HashMap<String, usize>,
    ) -> Result<Self, Error> {
        use std::io::{BufRead, BufReader};
        let f = std::fs::File::open(path)?;
        let lines = BufReader::new(f).lines().filter_map(|l| {
            let line = l.ok()?;
            let line = line.trim();
            (!line.is_empty() && !line.starts_with('#') && !line.starts_with("track"))
                .then(|| Ok(line.to_string()))
        });

        let per_chr = load_lappers(
            lines,
            name_to_id,
            |line: &String| line.split('\t').next().unwrap_or("").to_string(),
            |line, ref_id| {
                let cols: Vec<&str> = line.split('\t').collect();
                if cols.len() < 3 {
                    return Ok(None);
                }
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
                let region = ScoredRegion {
                    ref_id,
                    start,
                    end,
                    score,
                    strand,
                };
                Ok(Some(Interval {
                    start,
                    stop: end,
                    val: region,
                }))
            },
        )?;

        Ok(Self { per_chr })
    }

    pub(crate) fn overlapping(
        &self,
        ref_id: usize,
        start: usize,
        end: usize,
        read_is_reverse: bool,
    ) -> impl Iterator<Item = &ScoredRegion> {
        self.per_chr
            .get(&ref_id)
            .into_iter()
            .flat_map(move |lapper| lapper.find(start, end))
            .map(|iv| &iv.val)
            .filter(move |r| r.strand.matches(read_is_reverse))
    }

    pub(crate) fn any_overlap(&self, ref_id: usize, s: usize, e: usize, rev: bool) -> bool {
        self.overlapping(ref_id, s, e, rev).next().is_some()
    }

    pub(crate) fn overlaps(&self, ref_id: usize, s: usize, e: usize) -> bool {
        self.per_chr
            .get(&ref_id)
            .is_some_and(|lapper| lapper.find(s, e).next().is_some())
    }

    pub(crate) fn is_empty(&self) -> bool {
        self.per_chr.values().all(|lapper| lapper.is_empty())
    }

    pub(crate) fn empty() -> Self {
        Self {
            per_chr: HashMap::new(),
        }
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
    #[test]
    fn score_fn_from_str_table() {
        use std::str::FromStr;
        assert!(
            matches!(ScoreFn::from_str("linear"), Ok(ScoreFn::Linear(w)) if (w - 1.0).abs() < 1e-9)
        );
        assert!(
            matches!(ScoreFn::from_str("log:2.0"), Ok(ScoreFn::Log(w)) if (w - 2.0).abs() < 1e-9)
        );
        assert!(
            matches!(ScoreFn::from_str("constant:5"), Ok(ScoreFn::Constant(w)) if (w - 5.0).abs() < 1e-9)
        );
        assert!(matches!(
            ScoreFn::from_str("overlap_fraction:0.5"),
            Ok(ScoreFn::OverlapFraction(w)) if (w - 0.5).abs() < 1e-9
        ));
        assert!(ScoreFn::from_str("bogus").is_err());
        assert!(ScoreFn::from_str("linear:notanumber").is_err());
        assert!(ScoreFn::from_str("linear:").is_err());
    }
}
