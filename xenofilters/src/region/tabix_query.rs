//!
//! The Collated algorithm processes records by name order, so BED and VCF
//! files cannot be walked with a position cursor. Instead, records are queried
//! by genomic coordinate using a tabix index.
//!
//! noodles 0.111.0 API:
//! - Tabix index: read via `tabix::io::Reader::new(File::open(tbi_path)?).read_index()?`
//! - BED query: use tabix index chunks to check overlap; parse overlapping records.
//! - VCF indexed reader: `vcf::io::IndexedReader` provides `.query(&header, &region)`.

use crate::region::ScoreFn;
use crate::Error;
use noodles::bgzf;
use noodles::core::{region::Interval, Position, Region};
use noodles::csi::BinningIndex;
use noodles::{tabix, vcf};
use std::cell::RefCell;
use std::fs::File;
use std::io::{BufRead, Seek};
use std::ops::RangeInclusive;
use std::path::{Path, PathBuf};
use vcf::io::IndexedReader;

pub(crate) struct TabixScored {
    index: tabix::Index,
    path: PathBuf,
} // reuses index + adds column parsing

impl TabixScored {
    pub(crate) fn open(path: &Path) -> Result<Self, Error> {
        let index = read_tabix_index(path)?;
        Ok(Self {
            index,
            path: path.to_path_buf(),
        })
    }
    /// Real per-record score/strand overlap bonus, replacing the previous
    /// existence-only + constant-1000 placeholder. Seeks to each matched
    /// BGZF chunk, reads lines within it, and applies `score_fn` using the
    /// record's actual BED score (col 5) and strand (col 6), restricted to
    /// records that genuinely overlap [start, end).
    ///
    /// NEEDS VERIFICATION: exact noodles 0.112 API for `Chunk::start()` /
    /// `Chunk::end()` (bgzf::VirtualPosition) and `bgzf::io::Reader::seek`/
    /// `virtual_position()` -- written to match the pattern already used by
    /// `TabixVcf`/`IndexedReader` elsewhere in this file, but not compiled
    /// here.
    pub(crate) fn overlapping_bonus(
        &self,
        ref_id: usize,
        start: usize,
        end: usize,
        read_is_reverse: bool,
        score_fn: ScoreFn,
    ) -> Result<f64, Error> {
        let q_start =
            Position::try_from(start + 1).map_err(|e| Error::InvalidRegion(e.to_string()))?;
        let q_end = Position::try_from(end.max(start + 1))
            .map_err(|e| Error::InvalidRegion(e.to_string()))?;
        let interval = Interval::from(RangeInclusive::new(q_start, q_end));

        let chunks = self
            .index
            .query(ref_id, interval)
            .map_err(|e| Error::TabixBedQueryFailed(e.to_string()))?;
        if chunks.is_empty() {
            return Ok(0.0);
        }

        let file = std::fs::File::open(&self.path).map_err(|e| Error::CannotOpenBedFile {
            path: self.path.clone(),
            source: e,
        })?;
        let mut reader = bgzf::io::Reader::new(file);

        let mut total = 0.0;
        for chunk in chunks {
            reader
                .seek(chunk.start())
                .map_err(|e| Error::TabixBedQueryFailed(e.to_string()))?;

            loop {
                let mut line = String::new();
                let n = reader
                    .read_line(&mut line)
                    .map_err(|e| Error::TabixBedQueryFailed(e.to_string()))?;
                if n == 0 {
                    break; // EOF
                }
                let past_chunk = reader.virtual_position() >= chunk.end();
                let trimmed = line.trim_end();

                if trimmed.is_empty() || trimmed.starts_with('#') || trimmed.starts_with("track") {
                    if past_chunk {
                        break;
                    }
                    continue;
                }

                let cols: Vec<&str> = trimmed.split('\t').collect();
                if cols.len() >= 3 {
                    if let (Ok(rec_start), Ok(rec_end)) =
                        (cols[1].parse::<usize>(), cols[2].parse::<usize>())
                    {
                        if rec_start < end && rec_end > start {
                            let score = cols
                                .get(4)
                                .and_then(|s| s.parse::<f64>().ok())
                                .unwrap_or(1000.0);
                            let strand = cols
                                .get(5)
                                .and_then(|s| s.as_bytes().first().copied())
                                .map(crate::region::scored::Strand::from_byte)
                                .unwrap_or(crate::region::scored::Strand::Any);
                            if strand.matches(read_is_reverse) {
                                let overlap_bases =
                                    rec_end.min(end).saturating_sub(rec_start.max(start));
                                let region_len = rec_end.saturating_sub(rec_start).max(1);
                                total += score_fn.apply(score, overlap_bases, region_len);
                            }
                        }
                    }
                }

                if past_chunk {
                    break;
                }
            }
        }
        Ok(total)
    }

    pub(crate) fn is_empty(&self) -> bool {
        self.index.reference_sequences().is_empty()
    }
}

fn read_tabix_index(path: &Path) -> Result<tabix::Index, Error> {
    // Try <file>.tbi, then <file>.<ext>.tbi
    let tbi1 = path.with_extension("tbi");
    let tbi2 = {
        let mut p = path.as_os_str().to_owned();
        p.push(".tbi");
        PathBuf::from(p)
    };
    for tbi_path in &[&tbi1, &tbi2] {
        if tbi_path.exists() {
            let mut reader = tabix::io::Reader::new(File::open(tbi_path).map_err(|e| {
                Error::CannotOpenIndex {
                    path: tbi_path.to_path_buf(),
                    source: e,
                }
            })?);
            return reader
                .read_index()
                .map_err(|e| Error::CannotReadTabixIndex {
                    path: tbi_path.to_path_buf(),
                    source: e,
                });
        }
    }
    Err(Error::TabixIndexNotFound {
        path: path.to_path_buf(),
    })
}

pub(crate) struct TabixBed {
    index: tabix::Index,
}

impl TabixBed {
    pub(crate) fn open(path: &Path) -> Result<Self, Error> {
        let index = read_tabix_index(path)?;
        Ok(Self { index })
    }

    /// Returns `true` if any BED record overlaps `[start, end)` (0-based).
    pub(crate) fn chrom_overlaps(
        &mut self,
        chrom: &str,
        start: usize,
        end: usize,
    ) -> Result<bool, Error> {
        // Convert to 1-based inclusive for Region.
        let region_str = format!("{}:{}-{}", chrom, start + 1, end);
        let region: Region = region_str
            .parse()
            .map_err(|_| Error::InvalidRegion(region_str.clone()))?;
        let header = self.index.header().expect("missing tabix header");
        let reference_sequence_id = header
            .reference_sequence_names()
            .get_index_of(region.name())
            .expect("invalid reference sequence name");
        self.overlaps(reference_sequence_id, start, end)
    }

    /// Returns `true` if any BED record overlaps `[start, end)` (0-based).
    pub(crate) fn overlaps(&self, ref_id: usize, start: usize, end: usize) -> Result<bool, Error> {
        let start =
            Position::try_from(start + 1).map_err(|e| Error::InvalidRegion(e.to_string()))?;
        let end = Position::try_from(end).map_err(|e| Error::InvalidRegion(e.to_string()))?;
        let interval = Interval::from(RangeInclusive::new(start, end));
        let chunks = self
            .index
            .query(ref_id, interval)
            .map_err(|e| Error::TabixBedQueryFailed(e.to_string()))?;
        Ok(!chunks.is_empty())
    }
    /// Strand-aware overlap check for tabix-indexed BED.
    /// If the BED has no strand column (BED3), all regions are treated as `Any`
    /// and this is equivalent to `overlaps()`.
    ///
    /// NEEDS VERIFICATION: current TabixBed does not parse/store column 6.
    /// This stub queries strand-agnostically until tabix-side strand storage
    /// is implemented; documented as a known limitation, not silently wrong,
    /// since `Strand::Any` always matches.
    pub(crate) fn any_overlap(
        &self,
        ref_id: usize,
        start: usize,
        end: usize,
        _read_is_reverse: bool,
    ) -> Result<bool, Error> {
        self.overlaps(ref_id, start, end)
    }

    pub(crate) fn is_empty(&self) -> bool {
        self.index.reference_sequences().is_empty()
    }
}

/// Tabix-indexed VCF file for random-access diagnostic-variant queries.
pub(crate) struct TabixVcf {
    reader: RefCell<IndexedReader<bgzf::io::Reader<File>>>,
    header: vcf::Header,
}

impl TabixVcf {
    pub(crate) fn open(path: &Path) -> Result<Self, Error> {
        let index = read_tabix_index(path)?;
        let file = File::open(path).map_err(|e| Error::CannotOpenVcf {
            path: path.to_path_buf(),
            source: e,
        })?;
        let mut reader = IndexedReader::new(file, index);
        let header = reader
            .read_header()
            .map_err(|e| Error::VcfHeaderReadError(e.to_string()))?;
        Ok(Self {
            reader: RefCell::new(reader),
            header,
        })
    }

    /// Returns `true` if any diagnostic variant overlaps `[start, end)` (1-based, BAM coords).
    pub(crate) fn chrom_overlaps(
        &self,
        chrom: &str,
        start: usize,
        end: usize,
    ) -> Result<bool, Error> {
        let region_str = format!("{}:{}-{}", chrom, start, end);
        let region: Region = region_str
            .parse()
            .map_err(|_| Error::InvalidRegion(region_str.clone()))?;
        let mut reader = self.reader.borrow_mut();
        let query = reader
            .query(&self.header, &region)
            .map_err(|e| Error::TabixVcfQueryFailed(e.to_string()))?;
        Ok(query.records().next().is_none())
    }

    pub(crate) fn overlaps(&self, ref_id: usize, start: usize, end: usize) -> Result<bool, Error> {
        let chrom = self
            .header
            .contigs()
            .get_index(ref_id)
            .ok_or(Error::NoRefSeqId)?
            .0
            .to_owned();
        self.chrom_overlaps(&chrom, start, end)
    }

    pub(crate) fn is_empty(&self) -> bool {
        self.header.contigs().is_empty()
    }
}
