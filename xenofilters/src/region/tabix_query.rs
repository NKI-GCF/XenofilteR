//!
//! The Collated algorithm processes records by name order, so BED and VCF
//! files cannot be walked with a position cursor. Instead, records are queried
//! by genomic coordinate using a tabix index.
//!
//! noodles 0.111.0 API:
//! - Tabix index: read via `tabix::io::Reader::new(File::open(tbi_path)?).read_index()?`
//! - BED query: use tabix index chunks to check overlap; parse overlapping records.
//! - VCF indexed reader: `vcf::io::IndexedReader` provides `.query(&header, &region)`.

use noodles::bgzf;
use noodles::core::Region;
use noodles::csi::BinningIndex;
use noodles::{tabix, vcf};
use std::fs::File;
use std::path::{Path, PathBuf};
use crate::Error;

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
            let mut reader = tabix::io::Reader::new(
                File::open(tbi_path)
                    .map_err(|e| Error::CannotOpenIndex {
                        path: tbi_path.to_path_buf(),
                        source: e
                    })?,
            );
            return reader
                .read_index()
                .map_err(|e| Error::CannotReadTabixIndex {
                    path: tbi_path.to_path_buf(),
                    source: e
                });
        }
    }
    Err(Error::TabixIndexNotFound { path: path.to_path_buf() })
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
    pub(crate) fn overlaps(&mut self, chrom: &str, start: usize, end: usize) -> Result<bool, Error> {
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
        let chunks = self
            .index
            .query(reference_sequence_id, region.interval())
            .map_err(|e| Error::TabixBedQueryFailed(e.to_string()))?;
        Ok(!chunks.is_empty())
    }
}

/// Tabix-indexed VCF file for random-access diagnostic-variant queries.
pub(crate) struct TabixVcf {
    reader: vcf::io::IndexedReader<bgzf::io::Reader<File>>,
    header: vcf::Header,
}

impl TabixVcf {
    pub(crate) fn open(path: &Path) -> Result<Self, Error> {
        let index = read_tabix_index(path)?;
        let file = File::open(path)
            .map_err(|e| Error::CannotOpenVcf {
                path: path.to_path_buf(),
                source: e
            })?;
        let mut reader = vcf::io::IndexedReader::new(file, index);
        let header = reader
            .read_header()
            .map_err(|e| Error::VcfHeaderReadError(e.to_string()))?;
        Ok(Self { reader, header })
    }

    /// Returns `true` if any diagnostic variant overlaps `[start, end)` (1-based, BAM coords).
    pub(crate) fn overlaps(&mut self, chrom: &str, start: usize, end: usize) -> Result<bool, Error> {
        let region_str = format!("{}:{}-{}", chrom, start, end);
        let region: Region = region_str
            .parse()
            .map_err(|_| Error::InvalidRegion(region_str.clone()))?;
        let query = self
            .reader
            .query(&self.header, &region)
            .map_err(|e| Error::TabixVcfQueryFailed(e.to_string()))?;
        Ok(query.records().next().is_none())
    }
}

