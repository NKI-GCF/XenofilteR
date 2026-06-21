//! Tabix-backed random-access queries for the Collated algorithm.
//!
//! Because the Collated algorithm processes records by name order (not
//! position order), BED and VCF files cannot be walked with a cursor.
//! Instead they are queried by genomic coordinate using a tabix index.
//!
//! noodles 0.111.0 API notes:
//! - `tabix::io::Reader` wraps a `bgzf::io::Reader` and provides `query()`.
//! - Index is read via `tabix::io::read(path)`.
//! - VCF indexed reader: `vcf::io::IndexedReader::new(bgzf_reader, index)`.
//! - BED does not have a noodles IndexedReader; we use the tabix reader
//!   directly and parse BED lines as raw text.

use anyhow::{anyhow, Result};
use noodles::bgzf;
use noodles::core::Region;
use noodles::{tabix, vcf};
use std::fs::File;
use std::path::Path;

// ---------------------------------------------------------------------------
// TabixBed
// ---------------------------------------------------------------------------

/// Tabix-indexed BED file for random-access ambiguous-region queries.
/// We use the tabix reader directly since noodles-bed lacks an IndexedReader.
pub(crate) struct TabixBed {
    reader: bgzf::io::Reader<File>,
    index: tabix::Index,
}

impl TabixBed {
    pub(crate) fn open(path: &Path) -> Result<Self> {
        let tbi_path = path.with_extension(
            path.extension()
                .map(|e| format!("{}.tbi", e.to_string_lossy()))
                .unwrap_or_else(|| "tbi".into()),
        );
        let index = tabix::io::read(&tbi_path)
            .map_err(|e| anyhow!("Cannot read tabix index {}: {e}", tbi_path.display()))?;
        let reader = bgzf::io::Reader::new(
            File::open(path).map_err(|e| anyhow!("Cannot open BED {}: {e}", path.display()))?,
        );
        Ok(Self { reader, index })
    }

    /// Returns `true` if any BED record overlaps `[start, end)` (0-based).
    pub(crate) fn overlaps(&mut self, chrom: &str, start: usize, end: usize) -> Result<bool> {
        // Region uses 1-based inclusive coordinates.
        let region: Region = format!("{}:{}-{}", chrom, start + 1, end)
            .parse()
            .map_err(|e| anyhow!("Invalid region: {e}"))?;
        let chunks = self
            .index
            .query(&region)
            .map_err(|e| anyhow!("Tabix BED query: {e}"))?;
        // Any chunk means at least one record overlaps.
        Ok(!chunks.is_empty())
    }
}

// ---------------------------------------------------------------------------
// TabixVcf
// ---------------------------------------------------------------------------

/// Tabix-indexed VCF file for random-access diagnostic-variant queries.
pub(crate) struct TabixVcf {
    reader: vcf::io::IndexedReader<bgzf::io::Reader<File>>,
    header: vcf::Header,
}

impl TabixVcf {
    pub(crate) fn open(path: &Path) -> Result<Self> {
        let tbi_path = path.with_extension(
            path.extension()
                .map(|e| format!("{}.tbi", e.to_string_lossy()))
                .unwrap_or_else(|| "tbi".into()),
        );
        let index = tabix::io::read(&tbi_path)
            .map_err(|e| anyhow!("Cannot read tabix index {}: {e}", tbi_path.display()))?;
        let bgzf = bgzf::io::Reader::new(
            File::open(path).map_err(|e| anyhow!("Cannot open VCF {}: {e}", path.display()))?,
        );
        let mut reader = vcf::io::IndexedReader::new(bgzf, index);
        let header = reader
            .read_header()
            .map_err(|e| anyhow!("VCF header: {e}"))?;
        Ok(Self { reader, header })
    }

    /// Returns `true` if any diagnostic variant overlaps `[start, end)` (1-based).
    pub(crate) fn overlaps(&mut self, chrom: &str, start: usize, end: usize) -> Result<bool> {
        let region: Region = format!("{}:{}-{}", chrom, start, end)
            .parse()
            .map_err(|e| anyhow!("Invalid region: {e}"))?;
        let mut query = self
            .reader
            .query(&self.header, &region)
            .map_err(|e| anyhow!("Tabix VCF query: {e}"))?;
        // query() returns an iterator in noodles 0.111.0.
        Ok(query.next().is_some())
    }
}
