//! Tabix-backed random-access queries for the Collated algorithm.
//!
//! The Collated algorithm processes records by name order, so BED and VCF
//! files cannot be walked with a position cursor. Instead, records are queried
//! by genomic coordinate using a tabix index.
//!
//! noodles 0.111.0 API:
//! - Tabix index: read via `tabix::io::Reader::new(File::open(tbi_path)?).read_index()?`
//! - BED query: use tabix index chunks to check overlap; parse overlapping records.
//! - VCF indexed reader: `vcf::io::IndexedReader` provides `.query(&header, &region)`.

use anyhow::{anyhow, Result};
use noodles::bgzf;
use noodles::core::Region;
use noodles::csi::BinningIndex;
use noodles::{tabix, vcf};
use std::fs::File;
use std::path::Path;

fn read_tabix_index(path: &Path) -> Result<tabix::Index> {
    // Try <file>.tbi, then <file>.<ext>.tbi
    let tbi1 = path.with_extension("tbi");
    let tbi2 = {
        let mut p = path.as_os_str().to_owned();
        p.push(".tbi");
        std::path::PathBuf::from(p)
    };
    for tbi_path in &[&tbi1, &tbi2] {
        if tbi_path.exists() {
            let mut reader = tabix::io::Reader::new(
                File::open(tbi_path)
                    .map_err(|e| anyhow!("Cannot open index {}: {e}", tbi_path.display()))?,
            );
            return reader
                .read_index()
                .map_err(|e| anyhow!("Cannot read tabix index {}: {e}", tbi_path.display()));
        }
    }
    Err(anyhow!(
        "No tabix index found for {} (tried .tbi and .<ext>.tbi)",
        path.display()
    ))
}

// ---------------------------------------------------------------------------
// TabixBed
// ---------------------------------------------------------------------------

/// Tabix-indexed BED file for random-access ambiguous-region queries.
pub(crate) struct TabixBed {
    reader: bgzf::io::Reader<File>,
    index: tabix::Index,
    /// Reference sequence names from the tabix index header.
    ref_names: Vec<Vec<u8>>,
}

impl TabixBed {
    pub(crate) fn open(path: &Path) -> Result<Self> {
        let index = read_tabix_index(path)?;
        let ref_names: Vec<Vec<u8>> = index
            .header()
            .map(|h| {
                h.reference_sequence_names()
                    .iter()
                    .map(|n| n.to_vec())
                    .collect()
            })
            .unwrap_or_default();
        let reader = bgzf::io::Reader::new(
            File::open(path).map_err(|e| anyhow!("Cannot open BED {}: {e}", path.display()))?,
        );
        Ok(Self {
            reader,
            index,
            ref_names,
        })
    }

    /// Returns `true` if any BED record overlaps `[start, end)` (0-based).
    pub(crate) fn overlaps(&mut self, chrom: &str, start: usize, end: usize) -> Result<bool> {
        // Convert to 1-based inclusive for Region.
        let region_str = format!("{}:{}-{}", chrom, start + 1, end);
        let region: Region = region_str
            .parse()
            .map_err(|e| anyhow!("Invalid region {region_str}: {e}"))?;
        let chunks = self
            .index
            .query(region.name(), region.interval())
            .map_err(|e| anyhow!("Tabix BED query failed: {e}"))?;
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
        let index = read_tabix_index(path)?;
        let bgzf_reader = bgzf::io::Reader::new(
            File::open(path).map_err(|e| anyhow!("Cannot open VCF {}: {e}", path.display()))?,
        );
        let mut reader = vcf::io::IndexedReader::new(bgzf_reader, index);
        let header = reader
            .read_header()
            .map_err(|e| anyhow!("VCF header read error: {e}"))?;
        Ok(Self { reader, header })
    }

    /// Returns `true` if any diagnostic variant overlaps `[start, end)` (1-based, BAM coords).
    pub(crate) fn overlaps(&mut self, chrom: &str, start: usize, end: usize) -> Result<bool> {
        let region_str = format!("{}:{}-{}", chrom, start, end);
        let region: Region = region_str
            .parse()
            .map_err(|e| anyhow!("Invalid region {region_str}: {e}"))?;
        let mut query = self
            .reader
            .query(&self.header, &region)
            .map_err(|e| anyhow!("Tabix VCF query failed: {e}"))?;
        Ok(query.next().is_some())
    }
}
