//! Tabix-backed random-access queries for the Collated algorithm.
//!
//! Because the Collated algorithm processes records by name order (not
//! position order), BED and VCF files cannot be walked with a cursor.
//! Instead they are queried by genomic coordinate using a tabix index.
//!
//! This module wraps noodles tabix+bgzf to provide cheap overlap queries
//! against sorted, bgzipped, tabix-indexed BED and VCF files.

use anyhow::{anyhow, Result};
use noodles::bgzf;
use noodles::core::Region;
use noodles::tabix;
use noodles::{bed, bcf, vcf};
use std::fs::File;
use std::path::Path;

/// Tabix-indexed BED file for random-access ambiguous-region queries.
pub(crate) struct TabixBed {
    reader: bed::io::IndexedReader<bgzf::Reader<File>>,
}

impl TabixBed {
    pub(crate) fn open(path: &Path) -> Result<Self> {
        let reader = bed::io::IndexedReader::new(
            bgzf::Reader::new(
                File::open(path)
                    .map_err(|e| anyhow!("Cannot open BED {}: {}", path.display(), e))?,
            ),
            tabix::read(path.with_extension("bed.gz.tbi"))
                .or_else(|_| tabix::read(path.with_extension("tbi")))
                .map_err(|e| anyhow!("Cannot read tabix index for {}: {}", path.display(), e))?,
        );
        Ok(Self { reader })
    }

    /// Returns `true` if any BED record overlaps the given region.
    pub(crate) fn overlaps(&mut self, chrom: &str, start: usize, end: usize) -> Result<bool> {
        let region: Region = format!("{}:{}-{}", chrom, start, end)
            .parse()
            .map_err(|e| anyhow!("Invalid region {chrom}:{start}-{end}: {e}"))?;
        let query = self
            .reader
            .query(&region)
            .map_err(|e| anyhow!("Tabix BED query failed: {e}"))?;
        for result in query {
            let _record: bed::Record<3> = result
                .map_err(|e| anyhow!("BED record parse error: {e}"))?;
            return Ok(true); // any overlap suffices
        }
        Ok(false)
    }
}

/// Tabix-indexed VCF/BCF file for random-access diagnostic-variant queries.
pub(crate) struct TabixVcf {
    reader: vcf::io::IndexedReader<bgzf::Reader<File>>,
    header: vcf::Header,
}

impl TabixVcf {
    pub(crate) fn open(path: &Path) -> Result<Self> {
        let index = tabix::read(path.with_extension("vcf.gz.tbi"))
            .or_else(|_| tabix::read(path.with_extension("tbi")))
            .map_err(|e| anyhow!("Cannot read tabix index for {}: {}", path.display(), e))?;
        let mut reader = vcf::io::IndexedReader::new(
            bgzf::Reader::new(
                File::open(path)
                    .map_err(|e| anyhow!("Cannot open VCF {}: {}", path.display(), e))?,
            ),
            index,
        );
        let header = reader
            .read_header()
            .map_err(|e| anyhow!("VCF header read error: {e}"))?;
        Ok(Self { reader, header })
    }

    /// Returns `true` if any diagnostic variant overlaps `[start, end)`.
    pub(crate) fn overlaps(&mut self, chrom: &str, start: usize, end: usize) -> Result<bool> {
        let region: Region = format!("{}:{}-{}", chrom, start, end)
            .parse()
            .map_err(|e| anyhow!("Invalid region: {e}"))?;
        let query = self
            .reader
            .query(&self.header, &region)
            .map_err(|e| anyhow!("Tabix VCF query failed: {e}"))?;
        for result in query {
            let _record = result.map_err(|e| anyhow!("VCF record error: {e}"))?;
            return Ok(true);
        }
        Ok(false)
    }
}
