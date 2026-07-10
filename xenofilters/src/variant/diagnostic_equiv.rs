// src/variant/diagnostic_equiv.rs
//
// DiagnosticVariants indel equivalence expansion.
//
// DiagnosticSite is used purely for OVERLAP DETECTION (does the read touch a
// diagnostic position?), not for per-base NW scoring.  Expanding equivalents
// here means the overlap check fires whether the aligner placed the indel at
// the left-normalised or any right-shifted position.
//
// Because DiagnosticSite carries no scoring parameters (only pos and ref_len),
// with_alleles is trivial: just update pos and ref_len = ref_a.len().

use std::{
    collections::HashMap,
    io::{BufRead, Seek},
    path::Path,
};

use anyhow::{Context, Result};
use noodles::{bgzf, bcf, vcf};
use tracing::{debug, warn};

use crate::{
    region::diagnostic::{DiagnosticSite, DiagnosticVariants},
    variant::indel_equiv::{
        classify, enumerate_equivalents, IndelEquivalenceExpander, IndelKind,
    },
};

// ---------------------------------------------------------------------------
// WithAlleles for DiagnosticSite
// ---------------------------------------------------------------------------

impl crate::variant::indel_equiv::WithAlleles for DiagnosticSite {
    fn with_alleles(&self, pos_0based: usize, ref_a: &[u8], _alt_a: &[u8]) -> Self {
        DiagnosticSite {
            pos:     pos_0based,
            ref_len: ref_a.len(),
        }
    }
}

// ---------------------------------------------------------------------------
// AmbiguousRegions equivalence expansion
// ---------------------------------------------------------------------------
//
// AmbiguousRegions (BED-based) store intervals as (start, end).
// Indel-shifted intervals expand identically: the start slides by ±1 per
// shift and the end follows.  Since BED regions are for MASKING (any overlap
// forces full NW), we only need to expand the START positions.
//
// The expansion for BED is simpler than for VCF:
//   - For each BED interval (start, end, ref_id, score, strand):
//     - left-slide start by up to MAX_SHIFT if the context permits
//     - right-slide start by up to MAX_SHIFT
//     - Insert all expanded intervals into the AmbiguousRegions store
//
// This ensures that a read whose CIGAR places an indel ±k bases from the BED
// boundary still triggers full scoring.  Expansion is conservative: overlapping
// expanded intervals are merged by the existing BED merge step.
//
// NOTE: BED expansion does NOT require a reference FASTA because it operates
// on positional ranges, not allele sequences.  The "expansion" here simply
// widens the interval: new_start = start - MAX_SHIFT, new_end = end + MAX_SHIFT,
// and we rely on bedtools merge (or the in-memory sort-merge) to clean up.
//
// Implement as a flag: --expand-ambiguous-regions (default off).
// When on, each BED interval is padded by --indel-expand-padding (default 50).

pub(crate) const INDEL_EXPAND_PADDING_DEFAULT: usize = 50;

// ---------------------------------------------------------------------------
// build_diagnostic_store_expanded
// ---------------------------------------------------------------------------

/// Build `DiagnosticVariants` with indel equivalence expansion.
///
/// Parses a VCF/BCF of diagnostic positions and expands each indel variant
/// into all equivalent positions so that overlap detection works regardless
/// of aligner representation.
pub(crate) fn build_diagnostic_store_expanded<R: BufRead + Seek>(
    vcf_path:   &Path,
    expander:   &mut IndelEquivalenceExpander<R>,
    name_to_id: &HashMap<String, usize>,
    header:     &vcf::Header,
) -> Result<DiagnosticVariants> {
    use std::{fs::File, io::BufReader};

    let is_bcf = vcf_path.extension().map_or(false, |e| e == "bcf");
    let records: Box<dyn Iterator<Item = Result<vcf::variant::RecordBuf>>> = if is_bcf {
        let mut r = bcf::io::Reader::new(
            bgzf::io::Reader::new(BufReader::new(File::open(vcf_path)?))
        );
        let _h = r.read_header()?;
        Box::new(r.record_bufs(header).map(|r| r.map_err(anyhow::Error::from)))
    } else {
        let mut r = vcf::io::Reader::new(BufReader::new(File::open(vcf_path)?));
        let _h    = r.read_header()?;
        Box::new(r.record_bufs(header).map(|r| r.map_err(anyhow::Error::from)))
    };

    let mut per_ref: Vec<Vec<DiagnosticSite>> = Vec::new();
    let mut n_canonical = 0u64;
    let mut n_expanded  = 0u64;

    for rec_result in records {
        let mut rec  = rec_result?;
        let chrom    = rec.reference_sequence_name().to_string();
        let ref_id   = match name_to_id.get(&chrom) {
            Some(&id) => id,
            None => {
                warn!(chrom, "chromosome not in name_to_id map; skipping diagnostic record");
                continue;
            }
        };

        // Fetch alleles (ref + first alt).
        let ref_bytes: Vec<u8> = rec
            .reference_bases()
            .iter()
            .map(|b| b.to_ascii_uppercase())
            .collect();
        let alt_bytes: Vec<u8> = rec
            .alternate_bases()
            .iter()
            .next()
            .map(|a| a.unwrap_or(&[]).iter().map(|b| b.to_ascii_uppercase()).collect())
            .unwrap_or_default();

        if alt_bytes.is_empty() {
            warn!(chrom, "diagnostic record has no ALT allele; skipped");
            continue;
        }

        let pos_0based: usize = rec
            .variant_start()
            .transpose()
            .context("variant_start")?
            .map(|p| p.get() - 1)   // VCF 1-based → 0-based
            .unwrap_or(0);

        n_canonical += 1;

        let kind = classify(&ref_bytes, &alt_bytes);
        let sites: Vec<DiagnosticSite> = if matches!(kind, IndelKind::Snp | IndelKind::Complex) {
            if matches!(kind, IndelKind::Complex) {
                warn!(
                    chrom,
                    pos = pos_0based + 1,
                    "Complex diagnostic allele; stored as-is without expansion"
                );
            }
            vec![DiagnosticSite { pos: pos_0based, ref_len: ref_bytes.len() }]
        } else {
            let (ctx_start, ctx) = expander
                .fetch_context(&chrom, pos_0based, ref_bytes.len())
                .with_context(|| format!("diagnostic FASTA fetch at {chrom}:{}", pos_0based + 1))?;

            let equivalents =
                enumerate_equivalents(pos_0based, &ref_bytes, &alt_bytes, &ctx, ctx_start);

            debug!(
                chrom,
                pos = pos_0based + 1,
                n_equivalents = equivalents.len(),
                "Diagnostic indel expansion"
            );
            equivalents
                .iter()
                .map(|eq| DiagnosticSite { pos: eq.pos, ref_len: eq.ref_a.len() })
                .collect()
        };

        n_expanded += sites.len() as u64;

        if per_ref.len() <= ref_id {
            per_ref.resize_with(ref_id + 1, Vec::new);
        }
        per_ref[ref_id].extend(sites);
    }

    // Sort and dedup each chromosome bucket.
    for bucket in &mut per_ref {
        bucket.sort_unstable_by_key(|s| s.pos);
        bucket.dedup_by_key(|s| s.pos);
    }

    tracing::info!(
        vcf  = %vcf_path.display(),
        n_canonical,
        n_expanded,
        "Diagnostic variant store built with indel equivalence expansion"
    );

    Ok(DiagnosticVariants { per_ref })
}
