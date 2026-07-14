// src/variant/indel_equiv.rs
//
// Indel Equivalence Enumerator — Approach B implementation.
//
// Pre-computes all mathematically equivalent VCF representations of each
// indel by sliding it through local repeat regions.  The hot path
// (store.overlapping_multi / align_alt_to_read / wis_max_rescue_delta)
// is completely unchanged; the expansion happens once at startup.
//
// Coordinate contract (strict):
//   VCF input  — 1-based, inclusive (noodles Variant::variant_start()).
//   All output — 0-based (pos = vcf_pos - 1).

use crate::variant::name_to_id::header_name_to_id;
use crate::variant::Variant;
use crate::{
    variant::{
        population::{parse_population_record, Population},
        sample::{parse_sample_record, Sample as SampleVariant},
        store::Store,
    },
    Error,
};
use fasta::io::IndexedReader;
use noodles::{
    bgzf,
    core::{Position, Region},
    fasta,
    fasta::record::Sequence,
    vcf,
    vcf::variant::record::{
        samples::{Sample, Series},
        AlternateBases, ReferenceBases,
    },
    vcf::variant::RecordBuf,
};
use smallvec::SmallVec;
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader, Seek};
use std::path::Path;
use std::str::FromStr;
use tracing::{debug, warn};

// ---------------------------------------------------------------------------
// Constants
// ---------------------------------------------------------------------------

/// Hard limit on right-shift distance. Caps memory at ≤ MAX_SHIFT entries
/// per indel; handles homopolymers up to 100 bp.
pub(crate) const MAX_SHIFT: usize = 100;

/// Context fetched on each side of the indel for left-normalization and
/// right-shift enumeration. Slight over-fetch avoids off-by-one at limits.
const CTX_WINDOW: usize = MAX_SHIFT + 16;

// ---------------------------------------------------------------------------
// EquivalentAlleles — one equivalent VCF representation
// ---------------------------------------------------------------------------

/// One equivalent (POS, REF, ALT) representation of a biological indel,
/// stored in 0-based coordinates.
#[derive(Debug, Clone)]
pub(crate) struct EquivalentAlleles {
    /// 0-based anchor position (VCF POS − 1).
    pub(crate) pos: usize,
    /// REF allele bytes, ASCII uppercase, including anchor at index 0.
    pub(crate) ref_a: Vec<u8>,
    /// ALT allele bytes, ASCII uppercase, including anchor at index 0.
    pub(crate) alt_a: Vec<u8>,
}

// ---------------------------------------------------------------------------
// Allele classification
// ---------------------------------------------------------------------------

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum IndelKind {
    /// Both alleles are single bases; no sliding possible.
    Snp,
    /// len(ref) > len(alt) == 1 — pure deletion with VCF anchor base.
    Deletion,
    /// len(ref) == 1 < len(alt) — pure insertion with VCF anchor base.
    Insertion,
    /// Multi-base on both sides; sliding undefined.
    Complex,
}

pub(crate) fn classify(ref_a: &[u8], alt_a: &[u8]) -> IndelKind {
    match (ref_a.len(), alt_a.len()) {
        (1, 1) => IndelKind::Snp,
        (r, 1) if r > 1 => IndelKind::Deletion,
        (1, a) if a > 1 => IndelKind::Insertion,
        _ => IndelKind::Complex,
    }
}

// ---------------------------------------------------------------------------
// Core enumeration algorithm (pure function; no I/O)
// ---------------------------------------------------------------------------

/// Expand `(pos_0based, ref_a, alt_a)` into all equivalent VCF
/// representations by left-normalizing and then right-shifting through
/// the local repeat region.
///
/// # Parameters
/// - `pos_0based` : 0-based anchor position (VCF POS − 1).
/// - `ref_a`      : REF bytes including anchor; ASCII uppercase.
/// - `alt_a`      : ALT bytes including anchor; ASCII uppercase.
/// - `ref_ctx`    : Reference sequence bytes; `ref_ctx[0]` corresponds to
///                  chromosome position `ctx_start` (0-based).
///                  Must cover at least `pos_0based + len(ref_a) + MAX_SHIFT`.
/// - `ctx_start`  : 0-based chromosome position of `ref_ctx[0]`.
///
/// # Returns
/// `SmallVec` with inline capacity 8.  Always contains ≥ 1 entry
/// (the left-normalized original).  SNPs and complex alleles return
/// exactly 1 entry.
pub(crate) fn enumerate_equivalents(
    pos_0based: usize,
    ref_a: &[u8],
    alt_a: &[u8],
    ref_ctx: &[u8],
    ctx_start: usize,
) -> SmallVec<[EquivalentAlleles; 8]> {
    let mut out: SmallVec<[EquivalentAlleles; 8]> = SmallVec::new();
    let kind = classify(ref_a, alt_a);

    if matches!(kind, IndelKind::Snp | IndelKind::Complex) {
        // SNPs need no expansion; complex alleles cannot be safely slid.
        out.push(EquivalentAlleles {
            pos: pos_0based,
            ref_a: ref_a.to_vec(),
            alt_a: alt_a.to_vec(),
        });
        return out;
    }

    // -- Left-normalize ------------------------------------------------------
    let (mut cur_pos, mut cur_ref, mut cur_alt) =
        left_normalize(pos_0based, ref_a, alt_a, ref_ctx, ctx_start);

    out.push(EquivalentAlleles {
        pos: cur_pos,
        ref_a: cur_ref.clone(),
        alt_a: cur_alt.clone(),
    });

    // -- Right-shift enumeration ---------------------------------------------
    for _ in 0..MAX_SHIFT {
        let next = match kind {
            IndelKind::Deletion => {
                right_shift_deletion(cur_pos, &cur_ref, &cur_alt, ref_ctx, ctx_start)
            }
            IndelKind::Insertion => {
                right_shift_insertion(cur_pos, &cur_ref, &cur_alt, ref_ctx, ctx_start)
            }
            _ => unreachable!("Snp/Complex handled above"),
        };
        match next {
            Some((new_pos, new_ref, new_alt)) => {
                out.push(EquivalentAlleles {
                    pos: new_pos,
                    ref_a: new_ref.clone(),
                    alt_a: new_alt.clone(),
                });
                cur_pos = new_pos;
                cur_ref = new_ref;
                cur_alt = new_alt;
            }
            None => break,
        }
    }

    out
}

/// Left-normalize a simple indel.
///
/// Repeatedly shifts the anchor one position to the left as long as the
/// last base of the payload (deleted or inserted sequence) matches the
/// reference base immediately before the current anchor.
fn left_normalize(
    pos_0based: usize,
    ref_a: &[u8],
    alt_a: &[u8],
    ref_ctx: &[u8],
    ctx_start: usize,
) -> (usize, Vec<u8>, Vec<u8>) {
    let mut pos = pos_0based;
    let mut r = ref_a.to_vec();
    let mut a = alt_a.to_vec();

    loop {
        if pos == 0 {
            break;
        }

        // Chromosome index of the base preceding the current anchor.
        let prev_chrom = pos - 1;
        let prev_ctx = match prev_chrom.checked_sub(ctx_start) {
            Some(i) if i < ref_ctx.len() => ref_ctx[i],
            _ => break,
        };

        // Left-shift condition: payload's last base == preceding ref base.
        // For deletion payload is r[1:]; for insertion payload is a[1:].
        let last_payload = if r.len() > a.len() {
            *r.last().expect("ref non-empty")
        } else {
            *a.last().expect("alt non-empty")
        };

        if last_payload != prev_chrom_base_to_u8(prev_ctx) {
            break;
        }

        // Perform shift: prepend preceding base, drop the last base.
        r.insert(0, prev_ctx);
        r.pop();
        a.insert(0, prev_ctx);
        a.pop();
        pos -= 1;
    }
    (pos, r, a)
}

#[inline(always)]
const fn prev_chrom_base_to_u8(b: u8) -> u8 {
    b
} // identity; name for clarity

/// Attempt one right-shift of a deletion.
///
/// Condition: first deleted base == base immediately after the deletion
/// in the reference (i.e., `ref_ctx[anchor + 1] == ref_ctx[anchor + 1 + del_len]`).
///
/// New representation:
/// ```text
///   new_anchor   = cur_ref[1]               (first deleted base)
///   new_deleted  = cur_ref[2..] + ref_ctx[after_deletion]
///   new_ref      = [new_anchor] + new_deleted
///   new_alt      = [new_anchor]
/// ```
fn right_shift_deletion(
    cur_pos: usize,
    cur_ref: &[u8], // anchor + deleted sequence
    cur_alt: &[u8], // anchor only (len == 1)
    ref_ctx: &[u8],
    ctx_start: usize,
) -> Option<(usize, Vec<u8>, Vec<u8>)> {
    debug_assert_eq!(cur_alt.len(), 1, "deletion alt must be single anchor");
    let del_len = cur_ref.len().checked_sub(1)?; // number of deleted bases
    if del_len == 0 {
        return None;
    }

    // Chromosome coordinates (0-based):
    let first_del_chrom = cur_pos + 1;
    let after_del_chrom = cur_pos + 1 + del_len;

    // Map to ref_ctx slice indices.
    let first_del_idx = first_del_chrom.checked_sub(ctx_start)?;
    let after_del_idx = after_del_chrom.checked_sub(ctx_start)?;

    if first_del_idx >= ref_ctx.len() || after_del_idx >= ref_ctx.len() {
        return None;
    }

    // The shift is valid iff the first deleted base equals the base that
    // would "slide in" from the right.
    if ref_ctx[first_del_idx] != ref_ctx[after_del_idx] {
        return None;
    }

    let new_anchor = cur_ref[1]; // == ref_ctx[first_del_idx]

    let mut new_ref = Vec::with_capacity(cur_ref.len());
    new_ref.push(new_anchor);
    new_ref.extend_from_slice(&cur_ref[2..]); // remaining deleted bases
    new_ref.push(ref_ctx[after_del_idx]); // base that cycles in

    let new_alt = vec![new_anchor];

    Some((cur_pos + 1, new_ref, new_alt))
}

/// Attempt one right-shift of an insertion.
///
/// Condition: first inserted base == base immediately after the current anchor.
///
/// After a valid shift, the inserted sequence rotates left by one position.
/// ```text
///   new_anchor  = ref_ctx[cur_anchor + 1]
///   new_ins_seq = ins_seq[1..] + ins_seq[0]
///   new_alt     = [new_anchor] + new_ins_seq
///   new_ref     = [new_anchor]
/// ```
fn right_shift_insertion(
    cur_pos: usize,
    cur_ref: &[u8], // anchor only (len == 1)
    cur_alt: &[u8], // anchor + inserted sequence
    ref_ctx: &[u8],
    ctx_start: usize,
) -> Option<(usize, Vec<u8>, Vec<u8>)> {
    debug_assert_eq!(cur_ref.len(), 1, "insertion ref must be single anchor");
    let ins_len = cur_alt.len().checked_sub(1)?;
    if ins_len == 0 {
        return None;
    }

    let next_chrom = cur_pos + 1;
    let next_ctx_idx = next_chrom.checked_sub(ctx_start)?;
    if next_ctx_idx >= ref_ctx.len() {
        return None;
    }

    let first_ins = cur_alt[1]; // first inserted base (index 0 is anchor)
    let next_ref = ref_ctx[next_ctx_idx];

    // Shift valid iff first inserted base == next reference base.
    if first_ins != next_ref {
        return None;
    }

    let new_anchor = next_ref;

    let mut new_alt = Vec::with_capacity(cur_alt.len());
    new_alt.push(new_anchor);
    new_alt.extend_from_slice(&cur_alt[2..]); // ins_seq[1:]
    new_alt.push(cur_alt[1]); // ins_seq[0] rotated to end

    let new_ref = vec![new_anchor];

    Some((cur_pos + 1, new_ref, new_alt))
}

// ---------------------------------------------------------------------------
// IndelEquivalenceExpander — orchestrates FASTA fetch + enumeration
// ---------------------------------------------------------------------------

/// Wraps an indexed FASTA reader and expands each VCF indel into all
/// equivalent representations.
///
/// # Type parameter
/// `R` must implement `BufRead + Seek` so the FASTA index can random-access
/// the reference.  Typically `BufReader<File>`.
pub(crate) struct IndelEquivalenceExpander<R: BufRead + Seek> {
    fasta: IndexedReader<R>,
}

impl<R: BufRead + Seek> IndelEquivalenceExpander<R> {
    pub(crate) fn new(fasta: IndexedReader<R>) -> Self {
        Self { fasta }
    }

    /// Fetch the reference context window around a VCF indel.
    ///
    /// Returns `(ctx_start_0based, bytes)` where `bytes[0]` corresponds
    /// to chromosome position `ctx_start_0based` (0-based).
    pub(crate) fn fetch_context(
        &mut self,
        chrom: &str,
        pos_0based: usize,
        ref_len: usize,
    ) -> Result<(usize, Vec<u8>), Error> {
        // Start: ctx_start = max(0, pos_0based - CTX_WINDOW)
        let ctx_start = pos_0based.saturating_sub(CTX_WINDOW);
        // End: cover pos + ref_len + MAX_SHIFT (exclusive, 0-based).
        let ctx_end = pos_0based + ref_len + MAX_SHIFT + 1;

        let region_str = format!("{chrom}:{}-{}", ctx_start + 1, ctx_end);
        let region: Region = Region::from_str(&region_str)
            .map_err(|e| Error::InvalidRegion(format!("{region_str}: {e}")))?;
        let record = self.fasta.query(&region).map_err(|e| Error::FastaFetch {
            region: region_str,
            source: e,
        })?;

        let bytes: Vec<u8> = record
            .sequence()
            .as_ref()
            .iter()
            .map(|&b| b.to_ascii_uppercase())
            .collect();

        Ok((ctx_start, bytes))
    }

    /// Expand a single sample VCF record into all equivalent `SampleVariant`s.
    ///
    /// The original scoring parameters (`genotype_quality`, `alt_haplotype`,
    /// `phase_set`, `gamete_mode`) are replicated unchanged across all
    /// equivalent representations; only `pos`, `ref_a`, and `alt_a` change.
    pub(crate) fn expand_sample(
        &mut self,
        record: &mut RecordBuf,
        header: &vcf::Header,
        chrom: &str,
    ) -> Result<Vec<SampleVariant>, Error> {
        // Parse the canonical variant first (existing function; unchanged).
        let canonicals = match parse_sample_record(record, header) {
            Ok(v) => v,
            Err(e) => {
                warn!(chrom, "parse_sample_record failed: {e}; skipping record");
                return Ok(vec![]);
            }
        };

        let mut out: Vec<SampleVariant> = Vec::new();
        for &canon in &canonicals {
            let kind = classify(canon.ref_allele(), canon.alt_allele());
            if matches!(kind, IndelKind::Snp) {
                out.push(canon.clone());
                continue;
            }
            if matches!(kind, IndelKind::Complex) {
                warn!(
                    chrom,
                    pos = canon.pos() + 1, // report 1-based in warning
                    "Complex allele skipped for indel equivalence expansion"
                );
                out.push(canon.clone());
                continue;
            }

            let ref_len = canon.ref_allele().len();
            let (ctx_start, ctx) = self.fetch_context(chrom, canon.pos(), ref_len)?;

            let equivalents = enumerate_equivalents(
                canon.pos(),
                canon.ref_allele(),
                canon.alt_allele(),
                &ctx,
                ctx_start,
            );

            debug!(
                chrom,
                pos = canon.pos() + 1,
                n_equivalents = equivalents.len(),
                "Indel equivalence expansion"
            );

            for eq in &equivalents {
                // Clone with shifted coordinates; scoring parameters preserved.
                out.push(canon.with_alleles(eq.pos, &eq.ref_a, &eq.alt_a));
            }
        }
        Ok(out)
    }

    /// Expand a single population VCF record into all equivalent `Population`s.
    pub(crate) fn expand_population(
        &mut self,
        record: &mut RecordBuf,
        header: &vcf::Header,
        chrom: &str,
    ) -> Result<Vec<Population>, Error> {
        let canonicals = match parse_population_record(record, header) {
            Ok(v) => v,
            Err(e) => {
                warn!(chrom, "parse_population_record failed: {e}; skipping");
                return Ok(vec![]);
            }
        };

        let mut out: Vec<Population> = Vec::new();
        for &canon in &canonicals {
            let kind = classify(canon.ref_allele(), canon.alt_allele());
            if matches!(kind, IndelKind::Snp) {
                out.push(canon.clone());
                continue;
            }
            if matches!(kind, IndelKind::Complex) {
                warn!(
                    chrom,
                    pos = canon.pos() + 1,
                    "Complex allele skipped for indel equivalence expansion"
                );
                out.push(canon.clone());
                continue;
            }

            let ref_len = canon.ref_allele().len();
            let (ctx_start, ctx) = self.fetch_context(chrom, canon.pos(), ref_len)?;

            let equivalents = enumerate_equivalents(
                canon.pos(),
                canon.ref_allele(),
                canon.alt_allele(),
                &ctx,
                ctx_start,
            );

            debug!(
                chrom,
                pos = canon.pos() + 1,
                n_equivalents = equivalents.len(),
                "Population indel equivalence expansion"
            );

            for eq in &equivalents {
                out.push(canon.with_alleles(eq.pos, &eq.ref_a, &eq.alt_a));
            }
        }
        Ok(out)
    }
}

// ---------------------------------------------------------------------------
// `with_alleles` extension on variant types
// ---------------------------------------------------------------------------
// These must be implemented on SampleVariant and Population respectively.
// Shown here as a trait to avoid repeating the derive boilerplate.

pub(crate) trait WithAlleles: Sized {
    /// Clone self with new 0-based position and alleles; all scoring
    /// parameters (p_variant, GQ, AF, phase information) are preserved.
    fn with_alleles(&self, pos_0based: usize, ref_a: &[u8], alt_a: &[u8]) -> Self;
}

// ---------------------------------------------------------------------------
// Initialization loop (shown as a standalone function)
// ---------------------------------------------------------------------------
//
// Wire this into AlnStream::new() or the per-stream variant-store
// initialization block in main.rs where Store<Sample> / Store<Population>
// are currently constructed.

fn build_name_to_id(hdr: &vcf::Header) -> HashMap<String, usize> {
    hdr.contigs()
        .iter()
        .enumerate()
        .map(|(i, (name, _))| (name.to_string(), i))
        .collect()
}

/// Build a `Store<SampleVariant>` with indel equivalence expansion.
///
/// Replaces the bare `Store::from_vcf(path)` call in aln_stream.rs
/// when `--reference` is supplied.
///
/// # Coordinate contract
/// All positions in the returned `Store` are 0-based.  The existing
/// `store.overlapping_multi(tid, read_start, read_end)` hot path operates
/// on 0-based half-open intervals and is unchanged.
pub(crate) fn build_sample_store_expanded(
    vcf_path: &Path,
    fasta_path: &Path,
) -> Result<Store<SampleVariant>, Error> {
    use noodles::bgzf;
    use std::{fs::File, io::BufReader};
    let fasta_reader = fasta::io::indexed_reader::Builder::default().build_from_path(fasta_path)?;
    let mut expander = IndelEquivalenceExpander::new(fasta_reader);
    let mut store = Store::<SampleVariant>::new();

    // BCF (bgzf) or plain VCF — determine from extension.
    let is_bcf = vcf_path.extension().map_or(false, |e| e == "bcf");
    let hdr;

    // Declare uninitialized readers here so they live through the whole function
    let mut bcf_reader;
    let mut vcf_reader;
    // Add `+ '_` to tell Rust this Box borrows local variables
    let reader: Box<dyn Iterator<Item = Result<RecordBuf, Error>> + '_> = if is_bcf {
        use noodles::bcf;
        bcf_reader =
            bcf::io::Reader::new(bgzf::io::Reader::new(BufReader::new(File::open(vcf_path)?)));
        hdr = bcf_reader.read_header()?;
        Box::new(
            bcf_reader
                .record_bufs(&hdr)
                .map(|res| res.map_err(Error::from)),
        )
    } else {
        vcf_reader = vcf::io::Reader::new(BufReader::new(File::open(vcf_path)?));
        hdr = vcf_reader.read_header()?;
        Box::new(
            vcf_reader
                .record_bufs(&hdr)
                .map(|res| res.map_err(Error::from)),
        )
    };
    let name_to_id = build_name_to_id(&hdr);

    let mut n_canonical = 0u64;
    let mut n_expanded = 0u64;

    for rec_result in reader {
        let mut rec = rec_result?;
        // Derive chromosome name for FASTA fetch.
        let chrom = rec.reference_sequence_name().to_string();
        let ref_id = match name_to_id.get(&chrom) {
            Some(&id) => id,
            None => return Err(Error::NoReferenceSequenceId),
        };
        let variants = expander.expand_sample(&mut rec, &hdr, &chrom)?;
        n_expanded += variants.len() as u64;
        n_canonical += 1;
        store.insert_expanded(ref_id, variants);
    }

    tracing::info!(
        vcf  = %vcf_path.display(),
        n_canonical,
        n_expanded,
        "Sample variant store built with indel equivalence expansion"
    );
    Ok(store)
}

/// Build a `Store<Population>` with indel equivalence expansion.
pub(crate) fn build_population_store_expanded(
    vcf_path: &Path,
    fasta_path: &Path,
) -> Result<Store<Population>, Error> {
    use noodles::bgzf;
    use std::{fs::File, io::BufReader};

    // Validate that the .fai sidecar exists before opening the expander.
    let fasta_reader = fasta::io::indexed_reader::Builder::default().build_from_path(fasta_path)?;
    let mut expander = IndelEquivalenceExpander::new(fasta_reader);
    let mut store = Store::<Population>::new();

    let is_bcf = vcf_path.extension().map_or(false, |e| e == "bcf");
    let hdr;

    // Declare uninitialized readers here so they live through the whole function
    let mut bcf_reader;
    let mut vcf_reader;

    // Add `+ '_` to tell Rust this Box borrows local variables
    let reader: Box<dyn Iterator<Item = Result<RecordBuf, Error>> + '_> = if is_bcf {
        use noodles::bcf;
        bcf_reader =
            bcf::io::Reader::new(bgzf::io::Reader::new(BufReader::new(File::open(vcf_path)?)));
        hdr = bcf_reader.read_header()?;
        Box::new(
            bcf_reader
                .record_bufs(&hdr)
                .map(|res| res.map_err(Error::from)),
        )
    } else {
        vcf_reader = vcf::io::Reader::new(BufReader::new(File::open(vcf_path)?));
        hdr = vcf_reader.read_header()?;
        Box::new(
            vcf_reader
                .record_bufs(&hdr)
                .map(|res| res.map_err(Error::from)),
        )
    };
    let name_to_id = build_name_to_id(&hdr);

    let mut n_canonical = 0u64;
    let mut n_expanded = 0u64;
    for rec_result in reader {
        let mut rec = rec_result?;
        let chrom = rec.reference_sequence_name().to_string();
        let ref_id = match name_to_id.get(&chrom) {
            Some(&id) => id,
            None => return Err(Error::NoReferenceSequenceId),
        };
        let variants = expander.expand_population(&mut rec, &hdr, &chrom)?;
        n_expanded += variants.len() as u64;
        n_canonical += 1;
        store.insert_expanded(ref_id, variants);
    }

    tracing::info!(
        vcf  = %vcf_path.display(),
        n_canonical,
        n_expanded,
        "Population variant store built with indel equivalence expansion"
    );
    Ok(store)
}

#[cfg(test)]
mod tests;

#[cfg(test)]
mod integration;
