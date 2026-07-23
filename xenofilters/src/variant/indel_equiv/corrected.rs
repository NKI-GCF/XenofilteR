// src/variant/indel_equiv.rs  (replace the build_*_store_expanded functions
//   at the bottom with these; the enumeration core above is unchanged)
//
// Key corrections vs the previous version:
//   1. All `store.insert_expanded(variants)` calls now pass `ref_id` explicitly.
//   2. `expand_sample` / `expand_population` return `(ref_id, Vec<V>)` tuples.
//   3. `IndelEquivalenceExpander::fetch_context` is `pub(crate)` so
//      `diagnostic_equiv.rs` can call it.
//   4. The name->id mapping is passed into both build functions.

use std::{
    collections::HashMap,
    io::{BufRead, Seek},
    path::Path,
};

use noodles::{bcf, bgzf, fasta, vcf};
use tracing::{debug, warn};

use crate::variant::{Variant, population::Population, sample::Sample, store::Store};

use crate::Error;
use crate::variant::indel_equiv::{
    IndelEquivalenceExpander, IndelKind, classify, enumerate_equivalents,
};

// ---------------------------------------------------------------------------
// Corrected IndelEquivalenceExpander expand methods
// ---------------------------------------------------------------------------

impl<R: BufRead + Seek> IndelEquivalenceExpander<R> {
    /// Expand one VCF record into `(ref_id, Vec<Sample>)`.
    ///
    /// Returns `None` when the chromosome is absent from `name_to_id`
    /// (logged as a warning; record skipped without error).
    pub(crate) fn expand_sample_with_refid(
        &mut self,
        rec: &mut vcf::variant::RecordBuf,
        header: &vcf::Header,
        name_to_id: &HashMap<String, usize>,
        _gamete: bool,
    ) -> Result<Option<(usize, Vec<Sample>)>, Error> {
        let chrom = rec.reference_sequence_name().to_string();
        let ref_id = match name_to_id.get(&chrom) {
            Some(&id) => id,
            None => {
                warn!(
                    chrom,
                    "chromosome absent from reference header; skipping sample record"
                );
                return Ok(None);
            }
        };
        let canonicals = match crate::variant::sample::parse_sample_record(rec, header) {
            Ok(v) => v,
            Err(e) => {
                warn!(chrom, "parse_sample_record: {e}; skipping");
                return Ok(None);
            }
        };
        let out =
            self.expand_with_refid_inner(canonicals, &chrom, ref_id, "Sample indel expansion")?;
        Ok(Some((ref_id, out)))
    }

    /// Expand one VCF record into `(ref_id, Vec<Population>)`.
    pub(crate) fn expand_population_with_refid(
        &mut self,
        rec: &mut vcf::variant::RecordBuf,
        header: &vcf::Header,
        name_to_id: &HashMap<String, usize>,
    ) -> Result<Option<(usize, Vec<Population>)>, Error> {
        let chrom = rec.reference_sequence_name().to_string();
        let ref_id = match name_to_id.get(&chrom) {
            Some(&id) => id,
            None => {
                warn!(chrom, "chromosome absent; skipping population record");
                return Ok(None);
            }
        };
        let canonicals = match crate::variant::population::parse_population_record(rec, header) {
            Ok(v) => v,
            Err(e) => {
                warn!(chrom, "parse_population_record: {e}; skipping");
                return Ok(None);
            }
        };
        let out =
            self.expand_with_refid_inner(canonicals, &chrom, ref_id, "Population indel expansion")?;
        Ok(Some((ref_id, out)))
    }

    fn expand_with_refid_inner<V>(
        &mut self,
        canonicals: Vec<V>,
        chrom: &str,
        ref_id: usize,
        expand_label: &str,
    ) -> Result<Vec<V>, Error>
    where
        V: Clone + crate::variant::Variant + WithAllelesRefId,
    {
        let mut out = Vec::new();
        for canon in &canonicals {
            let kind = classify(canon.ref_allele(), canon.alt_allele());
            if matches!(kind, IndelKind::Complex) {
                warn!(chrom, pos = canon.pos() + 1, "Complex allele not expanded");
                out.push(canon.with_alleles_refid(
                    ref_id,
                    canon.pos(),
                    canon.ref_allele(),
                    canon.alt_allele(),
                ));
                continue;
            }
            if matches!(kind, IndelKind::Snp) {
                out.push(canon.with_alleles_refid(
                    ref_id,
                    canon.pos(),
                    canon.ref_allele(),
                    canon.alt_allele(),
                ));
                continue;
            }
            let (ctx_start, ctx) =
                self.fetch_context(chrom, canon.pos(), canon.ref_allele().len())?;
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
                n = equivalents.len(),
                "{expand_label}"
            );
            for eq in &equivalents {
                out.push(canon.with_alleles_refid(ref_id, eq.pos, &eq.ref_a, &eq.alt_a));
            }
        }
        Ok(out)
    }
}

fn build_expanded_store_corrected<V>(
    vcf_path: &Path,
    log_label: &str,
    mut expand: impl FnMut(
        &mut vcf::variant::RecordBuf,
        &vcf::Header,
    ) -> Result<Option<(usize, Vec<V>)>, Error>,
) -> Result<Store<V>, Error>
where
    V: Variant + Clone,
{
    let mut store = Store::<V>::new();
    let header = read_vcf_or_bcf_header(vcf_path)?;
    let mut n_canonical = 0u64;
    let mut n_expanded = 0u64;
    for_each_vcf_record(vcf_path, &header, |rec| {
        n_canonical += 1;
        if let Some((ref_id, variants)) = expand(rec, &header)? {
            n_expanded += variants.len() as u64;
            store.insert_expanded(ref_id, variants);
        }
        Ok(())
    })?;
    store.dedup();
    tracing::info!(vcf = %vcf_path.display(), n_canonical, n_expanded, "{log_label}");
    Ok(store)
}

// ---------------------------------------------------------------------------
// WithAllelesRefId -- variant of WithAlleles that also sets ref_id
// ---------------------------------------------------------------------------
//
// Avoids two separate clone+update steps. Implemented on Sample + Population
// in indel_equiv_impls.rs.

pub(crate) trait WithAllelesRefId: Sized {
    fn with_alleles_refid(
        &self,
        ref_id: usize,
        pos_0based: usize,
        ref_a: &[u8],
        alt_a: &[u8],
    ) -> Self;
}

impl WithAllelesRefId for Sample {
    fn with_alleles_refid(
        &self,
        ref_id: usize,
        pos_0based: usize,
        ref_a: &[u8],
        alt_a: &[u8],
    ) -> Self {
        Self {
            ref_id,
            pos: pos_0based,
            ref_a: ref_a.to_vec(),
            alt_a: alt_a.to_vec(),
            genotype_quality: self.genotype_quality,
            is_called: self.is_called,
        }
    }
}

impl WithAllelesRefId for Population {
    fn with_alleles_refid(
        &self,
        ref_id: usize,
        pos_0based: usize,
        ref_a: &[u8],
        alt_a: &[u8],
    ) -> Self {
        Self {
            ref_id,
            pos: pos_0based,
            ref_a: ref_a.to_vec(),
            alt_a: alt_a.to_vec(),
            allele_frequency: self.allele_frequency,
        }
    }
}

// ---------------------------------------------------------------------------
// build_sample_store_expanded  (corrected)
// ---------------------------------------------------------------------------

pub(crate) fn build_sample_store_expanded(
    vcf_path: &Path,
    fasta_path: &Path,
    name_to_id: &HashMap<String, usize>,
    gamete: bool,
) -> Result<Store<Sample>, Error> {
    let fasta_reader = fasta::io::indexed_reader::Builder::default().build_from_path(fasta_path)?;
    let mut expander = IndelEquivalenceExpander::new(fasta_reader);
    build_expanded_store_corrected(
        vcf_path,
        "Sample store built with indel equivalence expansion",
        |rec, hdr| expander.expand_sample_with_refid(rec, hdr, name_to_id, gamete),
    )
}

// ---------------------------------------------------------------------------
// build_population_store_expanded  (corrected)
// ---------------------------------------------------------------------------

pub(crate) fn build_population_store_expanded(
    vcf_path: &Path,
    fasta_path: &Path,
    name_to_id: &HashMap<String, usize>,
) -> Result<Store<Population>, Error> {
    let fasta_reader = fasta::io::indexed_reader::Builder::default().build_from_path(fasta_path)?;
    let mut expander = IndelEquivalenceExpander::new(fasta_reader);
    build_expanded_store_corrected(
        vcf_path,
        "Population store built with indel equivalence expansion",
        |rec, hdr| expander.expand_population_with_refid(rec, hdr, name_to_id),
    )
}

// ---------------------------------------------------------------------------
// Shared VCF/BCF iteration helpers (private)
// ---------------------------------------------------------------------------

/// Read header from either BCF or plain VCF (determined by extension).
pub(crate) fn read_vcf_or_bcf_header(path: &Path) -> Result<vcf::Header, Error> {
    use std::{fs::File, io::BufReader};
    if path.extension().is_some_and(|e| e == "bcf") {
        let mut r = bcf::io::Reader::new(bgzf::io::Reader::new(BufReader::new(File::open(path)?)));
        Ok(r.read_header()?)
    } else {
        let mut r = vcf::io::Reader::new(BufReader::new(File::open(path)?));
        Ok(r.read_header()?)
    }
}

/// Iterate over every record in a VCF or BCF file, calling `f` for each.
fn for_each_vcf_record<F>(path: &Path, header: &vcf::Header, mut f: F) -> Result<(), Error>
where
    F: FnMut(&mut vcf::variant::RecordBuf) -> Result<(), Error>,
{
    use std::{fs::File, io::BufReader};
    if path.extension().is_some_and(|e| e == "bcf") {
        let mut r = bcf::io::Reader::new(bgzf::io::Reader::new(BufReader::new(File::open(path)?)));
        let _h = r.read_header()?;
        for result in r.record_bufs(header) {
            let mut rec = result?;
            f(&mut rec)?;
        }
    } else {
        let mut r = vcf::io::Reader::new(BufReader::new(File::open(path)?));
        let _h = r.read_header()?;
        for result in r.record_bufs(header) {
            let mut rec = result?;
            f(&mut rec)?;
        }
    }
    Ok(())
}
