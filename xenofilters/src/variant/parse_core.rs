// src/variant/parse_core.rs
use crate::Error;
use noodles::vcf::{
    variant::record::AlternateBases,
    variant::RecordBuf,
    Header,
};

/// Fields common to every VCF-derived variant record, extracted once.
pub(crate) struct VariantCore {
    pub(crate) ref_id: usize,
    pub(crate) pos: usize, // 1-based, as returned by variant_start()
    pub(crate) ref_a: Vec<u8>,
    pub(crate) alts: Vec<Vec<u8>>,
}

pub(crate) fn parse_variant_core(
    record: &RecordBuf,
    header: &Header,
) -> Result<VariantCore, Error> {
    let chrom = record.reference_sequence_name().to_string();
    let ref_id = header.contigs().get_index_of(&chrom).ok_or(Error::NoRefSeqId)?;
    let pos = record.variant_start().map(|p| p.get()).unwrap_or(0);
    let ref_a = record.reference_bases().as_bytes().to_vec();
    let alts = record
        .alternate_bases()
        .iter()
        .map(|a| Ok(a?.as_bytes().to_vec()))
        .collect::<Result<Vec<_>, Error>>()?;
    Ok(VariantCore { ref_id, pos, ref_a, alts })
}
