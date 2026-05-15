mod population_variant;
mod sample_variant;
mod variant_store;
mod variant_eval;

pub(super) use population_variant::{PopulationVariant, parse_population_record};
pub(super) use sample_variant::{SampleVariant, parse_sample_record};
pub(super) use variant_store::{VariantStore, VariantStoreTrait};
use smallvec::SmallVec;
pub(crate) use variant_eval::VariantEval;

pub(crate) type DeltaPerRec<'a> = SmallVec<[VariantEval<'a>; 0]>;
pub(crate) type VntPerRec<'a> = SmallVec<[DeltaPerRec<'a>; 2]>;

/// Trait for any object that can be scored against an alignment.
pub(crate) trait Variant: Sync + Send {
    /// The 1-based reference position of the variant.
    fn pos(&self) -> usize;

    /// The reference allele
    fn ref_allele(&self) -> &[u8];

    /// The alternate allele
    fn alt_allele(&self) -> &[u8];

    /// End position (1-based, exclusive) — allows easy overlap check
    fn end(&self) -> usize {
        self.pos() + self.ref_allele().len()
    }

    /// Does this variant overlap [read_start, read_end) ?
    fn overlaps(&self, read_start: usize, read_end: usize) -> bool {
        let v_start = self.pos();
        let v_end   = self.end();
        v_start < read_end && v_end > read_start
    }
    fn p_variant(&self) -> f64;
}
