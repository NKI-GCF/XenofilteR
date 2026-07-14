pub mod eval;
pub(crate) mod population;
pub(crate) mod sample;
pub(crate) mod store;
pub(crate) mod indel_equiv;
pub(crate) mod indel_equiv_impls; // WithAlleles impls + Store::insert_expanded
pub(crate) mod indel_equiv_corrected;
pub(crate) mod diagnostic_equiv;
pub(crate) mod store_insert;
pub(crate) mod name_to_id;

pub(crate) use eval::Eval;
pub(crate) use population::{parse_population_record, Population};
pub(crate) use sample::{parse_sample_record, Sample};

pub(crate) use indel_equiv::{
    build_population_store_expanded, build_sample_store_expanded,
    enumerate_equivalents, IndelEquivalenceExpander, WithAlleles,
};
pub(crate) use diagnostic_equiv::build_diagnostic_store_expanded;
pub(super) use store::{Store, StoreTrait, VNT_CT};

use crate::filter_algorithm::line_by_line::READ_CT;
use smallvec::SmallVec;
use store::EvalVec;

pub(super) type FragEvalVec<'v> = SmallVec<[EvalVec<'v>; READ_CT]>;

/// Trait for any object that can be scored against an alignment.
pub trait Variant: Sync + Send {
    /// The reference ID of the variant.
    fn ref_id(&self) -> usize;

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
        let v_end = self.end();
        v_start < read_end && v_end > read_start
    }
    fn p_variant(&self) -> f64;
}
