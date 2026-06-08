mod population;
mod sample;
mod store;
mod eval;

pub(super) use population::{Population, parse_population_record};
pub(super) use sample::{Sample, parse_sample_record};
pub(super) use store::{Store, StoreTrait};
pub(crate) use eval::Eval;

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
