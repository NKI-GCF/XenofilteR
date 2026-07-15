pub(crate) mod ambiguous_vcf;
pub mod collated;
pub mod hash_lookup;
pub mod line_by_line;
pub mod strain;
pub mod viral_integration;

pub(crate) use line_by_line::core::COUNTER_STRIDE;
