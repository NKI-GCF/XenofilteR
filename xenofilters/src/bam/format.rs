use clap::ValueEnum;

#[derive(Copy, Clone, Debug, ValueEnum, Default, PartialEq)]
pub(crate) enum AlnFormat {
    #[default]
    Bam,
    Sam,
    Cram,
}
