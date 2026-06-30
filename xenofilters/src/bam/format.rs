use clap::ValueEnum;

#[derive(Copy, Clone, Debug, ValueEnum, Default, PartialEq)]
pub(crate) enum BamFormat {
    #[default]
    Bam,
    Sam,
    Cram,
}
