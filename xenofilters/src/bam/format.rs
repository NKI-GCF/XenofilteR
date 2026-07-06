use clap::ValueEnum;

#[derive(Copy, Clone, Debug, ValueEnum, Default, PartialEq)]
pub enum AlnFormat {
    #[default]
    Bam,
    Sam,
    Cram,
}
