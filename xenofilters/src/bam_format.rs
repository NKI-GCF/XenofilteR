use std::path::PathBuf;
use std::str::FromStr;

use anyhow::{Result, anyhow};
use clap::ValueEnum;
use rust_htslib::bam::{self, Format, HeaderView};

#[derive(Copy, Clone, Debug, ValueEnum)]
pub enum BamFormat {
    Bam,
    Sam,
    Cram,
}

impl FromStr for BamFormat {
    type Err = String;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s.to_lowercase().as_str() {
            "sam" => Ok(BamFormat::Sam),
            "cram" => Ok(BamFormat::Cram),
            _ => Ok(BamFormat::Bam),
        }
    }
}

impl From<BamFormat> for Format {
    fn from(f: BamFormat) -> Self {
        match f {
            BamFormat::Bam => Format::Bam,
            BamFormat::Sam => Format::Sam,
            BamFormat::Cram => Format::Cram,
        }
    }
}

pub(crate) fn out_from_file(f: &PathBuf, hdr_view: &HeaderView) -> Result<bam::Writer> {
    let f_str = f.extension().and_then(|e| e.to_str()).unwrap_or("bam");
    let fmt = <BamFormat as FromStr>::from_str(f_str)
        .map_err(|e| anyhow!(e))?
        .into();
    let hdr = bam::Header::from_template(hdr_view);
    Ok(bam::Writer::from_path(f, &hdr, fmt)?)
}

pub(crate) fn out_stdout(hdr_view: &HeaderView, fmt: Format) -> Result<bam::Writer> {
    let hdr = bam::Header::from_template(hdr_view);
    let mut ob = bam::Writer::from_stdout(&hdr, fmt)?;
    if let bam::Format::Bam = fmt {
        ob.set_compression_level(rust_htslib::bam::CompressionLevel::Fastest)?;
    }
    Ok(ob)
}
