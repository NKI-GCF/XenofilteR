use clap::ValueEnum;
use rust_htslib::errors::Error;
use rust_htslib::bam::{self, Format, HeaderView};
use std::path::{PathBuf, Path};

#[derive(Copy, Clone, Debug, ValueEnum, Default)]
pub enum BamFormat {
    #[default]
    Bam,
    Sam,
    Cram,
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

pub fn path_unicode_ok<'a, P: 'a + AsRef<Path>>(path: P) -> Result<(), Error> {
    path.as_ref()
        .to_str()
        .ok_or(Error::NonUnicodePath)?;
    Ok(())
}


fn format_from_extension(path: &PathBuf) -> Format {
    let ext = path.extension().and_then(|e| e.to_str()).unwrap_or("bam");
    match ext.to_lowercase().as_str() {
        "sam" => Format::Sam,
        "cram" => Format::Cram,
        _ => Format::Bam,
    }
}

pub(crate) fn out_from_file(f: &PathBuf, hdr_view: &HeaderView) -> Result<bam::Writer, Error> {
    let hdr = bam::Header::from_template(hdr_view);
    let fmt = format_from_extension(f);
    bam::Writer::from_path(f, &hdr, fmt)
}

pub(crate) fn out_stdout(hdr_view: &HeaderView, fmt: Format) -> Result<bam::Writer, Error> {
    let hdr = bam::Header::from_template(hdr_view);
    let mut ob: bam::Writer = bam::Writer::from_stdout(&hdr, fmt)?;
    match fmt {
        Format::Bam => ob.set_compression_level(rust_htslib::bam::CompressionLevel::Fastest)?,
        _ => {}
    }
    Ok(ob)
}
