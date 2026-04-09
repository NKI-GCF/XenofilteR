use clap::ValueEnum;
use rust_htslib::bam::{self, Format, Header, HeaderView, header::HeaderRecord};
use rust_htslib::errors::Error;
use std::path::{Path, PathBuf};

#[derive(Copy, Clone, Debug, ValueEnum, Default, PartialEq)]
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
    path.as_ref().to_str().ok_or(Error::NonUnicodePath)?;
    Ok(())
}

fn format_from_extension(path: &Path) -> Format {
    let ext = path.extension().and_then(|e| e.to_str()).unwrap_or("bam");
    match ext.to_lowercase().as_str() {
        "sam" => Format::Sam,
        "cram" => Format::Cram,
        _ => Format::Bam,
    }
}

pub(crate) fn add_pg_line(header: &mut Header) {
    let mut pg = HeaderRecord::new(b"PG");
    let s = String::from("xenofilter");
    pg.push_tag(b"ID", &s);
    pg.push_tag(b"PN", &s);
    pg.push_tag(b"VN", &env!("CARGO_PKG_VERSION").to_string());
    pg.push_tag(b"CL", &std::env::args().collect::<Vec<_>>().join(" "));

    header.push_record(&pg);
}

pub(crate) fn out_from_file(
    f: &PathBuf,
    hdr_view: &HeaderView,
    add_pg: bool,
) -> Result<bam::Writer, Error> {
    let mut header = bam::Header::from_template(hdr_view);
    if add_pg {
        add_pg_line(&mut header);
    }
    let fmt = format_from_extension(f);
    bam::Writer::from_path(f, &header, fmt)
}

pub(crate) fn out_stdout(
    hdr_view: &HeaderView,
    fmt: Format,
    add_pg: bool,
) -> Result<bam::Writer, Error> {
    let mut header = bam::Header::from_template(hdr_view);
    if add_pg {
        add_pg_line(&mut header);
    }
    let mut ob: bam::Writer = bam::Writer::from_stdout(&header, fmt)?;
    if let Format::Bam = fmt {
        ob.set_compression_level(rust_htslib::bam::CompressionLevel::Fastest)?
    }
    Ok(ob)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_format_from_extension() {
        match format_from_extension(Path::new("file.bam")) {
            Format::Bam => (),
            _ => panic!("Expected Bam format"),
        }
        match format_from_extension(Path::new("file.sam")) {
            Format::Sam => (),
            _ => panic!("Expected Sam format"),
        }
        match format_from_extension(Path::new("file.cram")) {
            Format::Cram => (),
            _ => panic!("Expected Cram format"),
        }
        match format_from_extension(Path::new("file.unknown")) {
            Format::Bam => (),
            _ => panic!("Expected Bam format for unknown extension"),
        }
        match format_from_extension(Path::new("file")) {
            Format::Bam => (),
            _ => panic!("Expected Bam format for file with no extension"),
        }
    }
    #[test]
    fn test_path_unicode_ok() {
        assert!(path_unicode_ok("file.bam").is_ok());
        assert!(path_unicode_ok("file with spaces.bam").is_ok());
        assert!(path_unicode_ok("file_with_üñîçødé.bam").is_ok());
    }
    #[test]
    fn test_add_pg_line() {
        let mut header = Header::new();
        add_pg_line(&mut header);
        let pg_lines: Vec<_> = header.to_hashmap().get("PG").cloned().unwrap_or_default();
        assert_eq!(pg_lines.len(), 1);
        let pg_line = &pg_lines[0];
        assert_eq!(pg_line.get("ID").unwrap(), "xenofilter");
        assert_eq!(pg_line.get("PN").unwrap(), "xenofilter");
        assert_eq!(pg_line.get("VN").unwrap(), env!("CARGO_PKG_VERSION"));
        assert!(pg_line.get("CL").is_some());
    }
}
