use anyhow::{Result, anyhow};
use noodles::bam::io::Writer as BamWriter;
use noodles::sam::{Header, header::record::value::{Map, map::program::tag}};
use noodles::bgzf::io::{Writer as BgzfWriter, writer::{CompressionLevel, Builder as BgzfBuilder}};
use std::{fs::File, path::Path};

pub(crate) fn path_unicode_ok<'a, P: 'a + AsRef<Path>>(path: P) -> Result<()> {
    path.as_ref().to_str()
        .ok_or_else(|| anyhow!("Path {} is not valid UTF-8", path.as_ref().display()))?;
    Ok(())
}

fn add_pg_line(header: &mut Header) -> Result<()>{
    let id = "xenofilter";

    let program = Map::builder()
        .insert(tag::NAME, "xenofilter")
        .insert(tag::VERSION, env!("CARGO_PKG_VERSION"))
        .insert(tag::COMMAND_LINE, std::env::args().collect::<Vec<_>>().join(" "))
        .build()
        .expect("Failed to build PG record");

    // Insert it into the header's program map
    header.programs_mut().add(id, program)?;
    Ok(())
}

pub(crate) fn out_from_file(
    f: &Path,
    header: &Header,
    add_pg: bool,
) -> Result<BamWriter<BgzfWriter<File>>> {
    let file = File::create(f)?;
    let mut modified_header = header.clone();
    if add_pg {
        add_pg_line(&mut modified_header)?;
    }

    let encoder = BgzfBuilder::default()
        .set_compression_level(CompressionLevel::FAST)
        .build_from_writer(file);
    let mut writer = BamWriter::from(encoder);
    writer.write_header(&modified_header)?;

    Ok(writer)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_path_unicode_ok() {
        assert!(path_unicode_ok("file.bam").is_ok());
        assert!(path_unicode_ok("file with spaces.bam").is_ok());
        assert!(path_unicode_ok("file_with_üñîçødé.bam").is_ok());
    }
    #[test]
    fn test_add_pg_line() {
        let mut header = Header::default();
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
