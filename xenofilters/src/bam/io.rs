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
    fn test_add_pg_line() -> Result<()> {
        let mut header = Header::default();
        add_pg_line(&mut header)?;
        let pg = header.programs();
        let mut roots = pg.roots();
        let (id, map): (&[u8], &_) = roots.next().map(|(id, map)| (id.as_ref(), map)).expect("No PG");
        assert_eq!(id, b"xenofilter");

        let of = map.other_fields();
        let vn = of.get(&tag::VERSION).expect("VN").to_string();
        let vn = vn.split('.').collect::<Vec<_>>();
        let env_vn = env!("CARGO_PKG_VERSION").split('.').collect::<Vec<_>>();
        assert_eq!(vn[0].parse::<u32>().expect("maj1"), env_vn[0].parse::<u32>().expect("maj2"));
        assert_eq!(vn[1].parse::<u32>().expect("min1"), env_vn[1].parse::<u32>().expect("min2"));

        let cl_value = of.get(&tag::COMMAND_LINE).expect("CL tag not found");
        assert_eq!(cl_value.to_string(), std::env::args().collect::<Vec<_>>().join(" "));
        assert_eq!(roots.next(), None);
        Ok(())
    }
}
