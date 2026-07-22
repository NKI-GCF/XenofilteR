use std::path::PathBuf;
use std::str::FromStr;

#[derive(Debug, Clone)]
pub(crate) struct FileSpec {
    pub idx: Option<usize>,
    pub path: PathBuf,
}

impl FromStr for FileSpec {
    type Err = String; // or a custom error type

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        if let Some((idx_str, path_str)) = s.split_once(':')
            && let Ok(idx) = idx_str.trim().parse::<usize>() {
                return Ok(FileSpec {
                    idx: Some(idx),
                    path: PathBuf::from(path_str.trim()),
                });
            }
        Ok(FileSpec {
            idx: None,
            path: PathBuf::from(s.trim()),
        })
    }
}

pub(crate) fn path_for_stream(specs: &[FileSpec], stream_idx: usize) -> Option<&PathBuf> {
    if let Some(spec) = specs.iter().find(|spec| spec.idx == Some(stream_idx)) {
        return Some(&spec.path);
    }
    specs
        .iter()
        .filter(|spec| spec.idx.is_none())
        .nth(stream_idx)
        .map(|spec| &spec.path)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn plain_relative_path_no_colon() {
        let fs: FileSpec = "a.vcf".parse().unwrap();
        assert_eq!(fs.idx, None);
        assert_eq!(fs.path, PathBuf::from("a.vcf"));
    }

    #[test]
    fn indexed_spec_simple_filename_works() {
        // Matches usage in tests.rs / config tests: "0:a.vcf"
        let fs: FileSpec = "0:a.vcf".parse().unwrap();
        assert_eq!(fs.idx, Some(0));
        assert_eq!(fs.path, PathBuf::from("a.vcf"));
    }

    #[test]
    fn invalid_index_is_rejected() {
        let fs = "abc:file.bam".parse::<FileSpec>();
        assert!(fs.is_ok(), "path should not cause parse error");
        assert!(fs.unwrap().idx.is_none(), "invalid index should be None");
    }

    #[test]
    fn absolute_plain_path_no_colon() {
        let fs: FileSpec = "/abs/path.bam".parse().unwrap();
        assert_eq!(fs.idx, None);
        assert_eq!(fs.path, PathBuf::from("/abs/path.bam"));
    }

    #[test]
    fn indexed_spec_with_absolute_path() {
        let fs: FileSpec = "1:/abs/path.bam".parse().unwrap();
        assert_eq!(fs.idx, Some(1), "indexed absolute path must set idx");
        assert_eq!(fs.path, PathBuf::from("/abs/path.bam"));
    }

    #[test]
    fn indexed_spec_with_relative_dot_path() {
        let fs: FileSpec = "1:./file.bam".parse().unwrap();
        assert_eq!(fs.idx, Some(1));
        assert_eq!(fs.path, PathBuf::from("./file.bam"));
    }

    #[test]
    fn path_for_stream_prefers_explicit_index_over_positional() {
        let specs = vec![
            FileSpec {
                idx: None,
                path: "positional0.vcf".into(),
            },
            FileSpec {
                idx: Some(2),
                path: "explicit2.vcf".into(),
            },
            FileSpec {
                idx: None,
                path: "positional1.vcf".into(),
            },
        ];
        assert_eq!(
            path_for_stream(&specs, 2),
            Some(&PathBuf::from("explicit2.vcf"))
        );
        // stream 1 has no explicit match; falls back to the 2nd unindexed spec (0-idx'd: nth(1))
        assert_eq!(
            path_for_stream(&specs, 1),
            Some(&PathBuf::from("positional1.vcf"))
        );
        assert_eq!(path_for_stream(&specs, 5), None);
    }
}
