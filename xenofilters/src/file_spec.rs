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
        // Try `IDX:FILE` first
        if let Some((idx_str, path_str)) = s.split_once(':') {
            let idx_str = idx_str.trim();
            let path_str = path_str.trim();

            if !path_str.is_empty() && (path_str.starts_with('/') || path_str.starts_with('.')) {
                // Treat as plain path, not an indexed spec
                return Ok(FileSpec {
                    idx: None,
                    path: PathBuf::from(s.trim()),
                });
            }

            let idx = idx_str
                .parse::<usize>()
                .map_err(|e| format!("Invalid stream index `{idx_str}`: {e}"))?;

            Ok(FileSpec {
                idx: Some(idx),
                path: PathBuf::from(path_str),
            })
        } else {
            // Plain FILE → positional spec
            Ok(FileSpec {
                idx: None,
                path: PathBuf::from(s.trim()),
            })
        }
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

