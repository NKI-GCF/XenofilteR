mod format;
mod io;
mod merged_output;

pub(crate) use format::BamFormat;
pub(crate) use io::{out_from_file, path_unicode_ok, BamOutput, MergedOutput};
pub(crate) use merged_output::{expand_header, rewrite_rg, SUFFIX_AMBIGUOUS, SUFFIX_FILTERED};
