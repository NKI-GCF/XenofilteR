// src/region/tabix_load.rs
use std::path::Path;
use crate::Error;
use crate::region::tabix_query::{TabixBed, TabixScored, TabixVcf};
use crate::region::{ScoreFn, load::load_pair};

pub(crate) fn open_tabix_bed(specs: &[String]) -> Result<[Option<TabixBed>; 2], Error> {
    load_pair(specs, |p| TabixBed::open(Path::new(p)))
}
pub(crate) fn open_tabix_vcf(specs: &[String]) -> Result<[Option<TabixVcf>; 2], Error> {
    load_pair(specs, |p| TabixVcf::open(Path::new(p)))
}
pub(crate) fn open_tabix_scored(
    specs: &[String],
    score_fn: ScoreFn,
) -> Result<[Option<(TabixScored, ScoreFn)>; 2], Error> {
    let mut out: [Option<(TabixScored, ScoreFn)>; 2] = [None, None];
    for (i, slot) in out.iter_mut().enumerate() {
        if let Some(s) = specs.get(i).filter(|s| !s.is_empty()) {
            *slot = Some((TabixScored::open(Path::new(s))?, score_fn));
        }
    }
    Ok(out)
}
