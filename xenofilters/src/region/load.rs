// src/region/load.rs

pub(crate) fn load_ambiguous_regions_memory(
    specs: &[String],
    name_to_id: &HashMap<String, usize>,
) -> Result<[Option<AmbiguousRegions>; 2], Error> {
    load_pair(specs, |p| {
        AmbiguousRegions::from_bed(Path::new(p), name_to_id)
    })
}

pub(crate) fn load_diagnostic_variants_memory(
    specs: &[String],
    name_to_id: &HashMap<String, usize>,
) -> Result<[Option<DiagnosticVariants>; 2], Error> {
    load_pair(specs, |p| {
        DiagnosticVariants::from_vcf(Path::new(p), name_to_id)
    })
}

pub(crate) fn load_positive_regions_memory(
    specs: &[String],
    name_to_id: &HashMap<String, usize>,
    score_fn: ScoreFn,
) -> Result<[Option<(ScoredRegions, ScoreFn)>; 2], Error> {
    let mut out: [Option<(ScoredRegions, ScoreFn)>; 2] = [None, None];
    for (i, slot) in out.iter_mut().enumerate() {
        if let Some(s) = specs.get(i).filter(|s| !s.is_empty()) {
            *slot = Some((ScoredRegions::from_bed(Path::new(s), name_to_id)?, score_fn));
        }
    }
    Ok(out)
}

/// Shared 2-slot loader: applies `f` to each of up to 2 file paths.
fn load_pair<T, F>(specs: &[String], f: F) -> Result<[Option<T>; 2], Error>
where
    F: Fn(&str) -> Result<T, Error>,
{
    let mut out: [Option<T>; 2] = [None, None];
    for (i, slot) in out.iter_mut().enumerate() {
        if let Some(s) = specs.get(i).filter(|s| !s.is_empty()) {
            *slot = Some(f(s)?);
        }
    }
    Ok(out)
}
