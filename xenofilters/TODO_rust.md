# TODO list — xenofilter (Rust port)

## High priority — make code maintainable

- [✓] Split large files (line_by_line.rs → core/scoring/io + tests/ subdir)
- [✓] Move most unit tests into sibling tests/ subdirs

## Medium priority — features

- [✓] Add @PG line to all output BAM headers, unless disabled.
  → helper `add_pg_line(&mut header)` in bam_format.rs
- [✓] Optional BAM aux tag with decision phred score XF:C:
  → new flag `--add-decision-tag` (default off)
  → new flag `--ambiguous-threshold <phred-score>`
- [ ] Implement variant rescue (-s / -p)
  → parse VCF → per-position alt match bonus in scoring
  → start with SNVs only
- add config flag --variant-rescue-bonus f64 to control strength
- unify error types.


## Testing

- [ ] Finish integration_test.rs (use tempfile + Command)
- [ ] Add 4–6 realistic small BAM fixtures
- [ ] Compare output BAMs via samtools view | sort | diff
- [ ] Test with synthetic data: create small BAM + VCF where one alignment has a rescue-able SNV/indel
  → check score delta changes decision.

## Later / nice-to-have

- [ ] Progress bar (indicatif)
- [ ] Final summary counters (fragments processed / % best / ambig / discard)
- [ ] Publish to crates.io (once stable)
- [ ] --help grouping / examples
- [ ] Benchmark vs original XenofilteR / xenome
- [ ] Reduce test duplication (table-driven where possible)
