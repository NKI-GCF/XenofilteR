# TODO list — xenofilter (Rust port)

## Completed

- [✓] Split large files (line_by_line.rs → core/score/io + tests/ subdir)
- [✓] Move most unit tests into sibling tests/ subdirs
- [✓] Add @PG line to all output BAM headers (`add_pg_line` in bam/io.rs)
- [✓] Optional BAM aux tag with decision phred score `XF:C:` (`--add-decision-tag`)
- [✓] `--ambiguous-threshold <phred>` flag
- [✓] Implement variant rescue (`-s` / `-p`): VCF → per-position alt match bonus
- [✓] Three backends wired and dispatching: `LineByLine`, `HashLookup`, `Collated`
- [✓] `EarlyKind` enum (`AllUnmapped`, `AllPerfect`) in `HashLookup` assembly
- [✓] `resolve_fragment` full 2×2 decision matrix over `EarlyKind` combinations
- [✓] Supplementary alignment structural penalties in sequential, parallel, and HashLookup paths
- [✓] CIGAR/MD subsumption pre-check (Tier 2.5) in `LineByLine::resolve`
- [✓] Parallel scoring pipeline: crossbeam-channel bounded channels, per-worker `Scratch`
- [✓] `OutputMode` enum (`MultiFile` / `Merged`) with `expand_header()` and `rewrite_rg()`
- [✓] `AmbiguousVcfEvaluator` trait stub (`src/filter_algorithm/ambiguous_vcf.rs`)
- [✓] `--ambiguous-regions` BED and `--diagnostic-variants` VCF CLI flags
- [✓] `AmbiguousRegions` and `DiagnosticVariants` wired into `HashLookup` and `Collated`
- [✓] Tabix-indexed BED/VCF (`TabixBed`, `TabixVcf`) for `Collated`
- [✓] Variant stores refactored to `Option<Arc<dyn StoreTrait>>` (`Send + Sync`)
- [✓] `dvnt` double-population bug fixed
- [✓] `finished: FragEvalVec` accumulator to prevent silent delta discard in `maximize_delta`
- [✓] `eff_ref_start` anchoring fix in `score_variant_in_seg`
- [✓] `MdCigFlags` precomputation threading (no redundant `Box<dyn Cigar>` allocations)
- [✓] Dead-code pruned: `BaseOp::Relocate`, `early_assign` module, `MergedOutput::header`,
       `alt` field from `DiagnosticSite`, `reader` field from `TabixBed`,
       `can_early_assign`/`is_complete` from `PendingFragment`, `primary_count` from `StreamKind`

## Testing

- [ ] Finish `tests/integration_test.rs` (tempfile + Command)
- [ ] Add 4–6 realistic small BAM fixtures (paired-end, supplementary, unmapped combos)
- [ ] Compare output BAMs via `samtools view | sort | diff`
- [ ] Synthetic BAM + VCF fixture exercising variant rescue decision change
- [ ] HashLookup integration test with position-sorted BAM pair
- [ ] Collated integration test with name-skewed BAM pair

## Medium priority — features

- [ ] `--variant-rescue-bonus f64` config flag to control rescue strength
- [ ] Unify error types: consolidate `AlignmentError`, `anyhow::Error`, `std::io::Error`
       into a single `xenofilters::Error` enum via `thiserror`
- [ ] Multiple ALT alleles: extend `parse_sample_record` and `parse_population_record`
       to iterate all ALT alleles and pick the highest-scoring one per read
- [ ] Implement `AmbiguousVcfEvaluator::evaluate_ambiguous_variants` (currently a stub)

## Later / nice-to-have

- [ ] Final summary counters (fragments processed / % best / ambig / discard)
- [ ] Publish to crates.io (once variant rescue and test suite are complete)
- [ ] `--help` grouping / examples (clap `next_help_heading`)
- [ ] Benchmark vs original XenofilteR / xenome
- [ ] Progress report (completion % from BAM file size vs position, updated every 100k fragments)
- [ ] bgzf read parallelism via `.set_worker_count()` on the bgzf reader builder
- [ ] Reduce `RecordBuf` allocations on the hot path (profile first)
- [ ] Table-driven unit tests (reduce duplication in fragment_state/tests.rs and ops/tests.rs)
- [ ] N-stream support beyond 2: Namesorted and Collated are architecturally ready;
       HashLookup is flagged high-memory risk for N > 2
- [ ] HashLookup seek threading: offload `fetch_by_virtual_offset` to a dedicated
       seek-IO thread via crossbeam_channel (stub documented in stage.rs)
- [ ] Collated worker pool: dispatch `score_pair` to a thread pool
       (stub documented in collated.rs)
- [ ] GitHub Actions CI (Linux x86_64, Linux aarch64, macOS)
- [ ] Conda / Bioconda packaging
