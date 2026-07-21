# Indel Equivalence Expansion — Implementation Checklist

All files listed below require changes.  Files marked **NEW** are fully
provided in the feature branch; files marked **PATCH** require the targeted
merge described in the inline comments of each file.

---

## New files (already written in full)

| File | Purpose |
|------|---------|
| `src/variant/indel_equiv.rs`           | Core enumeration algorithm + build functions |
| `src/variant/indel_equiv_impls.rs`     | `WithAllelesRefId` impls on Sample + Population |
| `src/variant/indel_equiv_corrected.rs` | Corrected build functions with ref_id threading |
| `src/variant/diagnostic_equiv.rs`      | SegregateVariants expansion |
| `src/variant/store_insert.rs`          | `Store::insert`, `insert_expanded`, `dedup` |
| `src/variant/variant_window_test_hook.rs` | Test-only scoring hook |
| `src/variant/name_to_id.rs`           | Header/FAI name→id mapping |
| `src/integration.rs`                  | `build_variant_stores`, BED padding wiring |
| `tests/indel_equiv_integration.rs`    | End-to-end integration tests |
| `src/variant/indel_equiv/tests.rs`    | Unit + property tests |

---

## Patches to existing files

### `src/variant/mod.rs`
- Add `pub(crate) mod indel_equiv;`
- Add `pub(crate) mod indel_equiv_impls;`
- Add `pub(crate) mod diagnostic_equiv;`
- Add `pub(crate) mod name_to_id;`
- Add `fn ref_id(&self) -> usize;` to the `Variant` trait

### `src/variant/sample.rs`
```rust
// ADD field to Sample struct:
pub(crate) ref_id: usize,

// ADD to Variant impl:
fn ref_id(&self) -> usize { self.ref_id }

// ADD derive:
#[derive(Clone)]
```

### `src/variant/population.rs`
```rust
// ADD field to Population struct:
pub(crate) ref_id: usize,

// ADD to Variant impl:
fn ref_id(&self) -> usize { self.ref_id }

// ADD derive:
#[derive(Clone)]
```

### `src/variant/store.rs`
- Add `pub(crate) per_ref: HashMap<usize, Vec<V>, RandomState>` if not present
- Replace or supplement `overlapping_multi` with the implementation in
  `store_and_diagnostic.rs` (binary search + scan-back)
- Add `insert`, `insert_expanded`, `dedup` from `store_insert.rs`

### `src/region/diagnostic.rs`
- Replace `SegregateVariants` with the version in `store_and_diagnostic.rs`
  (adds `is_reverse` parameter to `overlaps()`, `from_per_ref` constructor)
- Add `#[derive(Clone, Default)]` to `SegregateVariants`

### `src/config.rs`
- Add `expand_indels: bool` field to `CommonArgs`
- Add `indel_expand_padding: usize` field to `CommonArgs`
- Call `validate_indel_expansion(&common)` from `validate_common()`

### `src/error.rs`
- Add `ExpandIndelsRequiresReference`
- Add `FastaIndexMissing(PathBuf)`
- Add `FastaFetch { region: String, source: std::io::Error }`
- Add `ExpansionPositionOverflow { chrom: String, pos: usize }`

### `src/aln_stream.rs`
- Replace inline sample/population store construction with
  `build_variant_stores(opt, &header, i)?` from `src/integration.rs`
- For diagnostic stores (hashlookup / collated path):
  call `build_diagnostic_store_for_stream(opt, &header, i)?`

### `src/main.rs`
- In `run_hashlookup`: replace `load_ambiguous_regions` with
  `load_ambiguous_regions_with_padding(&args, &name_to_id, padding)`
- Pass `vcf` array built from `build_diagnostic_store_for_stream` to
  `HashLookup::new_from_hashlookup`

### `src/variant/variant_window.rs`
- Make `align_alt_to_read` `pub(crate)` if not already
  (required by `variant_window_test_hook.rs`)

### `Cargo.toml`
- Confirm `noodles` has `fasta` feature enabled (already listed)
- Add `proptest = "1"` to `[dev-dependencies]` if not present
- Ensure `[features] bench-internals = []` is declared

### `ROADMAP.md`
- Add the four deferred items from `store_and_diagnostic.rs` (MNP expansion,
  BED padding, MAX_SCAN_BACK, strand-aware diagnostic overlap)

---

## Coordinate contract (enforced throughout)

| Boundary | Convention |
|----------|-----------|
| VCF input `variant_start()` | 1-based, inclusive → subtract 1 on ingestion |
| All `Store<V>`, `DiagnosticSite`, `EquivalentAlleles` | 0-based |
| `overlapping_multi(tid, start, end)` query | 0-based, half-open `[start, end)` |
| FASTA fetch `Position::new(x)` | 1-based (`x = ctx_start + 1`) |
| BED intervals in `ScoredRegions` | 0-based, half-open (already) |

---

## Build verification sequence

```bash
# 1. Compile (no warnings expected)
cargo check --all-features 2>&1 | grep -v "^warning\[unused\]"

# 2. Unit tests (fast; no FASTA I/O)
cargo test -p xenofilters variant::indel_equiv

# 3. Integration tests (require --features bench-internals)
cargo test --features bench-internals --test indel_equiv_integration

# 4. Property tests
cargo test --features bench-internals variant::indel_equiv::tests::property

# 5. Full test suite
cargo test --all-features

# 6. Clippy
cargo clippy --all-features -- -D warnings
```

---

## Known deferred items (see ROADMAP.md)

1. MNP / complex allele expansion
2. Sequence-aware BED padding (replace fixed-bp with repeat-boundary detection)
3. `overlapping_multi` MAX_SCAN_BACK → interval-tree for SV VCFs
4. Strand-aware SegregateVariants overlap
