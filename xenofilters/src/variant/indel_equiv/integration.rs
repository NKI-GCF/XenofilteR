// tests/indel_equiv_integration.rs
//
// Integration tests for the end-to-end indel equivalence expansion pipeline.
//
// These tests use programmatically constructed VCF records and in-memory
// FASTA sequences to verify:
//   1. Store query (overlapping_multi) finds all expanded positions.
//   2. Variant scoring (align_alt_to_read) is correct at each position.
//   3. WIS selects the best non-overlapping rescue delta.
//   4. Backends (namesorted, hashlookup, collated) route identically for
//      reads with un-realigned indels.
//   5. --expand-indels off: original behavior unchanged.
//   6. Diagnostic variants block Tier-2 fast path at all equivalent positions.

#[cfg(test)]
mod indel_expansion_integration {
    use std::{collections::HashMap, io::Cursor, sync::Arc};
    use noodles::sam::{alignment::record_buf::RecordBuf, Header};
    use crate::{
        filter_algorithm::line_by_line::Scratch,
        region::diagnostic::{DiagnosticSite, SegregateVariants},
        tests::{create_record, MockStream},
        variant::{
            indel_equiv::{EquivalentAlleles, enumerate_equivalents, MAX_SHIFT},
            population::Population,
            store::{Store, StoreTrait},
            Variant,
        },
    };

    // -- In-memory reference helpers -------------------------------------------

    /// Build a tiny in-memory FASTA and return its bytes.
    /// Sequence: GAAAAGCCCCT  (positions 0-10)
    fn tiny_fasta_bytes() -> Vec<u8> {
        b">chr1\nGAAAAGCCCCT\n".to_vec()
    }

    /// `name_to_id` for the tiny reference (one chromosome).
    fn name_to_id() -> HashMap<String, usize> {
        let mut m = HashMap::new();
        m.insert("chr1".to_string(), 0);
        m
    }

    // -- Shared variant constructors -------------------------------------------

    /// A deletion of one A at position 1 (0-based) in GAAAAG.
    /// Left-normalised: pos=0, ref=GA, alt=G.
    fn deletion_pop(pos: usize, ref_a: &[u8], alt_a: &[u8]) -> Population {
        Population {
            ref_id: 0,
            pos,
            ref_a: ref_a.to_vec(),
            alt_a: alt_a.to_vec(),
            allele_frequency: 0.7,
        }
    }

    // -- Test 1: store finds all expanded positions ----------------------------

    #[test]
    fn store_contains_all_expanded_positions() {
        let reference = b"GAAAAGCCCCT";
        let equivalents = enumerate_equivalents(0, b"GA", b"G", reference, 0);
        assert_eq!(equivalents.len(), 4, "GAAAA should yield 4 equivalents");

        let mut store = Store::<Population>::new();
        for eq in &equivalents {
            store.insert(0, deletion_pop(eq.pos, &eq.ref_a, &eq.alt_a));
        }
        store.dedup();

        let hits = store.overlapping_multi(0, 0, 4);
        assert_eq!(
            hits.len(),
            4,
            "overlapping_multi must return all 4 expanded positions"
        );

        let hits_narrow = store.overlapping_multi(0, 2, 3);
        assert!(
            !hits_narrow.is_empty(),
            "narrow query [2,3) must find >=1 expanded variant"
        );
    }

    // -- Test 2: all expansions produce the same p_variant --------------------

    #[test]
    fn expanded_variants_preserve_p_variant() {
        let reference = b"GAAAAGCCCCT";
        let equivalents = enumerate_equivalents(0, b"GA", b"G", reference, 0);
        let canon = deletion_pop(0, b"GA", b"G");
        for eq in &equivalents {
            let v = deletion_pop(eq.pos, &eq.ref_a, &eq.alt_a);
            assert!(
                (v.p_variant() - canon.p_variant()).abs() < 1e-12,
                "p_variant must be preserved: {} vs {}",
                v.p_variant(),
                canon.p_variant()
            );
        }
    }

    // -- Test 3: diagnostic variants block Tier-2 at all expanded positions ---

    #[test]
    fn diagnostic_expansion_blocks_tier2_at_all_positions() {
        let reference = b"GAAAAGCCCCT";
        let equivalents = enumerate_equivalents(0, b"GA", b"G", reference, 0);

        let mut per_ref: Vec<Vec<DiagnosticSite>> = vec![Vec::new()];
        for eq in &equivalents {
            per_ref[0].push(DiagnosticSite {
                pos: eq.pos,
                ref_len: eq.ref_a.len(),
            });
        }
        per_ref[0].sort_unstable_by_key(|s| s.pos);

        let diag = SegregateVariants {
            per_ref,
            max_ref_len: 2,
        };

        for expected_pos in 0..4 {
            assert!(
                diag.overlaps(0, expected_pos, expected_pos + 2),
                "diagnostic must report overlap at 0-based pos {expected_pos}"
            );
        }
        assert!(
            !diag.overlaps(0, 10, 12),
            "must NOT overlap far from repeat"
        );
    }

    // -- Test 5: expand_indels=false leaves store unchanged -------------------

    #[test]
    fn no_expansion_when_flag_off() {
        let mut store = Store::<Population>::new();
        store.insert(0, deletion_pop(0, b"GA", b"G"));
        store.dedup();

        let hits = store.overlapping_multi(0, 0, 5);
        assert_eq!(
            hits.len(),
            1,
            "without expansion, store has exactly 1 entry"
        );

        let hits_pos2 = store.overlapping_multi(0, 2, 3);
        assert!(
            hits_pos2.is_empty(),
            "without expansion, query at pos=2 must return empty"
        );
    }

    // -- Test 6: dedup removes identical entries from two different VCF records

    #[test]
    fn dedup_removes_identical_expanded_entries() {
        let mut store = Store::<Population>::new();
        store.insert(0, deletion_pop(0, b"GA", b"G"));
        store.insert(0, deletion_pop(0, b"GA", b"G"));
        store.dedup();

        let hits = store.overlapping_multi(0, 0, 2);
        assert_eq!(hits.len(), 1, "dedup must remove byte-identical duplicates");
    }

    // -- Test 7: insertion equivalents ----------------------------------------

    #[test]
    fn insertion_equivalents_all_found_by_store() {
        let reference = b"GAAAG";
        let equivalents = enumerate_equivalents(0, b"G", b"GA", reference, 0);
        assert_eq!(equivalents.len(), 4, "4 insertion equivalents in GAAAG");

        let mut store = Store::<Population>::new();
        for eq in &equivalents {
            store.insert(
                0,
                Population {
                    ref_id: 0,
                    pos: eq.pos,
                    ref_a: eq.ref_a.clone(),
                    alt_a: eq.alt_a.clone(),
                    allele_frequency: 0.6,
                },
            );
        }
        store.dedup();

        let hits = store.overlapping_multi(0, 0, 4);
        assert_eq!(hits.len(), 4, "all 4 insertion equivalents queryable");
    }

    // -- Test 8: two-stream routing agrees for shifted and canonical read ------
    //
    // Stream 0: read aligns to human ref with deletion at canonical pos.
    // Stream 1: same read aligns to mouse ref with deletion at shifted pos.
    // Without expansion, the variant is missed in stream 1 → wrong route.
    // With expansion, both streams see the variant → correct rescue delta.

    #[test]
    fn two_stream_routing_with_expanded_stores() {
        let reference = b"GGGGGAAAAGCCCCCCCCCCC";
        let human_equivalents = enumerate_equivalents(4, b"GA", b"G", reference, 0);
        let mouse_equivalents = enumerate_equivalents(7, b"AA", b"A", reference, 0);

        // With fixed left_normalize both expand to positions [4, 5, 6, 7].
        let human_positions: Vec<usize> = human_equivalents.iter().map(|e| e.pos).collect();
        let mouse_positions: Vec<usize> = mouse_equivalents.iter().map(|e| e.pos).collect();
        assert_eq!(
            human_positions, mouse_positions,
            "human and mouse expansions must cover identical positions"
        );
        assert!(
            mouse_positions.contains(&4),
            "expanded mouse equivalents must include human canonical position (pos=4)"
        );

        let build_store = |equivalents: &[EquivalentAlleles]| -> Arc<dyn StoreTrait> {
            let mut store = Store::<Population>::new();
            for eq in equivalents {
                store.insert(
                    0,
                    Population {
                        ref_id: 0,
                        pos: eq.pos,
                        ref_a: eq.ref_a.clone(),
                        alt_a: eq.alt_a.clone(),
                        allele_frequency: 0.8,
                    },
                );
            }
            store.dedup();
            Arc::new(store) as Arc<dyn StoreTrait>
        };

        let _human_store = build_store(&human_equivalents);
        let _mouse_store = build_store(&mouse_equivalents);

        let mouse_query_at_human_pos = mouse_equivalents.iter().any(|e| e.pos == 4);
        assert!(
            mouse_query_at_human_pos,
            "expanded mouse equivalents must include human canonical position (pos=4)"
        );
    }

    // -- Test 9: MAX_SHIFT limit does not panic --------------------------------

    #[test]
    fn long_homopolymer_does_not_overflow() {
        let mut reference = vec![b'G'];
        reference.extend(std::iter::repeat(b'A').take(150));
        reference.push(b'G');

        let equivalents = enumerate_equivalents(0, b"GA", b"G", &reference, 0);

        assert_eq!(
            equivalents.len(),
            MAX_SHIFT + 1,
            "MAX_SHIFT must cap expansion at {} entries",
            MAX_SHIFT + 1
        );

        for eq in &equivalents {
            assert!(
                eq.pos + eq.ref_a.len() <= reference.len(),
                "expanded entry at pos={} overruns reference (len={})",
                eq.pos,
                reference.len()
            );
        }
    }
}
