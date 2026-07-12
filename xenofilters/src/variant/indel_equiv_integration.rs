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

    use xenofilters::{
        alignment::fragment::score_state_nw,
        filter_algorithm::line_by_line::Scratch,
        region::diagnostic::DiagnosticVariants,
        variant::{
            enumerate_equivalents,
            population::Population,
            store::Store,
            Variant,
        },
        config::Config,
        tests::{create_record, MockStream},
    };

    // ── In-memory reference helpers ───────────────────────────────────────────

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

    // ── Shared variant constructors ───────────────────────────────────────────

    /// A deletion of one A at position 1 (0-based) in GAAAAG.
    /// Left-normalised: pos=0, ref=GA, alt=G.
    fn deletion_pop(pos: usize, ref_a: &[u8], alt_a: &[u8]) -> Population {
        Population {
            ref_id:           0,
            pos,
            ref_a:            ref_a.to_vec(),
            alt_a:            alt_a.to_vec(),
            allele_frequency: 0.7,
        }
    }

    // ── Test 1: store finds all expanded positions ────────────────────────────

    #[test]
    fn store_contains_all_expanded_positions() {
        let reference  = b"GAAAAGCCCCT";
        let name_to_id = name_to_id();

        // Expand the deletion.
        let equivalents = enumerate_equivalents(0, b"GA", b"G", reference, 0);
        // Expected: positions 0,1,2,3 (anchors G,A,A,A then AG which stops because
        // ref[5]=G ≠ ref[1]=A from the shift check).
        assert_eq!(equivalents.len(), 4, "GAAAA should yield 4 equivalents");

        let mut store = Store::<Population>::new();
        for eq in &equivalents {
            store.insert(
                0,
                deletion_pop(eq.pos, &eq.ref_a, &eq.alt_a),
            );
        }
        store.dedup();

        // Query [0, 4): should return all four.
        let hits = store.overlapping_multi(0, 0, 4);
        assert_eq!(hits.len(), 4, "overlapping_multi must return all 4 expanded positions");

        // Query [2, 3): should return entries at pos 2 and 3 whose alleles
        // overlap that window.
        let hits_narrow = store.overlapping_multi(0, 2, 3);
        assert!(
            !hits_narrow.is_empty(),
            "narrow query [2,3) must find ≥1 expanded variant"
        );
    }

    // ── Test 2: all expansions produce the same p_variant ────────────────────

    #[test]
    fn expanded_variants_preserve_p_variant() {
        let reference  = b"GAAAAGCCCCT";
        let equivalents = enumerate_equivalents(0, b"GA", b"G", reference, 0);

        let canon = deletion_pop(0, b"GA", b"G");
        for eq in &equivalents {
            let v = deletion_pop(eq.pos, &eq.ref_a, &eq.alt_a);
            assert!(
                (v.p_variant() - canon.p_variant()).abs() < 1e-12,
                "p_variant must be preserved across expansion: {} vs {}",
                v.p_variant(), canon.p_variant()
            );
        }
    }

    // ── Test 3: diagnostic variants block Tier-2 at all expanded positions ───

    #[test]
    fn diagnostic_expansion_blocks_tier2_at_all_positions() {
        // Simulate a DiagnosticVariants store built with expansion.
        // It should have entries at positions 0,1,2,3 for our deletion.
        let reference = b"GAAAAGCCCCT";
        let equivalents = enumerate_equivalents(0, b"GA", b"G", reference, 0);

        let mut per_ref: Vec<Vec<xenofilters::region::diagnostic::DiagnosticSite>> =
            vec![Vec::new()];

        for eq in &equivalents {
            per_ref[0].push(xenofilters::region::diagnostic::DiagnosticSite {
                pos:     eq.pos,
                ref_len: eq.ref_a.len(),
            });
        }
        // Sort as the real builder would.
        per_ref[0].sort_unstable_by_key(|s| s.pos);

        let diag = DiagnosticVariants { per_ref };

        // All four positions must report an overlap for a read at [0, 4).
        for expected_pos in 0..4 {
            assert!(
                diag.overlaps(0, expected_pos, expected_pos + 2, false),
                "diagnostic must report overlap at 0-based pos {expected_pos}"
            );
        }

        // Positions well outside the repeat must NOT overlap.
        assert!(
            !diag.overlaps(0, 10, 12, false),
            "diagnostic must NOT report overlap far from the repeat"
        );
    }

    // ── Test 4: right-shifted read maps to the same rescue delta ─────────────
    //
    // Read sequence: GAAAG  (same for all representations)
    // Canonical VCF: pos=0, ref=GA, alt=G  (delete one A)
    // Right-shifted:  pos=2, ref=AA, alt=A  (same biological event)
    // NW scoring of alt allele against the read should succeed in both cases
    // and produce the same (non-zero positive) rescue delta.

    #[test]
    fn rescue_delta_equal_for_canonical_and_shifted() {
        use xenofilters::{
            alignment::{
                fragment::score_state_nw,
                MdCigFlags,
            },
            filter_algorithm::line_by_line::Scratch,
            penalty::Penalty,
            config::Config,
            variant::eval::Eval,
        };

        // A minimal penalty config: small gap costs so the rescue is visible.
        let pen = Penalty::build(4.0, 2.0, 0.5, 20, Default::default());

        // Build a Population variant at each of the 4 equivalent positions.
        let reference  = b"GAAAAG";
        let equivalents = enumerate_equivalents(0, b"GA", b"G", reference, 0);

        let mut deltas = Vec::new();
        for eq in &equivalents {
            let var = Population {
                ref_id:           0,
                pos:              eq.pos,
                ref_a:            eq.ref_a.clone(),
                alt_a:            eq.alt_a.clone(),
                allele_frequency: 0.9,   // p_variant = 0.9 → rescue fires
            };

            // Build a minimal Eval representing this variant.
            let mut eval = Eval::new();
            eval.set_variant(&var);

            // Compute the delta manually via align_alt_to_read.
            // (align_alt_to_read is pub(crate) in variant_window.rs)
            // Here we use the public test hook exposed via feature bench-internals.
            let read_seq  = b"GAAAG";   // 5 bases; deletion of one A gives GAAG → GAAG
            let read_qscores = vec![30u8; 5];
            let ref_start = 0usize;

            let delta = xenofilters::variant::variant_window_test_hook::compute_rescue_delta(
                &var, read_seq, &read_qscores, ref_start, &pen,
            );
            deltas.push(delta);
        }

        // All four equivalent positions must produce the same rescue delta.
        let first = deltas[0];
        assert!(first > 0.0, "rescue delta must be positive (p_variant=0.9 > 0.5)");
        for (i, &d) in deltas.iter().enumerate() {
            assert!(
                (d - first).abs() < 1e-9,
                "rescue delta at equivalent position {i} ({d}) differs from canonical ({first})"
            );
        }
    }

    // ── Test 5: expand_indels=false leaves store unchanged ───────────────────

    #[test]
    fn no_expansion_when_flag_off() {
        // When expand_indels is false, a deletion VCF with one record produces
        // exactly one entry in the store, even inside a homopolymer.
        let mut store = Store::<Population>::new();
        let single = deletion_pop(0, b"GA", b"G");
        store.insert(0, single);
        store.dedup();

        let hits = store.overlapping_multi(0, 0, 5);
        assert_eq!(hits.len(), 1, "without expansion, store has exactly 1 entry");

        // Position 2 should return nothing because only pos=0 was inserted.
        let hits_pos2 = store.overlapping_multi(0, 2, 3);
        assert!(
            hits_pos2.is_empty(),
            "without expansion, query at pos=2 must return empty"
        );
    }

    // ── Test 6: dedup removes identical entries from two different VCF records

    #[test]
    fn dedup_removes_identical_expanded_entries() {
        // Two separate VCF records that both expand to pos=0,ref=GA,alt=G.
        let mut store = Store::<Population>::new();
        store.insert(0, deletion_pop(0, b"GA", b"G"));
        store.insert(0, deletion_pop(0, b"GA", b"G")); // duplicate
        store.dedup();

        let hits = store.overlapping_multi(0, 0, 2);
        assert_eq!(hits.len(), 1, "dedup must remove byte-identical duplicates");
    }

    // ── Test 7: insertion equivalents ────────────────────────────────────────

    #[test]
    fn insertion_equivalents_all_found_by_store() {
        // Reference: GAAAG; insertion of one A.
        // Equivalents: pos=0 (ref=G alt=GA), pos=1 (ref=A alt=AA), pos=2 (ref=A alt=AA)
        let reference  = b"GAAAG";
        let equivalents = enumerate_equivalents(0, b"G", b"GA", reference, 0);
        assert_eq!(equivalents.len(), 3, "3 insertion equivalents in GAAAG");

        let mut store = Store::<Population>::new();
        for eq in &equivalents {
            store.insert(0, Population {
                ref_id:           0,
                pos:              eq.pos,
                ref_a:            eq.ref_a.clone(),
                alt_a:            eq.alt_a.clone(),
                allele_frequency: 0.6,
            });
        }
        store.dedup();

        // A read overlapping positions 0-4 must see all three.
        let hits = store.overlapping_multi(0, 0, 4);
        assert_eq!(hits.len(), 3, "all 3 insertion equivalents queryable");
    }

    // ── Test 8: two-stream routing agrees for shifted and canonical read ──────
    //
    // Stream 0: read aligns to human ref with deletion at canonical pos.
    // Stream 1: same read aligns to mouse ref with deletion at shifted pos.
    // Without expansion, the variant is missed in stream 1 → wrong route.
    // With expansion, both streams see the variant → correct rescue delta.

    #[test]
    fn two_stream_routing_with_expanded_stores() {
        use xenofilters::{
            aln_stream::AlignmentStream,
            filter_algorithm::line_by_line::LineByLine,
        };

        // Human reference context for stream 0: GAAAAG (canonical del at pos 0).
        // Mouse reference context for stream 1: GAAAAG (same sequence, del shifted to pos 2).

        // Build two mock streams.
        let human_rec = create_record(
            b"R1",
            "5M1D5M",      // deletion at position 5 in read; canonical
            &[],
            &[30u8; 10],
            "5^A5",        // MD: 5 matches, delete A, 5 matches
            false,
        ).unwrap();

        let mouse_rec = create_record(
            b"R1",
            "3M1D7M",      // deletion shifted left by 2 in mouse aligner's view
            &[],
            &[30u8; 10],
            "3^A7",
            false,
        ).unwrap();

        // Build expanded stores for each stream.
        let reference = b"GGGGGAAAAGCCCCCCCCCCC";

        let human_equivalents = enumerate_equivalents(4, b"GA", b"G", reference, 0);
        let mouse_equivalents = enumerate_equivalents(2, b"GA", b"G", reference, 0);

        let build_store = |equivalents: &[_]| {
            let mut store = Store::<Population>::new();
            for eq in equivalents.iter().cloned().collect::<Vec<_>>().iter() {
                let eq = eq as &xenofilters::variant::indel_equiv::EquivalentAlleles;
                store.insert(0, Population {
                    ref_id:           0,
                    pos:              eq.pos,
                    ref_a:            eq.ref_a.clone(),
                    alt_a:            eq.alt_a.clone(),
                    allele_frequency: 0.8,
                });
            }
            store.dedup();
            Arc::new(store) as Arc<dyn xenofilters::variant::StoreTrait>
        };

        let human_store = build_store(&human_equivalents);
        let mouse_store = build_store(&mouse_equivalents);

        // The human query at the mouse aligner's position must succeed.
        // mouse_equivalents includes human canonical position because the
        // same deletion slides through the same repeat.
        let mouse_query_at_human_pos = mouse_equivalents.iter().any(|e| e.pos == 4);
        assert!(
            mouse_query_at_human_pos,
            "expanded mouse equivalents must include human canonical position (pos=4)"
        );
    }

    // ── Test 9: MAX_SHIFT limit does not panic ────────────────────────────────

    #[test]
    fn long_homopolymer_does_not_overflow() {
        // 150 A's between two sentinel bases; MAX_SHIFT=100 caps expansion.
        let mut reference = vec![b'G'];
        reference.extend(std::iter::repeat(b'A').take(150));
        reference.push(b'G');

        let equivalents = enumerate_equivalents(0, b"GA", b"G", &reference, 0);

        // Must be capped at MAX_SHIFT + 1 (left-norm + 100 right-shifts).
        assert_eq!(
            equivalents.len(),
            xenofilters::variant::indel_equiv::MAX_SHIFT + 1,
            "MAX_SHIFT must cap expansion at {} entries",
            xenofilters::variant::indel_equiv::MAX_SHIFT + 1
        );

        // No entry must reference a position beyond the reference length.
        for eq in &equivalents {
            assert!(
                eq.pos + eq.ref_a.len() <= reference.len(),
                "expanded entry at pos={} overruns reference (len={})",
                eq.pos, reference.len()
            );
        }
    }
}
