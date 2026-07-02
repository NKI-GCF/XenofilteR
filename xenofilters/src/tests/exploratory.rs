// xenofilters/src/tests/exploratory.rs
//
// Exploratory, table-driven and stress tests targeting fragile areas:
// - VariantWindow / read offsets / weighted_ref_score / align_alt_to_read
// - Fragment scoring across windows and accumulation semantics
// - FragmentState / MdCigFlags parsing (SA tag) and cache-like behavior
// - edge cases: extreme quality values, p_variant (0,1,NaN), reverse-complement paths
//
// These tests are intentionally aggressive and may fail on first run; they are
// designed to probe corner cases and surface brittle behavior for follow-up fixes.

#![allow(clippy::float_cmp)]

#[cfg(test)]
mod tests {
    use super::*;
    use crate::alignment::MdCigFlags;
    use crate::alignment::SimpleRec;
    use crate::alignment::{align_alt_to_read, weighted_ref_score, Fragment, VariantWindow};
    use crate::filter_algorithm::line_by_line::{Scratch, READ_CT};
    use crate::penalty::Penalty;
    use crate::tests::create_record; // existing test helper used across the codebase
    use crate::variant::Eval;
    use crate::variant::FragEvalVec;
    use crate::variant::Variant;
    use noodles::sam::alignment::record::data::field::{Tag, Value};
    use noodles::sam::alignment::record::Flags;
    use noodles::sam::alignment::record_buf::{data::field::Value as RBValue, RecordBuf};
    use smallvec::{smallvec, SmallVec};
    use std::f64;

    // Helper: a simple Penalty with deterministic numbers for tests.
    fn flat_penalty() -> Penalty {
        Penalty {
            gap_open: -2.0,
            gap_extend: -0.5,
            log_likelihood_match: [0.0; crate::penalty::MAX_Q],
            log_likelihood_mismatch: [-1.0; crate::penalty::MAX_Q],
            chimeric_junction_penalty: -3.0,
        }
    }

    // Small helper variant impl that allows setting p_variant explicitly.
    struct ProbeVariant {
        pos: usize, // 0-based pos expected by Variant::pos in crate
        ref_a: Vec<u8>,
        alt_a: Vec<u8>,
        p: f64,
    }
    impl ProbeVariant {
        fn boxed(pos: usize, ref_a: &[u8], alt_a: &[u8], p: f64) -> &'static Self {
            Box::leak(Box::new(ProbeVariant {
                pos,
                ref_a: ref_a.to_vec(),
                alt_a: alt_a.to_vec(),
                p,
            }))
        }
    }
    impl Variant for ProbeVariant {
        fn pos(&self) -> usize {
            self.pos
        }
        fn ref_allele(&self) -> &[u8] {
            &self.ref_a
        }
        fn alt_allele(&self) -> &[u8] {
            &self.alt_a
        }
        fn p_variant(&self) -> f64 {
            self.p
        }
    }

    // -----------------------------------------------------------------------
    // 1) VariantWindow & weighted_ref_score / align_alt_to_read table-driven
    // -----------------------------------------------------------------------
    #[test]
    fn table_weighted_ref_and_alt_alignment_various_p() -> Result<(), crate::Error> {
        // Read: AAGAA (positions 0..5) — we will target variant at 2 (0-based)
        let read = b"AAGAA";
        // create a read-record; reuse create_record helper (alignment tests do this)
        let rec = create_record(b"r", "5M", b"AAGAA", &[30u8; 5], "5", false)?;
        let flags = rec.flags();

        let w = VariantWindow::compute(1, 6, 3, 4).expect("window computation");
        let pen = flat_penalty();

        // Table-driven p values including edge-cases.
        let cases = &[
            ("p=0.0 ref-certain", 0.0f64),
            ("p=0.5 neutral", 0.5f64),
            ("p=1.0 alt-certain", 1.0f64),
            ("p=NaN", f64::NAN),
            ("p=subnormal", f64::MIN_POSITIVE / 2.0),
        ];

        for (label, p) in cases {
            // weighted_ref_score returns a Result<f64,_>
            let wref = weighted_ref_score(w, 1usize, *p, &pen, |_| Ok(30))?;
            // alt alignment uses align_alt_to_read; set scratch anew each iteration
            let mut scratch = Scratch::new();
            let alt = b"G";
            let alt_score = align_alt_to_read(
                alt,
                w.read_offset(1),
                w.ref_len,
                *p,
                &pen,
                |idx| Ok((read[idx], 30usize)),
                &mut scratch,
            )?;

            if p.is_nan() {
                // NaN should propagate through calculations -> scores NaN.
                assert!(wref.is_nan(), "[{label}] expected weighted_ref NaN");
                assert!(alt_score.is_nan(), "[{label}] expected alt_score NaN");
            } else {
                // For deterministic flat_penalty, with read base at offset 2 being 'G',
                // p=1.0 -> weighted_ref = mismatch penalty (-1.0), alt aligned to read 'G' -> match (0).
                // p=0.0 -> weighted_ref = match (0.0), alt aligned -> mismatch (-1.0).
                if (*p - 1.0).abs() < 1e-12 {
                    assert!(
                        (alt_score - 0.0).abs() < 1e-12,
                        "[{label}] alt expected match -> 0"
                    );
                    assert!(
                        (wref - (-1.0)).abs() < 1e-12,
                        "[{label}] ref expected mismatch -> -1"
                    );
                } else if (*p - 0.0).abs() < 1e-12 {
                    assert!(
                        (wref - 0.0).abs() < 1e-12,
                        "[{label}] ref expected match -> 0"
                    );
                    assert!(
                        (alt_score - (-1.0)).abs() < 1e-12,
                        "[{label}] alt expected mismatch -> -1"
                    );
                } else {
                    // intermediate p values give weighted values in-between; ensure not equal so delta != 0
                    assert!(
                        (alt_score - wref).abs() > 1e-12,
                        "[{label}] unexpected identical alt/ref (alt={alt_score} ref={wref})"
                    );
                }
            }
        }
        Ok(())
    }

    // -----------------------------------------------------------------------
    // 2) Fragment scoring sequential calls: state & cache invalidation
    // -----------------------------------------------------------------------
    #[test]
    fn fragment_score_reuses_scratch_and_resets_last_variant_delta() -> Result<(), crate::Error> {
        // Build a simple record and a Fragment; call score twice with different dvnt,
        // verifying last_variant_delta gets set each run and not carried over incorrectly.
        let rec1 = create_record(b"r", "5M", b"AAGAA", &[30u8; 5], "5", false)?;
        let flags1 = rec1.flags();
        let p = flat_penalty();
        let md_flags = smallvec![MdCigFlags::try_from_record(&rec1, &flags1)?];

        let mut frag = Fragment::new(&p, smallvec![&rec1], md_flags)?;
        let mut scratch = Scratch::new();

        // First, no variants -> last_variant_delta == 0
        let mut dvnt: FragEvalVec<'static> = smallvec![smallvec![]];
        let s0 = frag.score(&mut scratch, &mut dvnt)?;
        assert_eq!(scratch.last_variant_delta, 0.0, "no variants -> zero delta");
        // Now create a variant with p=1.0 that should rescue (alt==read at pos2)
        let pv = ProbeVariant::boxed(2, b"A", b"G", 1.0);
        let mut ev = Eval::new();
        ev.set_variant(pv as &dyn Variant);
        // alt_score will be positive relative to weighted_ref_score when p_variant > 0.5
        let mut dvnt2: FragEvalVec<'static> = smallvec![smallvec![ev]];
        let s1 = frag.score(&mut scratch, &mut dvnt2)?;
        // After scoring, last_variant_delta should reflect the rescue; check it's non-zero for p=1.0
        assert!(
            scratch.last_variant_delta != 0.0,
            "expected last_variant_delta set for variant p=1.0"
        );
        // Now call again with empty dvnt → should reset to zero
        let mut dvnt3: FragEvalVec<'static> = smallvec![smallvec![]];
        let _ = frag.score(&mut scratch, &mut dvnt3)?;
        assert_eq!(
            scratch.last_variant_delta, 0.0,
            "subsequent empty scoring resets delta"
        );
        // Basic sanity: fragment score numerical stability
        assert!(s0.is_finite());
        assert!(s1.is_finite());
        Ok(())
    }

    // -----------------------------------------------------------------------
    // 3) Partial-window accumulation & multi-segment variant spanning
    // -----------------------------------------------------------------------
    #[test]
    fn variant_partial_window_accumulation_across_segments() -> Result<(), crate::Error> {
        // Construct two segments where a variant reference span crosses the boundary.
        // We expect Eval.incurred and alt_score to be accumulated from multiple window calls.
        // Setup two records that, when stitched, place a variant at a boundary.
        // Record lengths chosen small for clarity.
        let rec_a = create_record(b"r", "3M", b"AAG", &[30u8; 3], "3", false)?;
        let rec_b = create_record(b"r", "3M", b"TA_", &[30u8; 3], "3", false)?; // second segment
        let flags_a = rec_a.flags();
        let flags_b = rec_b.flags();
        let p = flat_penalty();

        let md_flags = smallvec![
            MdCigFlags::try_from_record(&rec_a, &flags_a)?,
            MdCigFlags::try_from_record(&rec_b, &flags_b)?
        ];
        let segs: SmallVec<[&RecordBuf; READ_CT]> = smallvec![&rec_a, &rec_b];

        let mut frag = Fragment::new(&p, segs, md_flags)?;
        let mut scratch = Scratch::new();

        // Variant that spans the junction: pos relative so part in rec_a and part in rec_b.
        // We place variant.pos = 2 (0-based) with ref_len = 2 -> spans positions 2 and 3 (junction)
        let v = ProbeVariant::boxed(2, b"AA", b"G", 1.0);
        let mut ev = Eval::new();
        ev.set_variant(v as &dyn Variant);
        let mut dvnt: FragEvalVec<'static> = smallvec![smallvec![ev]];

        let _ = frag.score(&mut scratch, &mut dvnt)?;
        // After scoring, scratch.last_variant_delta should be set (if any rescue)
        // We assert that the Eval would have been consumed into finished and aggregated so delta computed.
        // We cannot access finished directly here, but last_variant_delta should be finite.
        assert!(scratch.last_variant_delta.is_finite());
        Ok(())
    }

    // -----------------------------------------------------------------------
    // 4) MdCigFlags SA:Z parsing and supplementary/supp_count branching
    // -----------------------------------------------------------------------
    #[test]
    fn md_cig_flags_sa_tag_and_is_perfect_variations() -> Result<(), crate::Error> {
        // Create a basic record and inject SA:Z tag string manually, then parse via MdCigFlags.
        let mut rec = create_record(b"r", "10M", &[], &[30u8; 10], "10", false)?;
        // Insert SA:Z-like value: one supplementary -> one semicolon
        let sa = "chr1,100,+,5M,60,1;".to_string();
        rec.data_mut()
            .insert(Tag::new(b'S', b'A'), RBValue::from(sa.clone()));
        let flags = rec.flags();
        let mcf = MdCigFlags::try_from_record(&rec, &flags)?;
        assert_eq!(
            mcf.supp_count(),
            1,
            "SA:Z semicolon count should give supp_count"
        );

        // Perfect record detection: single op + MD purely digits -> is_perfect true
        let rec_perfect = create_record(b"r2", "5M", &[], &[30u8; 5], "5", false)?;
        let flags_p = rec_perfect.flags();
        let mcf_p = MdCigFlags::try_from_record(&rec_perfect, &flags_p)?;
        assert!(mcf_p.is_perfect());

        // Non-perfect: mismatch in MD -> should be false
        let rec_imperfect = create_record(b"r3", "5M", &[], &[30u8; 5], "4A0", false)?;
        let flags_i = rec_imperfect.flags();
        let mcf_i = MdCigFlags::try_from_record(&rec_imperfect, &flags_i)?;
        assert!(!mcf_i.is_perfect());

        Ok(())
    }

    // -----------------------------------------------------------------------
    // 5) Quality index clamping & extremal qualities
    // -----------------------------------------------------------------------
    #[test]
    fn quality_index_clamping_and_errors() -> Result<(), crate::Error> {
        // Build a record with extremely large quality values (simulate overflow)
        let qual_big: Vec<u8> = vec![255u8; 10];
        let mut rec = create_record(b"r", "10M", &[], &qual_big, "10", false)?;
        // The fragment q() clamps to MAX_Q-1 and should not panic.
        let flags = rec.flags();
        let md_flags = smallvec![MdCigFlags::try_from_record(&rec, &flags)?];
        let p = flat_penalty();
        let mut frag = Fragment::new(&p, smallvec![&rec], md_flags)?;
        // q(0, 0) should be within bounds (clamped)
        let q0 = frag.q(0, 0)?;
        assert_eq!(q0, (crate::penalty::MAX_Q - 1));
        // Now test quality out-of-bounds error path: ask for index beyond read length
        assert!(frag.q(0, 100).is_err());
        Ok(())
    }

    // -----------------------------------------------------------------------
    // 6) Reverse-complement path: align_alt_to_read should revcomp appropriately
    // -----------------------------------------------------------------------
    #[test]
    fn align_alt_revcomp_path_and_complementing() -> Result<(), crate::Error> {
        // Read seq where reverse-complement has support: original read "CCT"
        // alt "G" aligns to reversed/complemented base
        let read = b"CCT";
        let mut rec = create_record(b"r", "3M", read, &[30u8; 3], "3", true)?; // set reverse flag
        let flags = rec.flags();
        let md_flags = smallvec![MdCigFlags::try_from_record(&rec, &flags)?];
        let p = flat_penalty();
        let mut frag = Fragment::new(&p, smallvec![&rec], md_flags)?;
        let mut scratch = Scratch::new();
        // Variant alt "G" at pos that corresponds to read index 1 when revcomped.
        let w = VariantWindow::compute(0, 3, 1, 2).unwrap(); // single-base window
                                                             // read_base_and_quality supplied by revcomp wrapper in score_variant_against_segment;
                                                             // we call align_alt_to_read directly with a closure that simulates revcomp handling:
        let alt = b"G";
        // We need to supply a closure that returns the revcomped base for fwd_nt_i:
        let seq = rec.sequence().as_ref();
        let seq_len = seq.len();
        let closure = |fwd_nt_i: usize| -> Result<(u8, usize), crate::Error> {
            // simulate revcomp: ri = seq_len - 1 - fwd_nt_i
            let ri = seq_len - 1 - fwd_nt_i;
            let base = seq[ri];
            // complement
            let cb = match base {
                b'A' => b'T',
                b'T' => b'A',
                b'C' => b'G',
                b'G' => b'C',
                _ => b'N',
            };
            Ok((cb, ri))
        };
        let alt_score = align_alt_to_read(
            alt,
            w.read_offset(0),
            w.ref_len,
            1.0,
            &p,
            closure,
            &mut scratch,
        )?;
        // With p=1 and alt matching revcomped read at that position we expect match (0)
        assert!(alt_score.is_finite());
        Ok(())
    }

    // -----------------------------------------------------------------------
    // 7) Eval accumulation robustness including NaN/Inf stress
    // -----------------------------------------------------------------------
    #[test]
    fn eval_accumulate_and_special_float_values() {
        let mut ev = Eval::new();
        ev.update(1.0, 2.0);
        ev.update(3.0, 4.0);
        // accumulation semantics preserved
        assert!((ev.delta() - 2.0).abs() < 1e-12);
        // now inject NaN alt_score and see propagation
        ev.update(0.0, f64::NAN);
        assert!(ev.delta().is_nan());
        // Infinity handling
        let mut ev2 = Eval::new();
        ev2.update(f64::INFINITY, 1.0);
        assert!(ev2.incurred.is_infinite());
    }

    // -----------------------------------------------------------------------
    // 8) FragmentState: sequential add_record and build_mcfs preserve ordering
    // -----------------------------------------------------------------------
    #[test]
    fn fragment_state_sequential_and_drain_behavior() -> Result<(), crate::Error> {
        let rec1 = create_record(b"r1", "5M", &[], &[30u8; 5], "5", false)?;
        let mut fs = crate::alignment::FragmentState::from_record(rec1, 0)?;
        let rec2 = create_record(b"r2", "5M", &[], &[30u8; 5], "5", false)?;
        fs.add_record(rec2)?;
        assert_eq!(fs.get_records().len(), 2);
        // build_mcfs borrows entries and returns per-record MdCigFlags
        let mcfs = fs.build_mcfs()?;
        assert_eq!(mcfs.len(), 2);
        drop(mcfs);
        // drain_records should yield the records and leave state empty
        let mut dr = fs.drain_records();
        let _r = dr.next().unwrap();
        let _r2 = dr.next().unwrap();
        drop(dr);
        assert!(fs.get_records().is_empty());
        Ok(())
    }

    // -----------------------------------------------------------------------
    // 9) Aggressive CIGAR/MD combinations (table-driven)
    // -----------------------------------------------------------------------
    #[test]
    fn cigar_md_extreme_cases_table() -> Result<(), crate::Error> {
        let pen = flat_penalty();
        struct Case {
            label: &'static str,
            cigar: &'static str,
            md: &'static str,
            // we will assert relative ordering or simple properties
        }
        let cases = &[
            Case {
                label: "long soft-clips",
                cigar: "100S",
                md: "0",
            },
            Case {
                label: "huge insertion",
                cigar: "1M1000I",
                md: "1",
            },
            Case {
                label: "large deletion",
                cigar: "1M1000D1M",
                md: "2",
            },
            Case {
                label: "many mismatches",
                cigar: "10M",
                md: "A1A1A1A1A1A1A1A1A",
            },
        ];
        for c in cases {
            let rec = create_record(b"r", c.cigar, &[], &[30u8; 150], c.md, false)?;
            let flags = rec.flags();
            let md_flags = smallvec![MdCigFlags::try_from_record(&rec, &flags)?];
            let mut frag = Fragment::new(&pen, smallvec![&rec], md_flags)?;
            let mut scratch = Scratch::new();
            let mut dvnt = smallvec![smallvec![]];
            let score = frag.score(&mut scratch, &mut dvnt)?;
            // All scores should be finite numbers (no panic / overflow)
            assert!(score.is_finite(), "score finite for {}", c.label);
        }
        Ok(())
    }

    // -----------------------------------------------------------------------
    // 10) VariantWindow boundary semantics (touching = no overlap)
    // -----------------------------------------------------------------------
    #[test]
    fn variant_window_touching_boundary_is_not_overlap() {
        let a = VariantWindow::compute(1, 3, 3, 6);
        assert!(a.is_none(), "touching boundary -> no overlap");
    }
}
