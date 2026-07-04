//! Property-based tests targeting parser correctness.
//! These must never panic regardless of input.
//! Expected to be included from src/main.rs #[cfg(test)] tests module.

use proptest::prelude::*;
use crate::alignment::pre_assess::{md_mismatches, match_count_raw};

/// Strategy: arbitrary MD string containing only characters that can appear
/// in real MD tags. Does NOT include truly malformed bytes (tested separately).
fn arb_md() -> impl Strategy<Value = Vec<u8>> {
    prop::collection::vec(
        prop_oneof![
            (b'0'..=b'9').prop_map(|b| b),
            Just(b'A'), Just(b'C'), Just(b'G'), Just(b'T'), Just(b'N'),
            Just(b'^'),
        ],
        0..300,
    )
}

/// Strategy: arbitrary BAM-encoded CIGAR bytes (4 bytes per op, LE u32).
/// Opcodes 0–8 only (known valid). Length 1–300.
fn arb_cigar_bytes() -> impl Strategy<Value = Vec<u8>> {
    prop::collection::vec(
        (0u32..9u32, 1u32..=300u32)
            .prop_map(|(op, len)| ((len << 4) | op).to_le_bytes()),
        0..=20,
    )
    .prop_map(|chunks| chunks.into_iter().flatten().collect())
}

proptest! {
    /// md_mismatches must never panic and must return a value ≤ MD length.
    #[test]
    fn md_mismatches_no_panic_bounded(md in arb_md()) {
        let result = md_mismatches(&md);
        prop_assert!(result <= md.len());
    }

    /// md_mismatches on arbitrary bytes (including control chars) must not panic.
    #[test]
    fn md_mismatches_arbitrary_bytes_no_panic(md in prop::collection::vec(any::<u8>(), 0..300)) {
        let _ = md_mismatches(&md);
    }

    /// match_count_raw must never panic and must return a count ≤
    /// total bases implied by M/X/= ops in the CIGAR.
    #[test]
    fn match_count_raw_no_panic(
        cigar in arb_cigar_bytes(),
        md    in arb_md(),
    ) {
        let result = match_count_raw(&cigar, &md);
        // Upper bound: sum of M/X/= op lengths
        let cigar_m_total: usize = cigar.chunks_exact(4).map(|c| {
            let enc = u32::from_le_bytes([c[0], c[1], c[2], c[3]]);
            let op  = enc & 0xF;
            let len = (enc >> 4) as usize;
            if matches!(op, 0 | 7 | 8) { len } else { 0 }
        }).sum();
        prop_assert!(result <= cigar_m_total);
    }

    /// ScoreOpIter built from arbitrary CIGAR+MD must not panic.
    /// Errors (Err variant) are acceptable; panic is not.
    #[test]
    fn score_op_iter_no_panic(
        cigar_str in "[0-9]{1,3}[MIDNSHP=X]([0-9]{1,3}[MIDNSHP=X]){0,9}",
        md        in "[0-9A-Z\\^]{0,50}",
        is_rev    in any::<bool>(),
    ) {
        use crate::{alignment::{MdCigFlags, ScoreOpIter}, tests::create_record};
        let qual = vec![30u8; 200];
        // create_record may itself fail (malformed CIGAR); that is acceptable.
        let Ok(rec) = create_record(b"r", &cigar_str, &[], &qual, &md, is_rev) else {
            return Ok(());
        };
        let flags = rec.flags();
        let Ok(mcf) = MdCigFlags::try_from_record(&rec, &flags) else {
            return Ok(());
        };
        // Collecting all ops must not panic; errors are fine.
        let _ops: Vec<_> = ScoreOpIter::new(&mcf).collect();
    }

    /// is_perfect must never panic.
    #[test]
    fn is_perfect_no_panic(
        cigar_str in "[0-9]{1,3}[MIDNSHP=X]([0-9]{1,3}[MIDNSHP=X]){0,5}",
        md        in "[0-9A-Z^]{0,30}",
    ) {
        use crate::{alignment::MdCigFlags, tests::create_record};
        let qual = vec![30u8; 200];
        let Ok(rec) = create_record(b"r", &cigar_str, &[], &qual, &md, false) else {
            return Ok(());
        };
        let flags = rec.flags();
        let Ok(mcf) = MdCigFlags::try_from_record(&rec, &flags) else {
            return Ok(());
        };
        let _ = mcf.is_perfect();
    }

    /// wis_max_rescue_delta on arbitrary delta values must not panic
    /// and must return a non-negative value.
    #[test]
    fn wis_returns_nonnegative(deltas in prop::collection::vec(-1e6f64..1e6f64, 0..50)) {
        use crate::{
            alignment::fragment::wis_max_rescue_delta,
            variant::{Eval, Variant},
        };
        struct V { pos: usize, delta: f64 }
        impl Variant for V {
            fn pos(&self) -> usize { self.pos }
            fn ref_allele(&self) -> &[u8] { b"A" }
            fn alt_allele(&self) -> &[u8] { b"G" }
            fn p_variant(&self) -> f64 { 0.9 }
        }
        let vs: Vec<&'static V> = deltas
            .iter()
            .enumerate()
            .map(|(i, _)| Box::leak(Box::new(V { pos: i * 5, delta: 0.0 })) as &'static V)
            .collect();
        let mut evals: smallvec::SmallVec<[Eval<'_>; 4]> = vs
            .iter()
            .zip(deltas.iter())
            .map(|(v, &d)| {
                let mut e = Eval::new();
                e.set_variant(*v as &dyn Variant);
                e.update(0.0, d);
                e
            })
            .collect();
        let mut dvnt = smallvec::smallvec![evals];
        let mut dp   = smallvec::smallvec![];
        let result = wis_max_rescue_delta(&mut dvnt, &mut dp);
        prop_assert!(result >= 0.0, "WIS must return non-negative: got {result}");
        prop_assert!(!result.is_nan(), "WIS must not return NaN");
    }
}
