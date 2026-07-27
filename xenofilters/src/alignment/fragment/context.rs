//! Scoring window and variant-within-window coordinate contexts.
//! These are purely data-carrying structs; all scoring logic lives in `mod.rs`.

/// Identifies a scoring window within a segment.
pub(super) struct WindowCtx {
    pub(super) dvnt_i: usize,
    pub(super) seg_i: usize,
    pub(super) ref_start: usize,
    pub(super) ref_end: usize,
}

impl WindowCtx {
    pub(super) fn new(dvnt_i: usize, seg_i: usize, ref_start: usize, ref_end: usize) -> Self {
        Self {
            dvnt_i,
            seg_i,
            ref_start,
            ref_end,
        }
    }

    pub(super) fn to_variant_ctx(&self, dvnt_j: usize) -> VariantCtx {
        VariantCtx {
            dvnt_i: self.dvnt_i,
            dvnt_j,
            seg_i: self.seg_i,
            ref_start: self.ref_start,
            ref_end: self.ref_end,
        }
    }
}

/// Identifies a specific variant within a scoring window.
pub(super) struct VariantCtx {
    pub(super) dvnt_i: usize,
    pub(super) dvnt_j: usize,
    pub(super) seg_i: usize,
    pub(super) ref_start: usize,
    pub(super) ref_end: usize,
}
