use crate::variant::Variant;

pub(crate) struct VariantEval<'a> {
    vnt: Option<&'a dyn Variant>,
    // FIXME: a variant could have multiple alleles, so these should become smallvecs:
    alt_score: f64,
    incurred: f64,
}

impl<'a> VariantEval<'a> {
    pub(crate) fn new() -> Self {
        VariantEval {
            incurred: 0.0,
            alt_score: 0.0,
            vnt: None
        }
    }
    pub(crate) fn delta(&self) -> f64 {
        self.alt_score - self.incurred
    }
    pub(crate) fn end(&self) -> i64 {
        self.ref_end().max(self.alt_end())
    }
    pub(crate) fn vnt(&self) -> &'a dyn Variant {
        self.vnt.expect("VariantEval should always have a variant reference")
    }
    pub(crate) fn set_variant(&mut self, vnt: &'a dyn Variant) {
        self.vnt = Some(vnt);
    }
    pub(crate) fn start(&self) -> i64 {
        self.vnt().pos()
    }
    pub(crate) fn ref_end(&self) -> i64 {
        self.vnt().pos() + self.vnt().ref_allele().len() as i64
    }
    pub(crate) fn alt_end(&self) -> i64 {
        self.vnt().pos() + self.vnt().alt_allele().len() as i64
    }
    pub(crate) fn update(&mut self, add: f64, alt_score: f64) {
        self.incurred += add;
        self.alt_score += alt_score;
    }
}
