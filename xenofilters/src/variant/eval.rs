use crate::variant::Variant;

#[derive(Clone, Copy, Default)]
pub struct Eval<'a> {
    vnt: Option<&'a dyn Variant>,
    // FIXME: a variant could have multiple alleles, so these should become smallvecs:
    alt_score: f64,
    pub(crate) incurred: f64,
}

impl<'a> Eval<'a> {
    pub fn new() -> Self {
        Eval {
            incurred: 0.0,
            alt_score: 0.0,
            vnt: None,
        }
    }
    pub(crate) fn delta(&self) -> f64 {
        self.alt_score - self.incurred
    }
    pub(crate) fn end(&self) -> usize {
        self.ref_end().max(self.alt_end())
    }
    pub(crate) fn vnt(&self) -> &'a dyn Variant {
        self.vnt
            .expect("VariantEval should always have a variant reference")
    }
    pub fn set_variant(&mut self, vnt: &'a dyn Variant) {
        self.vnt = Some(vnt);
    }
    pub(crate) fn start(&self) -> usize {
        self.vnt().pos()
    }
    pub(crate) fn ref_end(&self) -> usize {
        self.vnt().pos() + self.vnt().ref_allele().len()
    }
    pub(crate) fn alt_end(&self) -> usize {
        self.vnt().pos() + self.vnt().alt_allele().len()
    }
    pub fn update(&mut self, add: f64, alt_score: f64) {
        self.incurred += add;
        self.alt_score += alt_score;
    }
}

#[cfg(test)]
mod tests;
