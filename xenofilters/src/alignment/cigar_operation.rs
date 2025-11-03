use rust_htslib::bam::record::Cigar;

/// Our new, owned CIGAR enum
#[derive(Debug, Clone, PartialEq)]
pub enum CigarOp {
    Match(u32),
    Ins(u32),
    Del(u32),
    RefSkip(u32),
    SoftClip(u32),
    HardClip(u32),
    Pad(u32),
    Equal(u32),
    Diff(u32),
    /// Our new op for supplementary jumps.
    Translocate {
        new_tid: u32,
        new_pos: i64,
        new_strand: bool,
        // We can add a score/penalty here if we want
        gap_score: f64,
    },
}

impl From<&Cigar> for CigarOp {
    fn from(cigar: &Cigar) -> Self {
        match *cigar {
            Cigar::Match(l) => CigarOp::Match(l),
            Cigar::Ins(l) => CigarOp::Ins(l),
            Cigar::Del(l) => CigarOp::Del(l),
            Cigar::RefSkip(l) => CigarOp::RefSkip(l),
            Cigar::SoftClip(l) => CigarOp::SoftClip(l),
            Cigar::HardClip(l) => CigarOp::HardClip(l),
            Cigar::Pad(l) => CigarOp::Pad(l),
            Cigar::Equal(l) => CigarOp::Equal(l),
            Cigar::Diff(l) => CigarOp::Diff(l),
        }
    }
}

impl CigarOp {
    pub fn len(&self) -> u32 {
        match self {
            CigarOp::Match(l)
            | CigarOp::Ins(l)
            | CigarOp::Del(l)
            | CigarOp::RefSkip(l)
            | CigarOp::SoftClip(l)
            | CigarOp::HardClip(l)
            | CigarOp::Pad(l)
            | CigarOp::Equal(l)
            | CigarOp::Diff(l) => *l,
            CigarOp::Translocate { .. } => 0,
        }
    }
    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }
}

pub struct CigarOpIterator<'a> {
    cigars: &'a [CigarOp],
    index: usize,
}

impl<'a> CigarOpIterator<'a> {
    pub fn new(cigars: &'a [CigarOp]) -> Self {
        CigarOpIterator {
            cigars,
            index: 0,
        }
    }
}

impl<'a> Iterator for CigarOpIterator<'a> {
    type Item = &'a CigarOp;

    fn next(&mut self) -> Option<Self::Item> {
        if self.index < self.cigars.len() {
            let item = &self.cigars[self.index];
            self.index += 1;
            Some(item)
        } else {
            None
        }
    }
}
