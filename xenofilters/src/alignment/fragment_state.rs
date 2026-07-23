use crate::Error;
use crate::alignment::mate_kind::{MateClassifiable, MateKind, mate_slot, segment_id};
use crate::alignment::{Fragment, MdCigFlags, SimpleRec, tie_break_bool};
use crate::filter_algorithm::line_by_line::{READ_CT, Scratch};
use crate::penalty::Penalty;
use crate::region::{PositiveRegions, ScoreFn};
use crate::variant::{FragEvalVec, StoreTrait};
use noodles::sam::alignment::record::Cigar;
use noodles::sam::alignment::record::Flags;
use smallvec::{SmallVec, smallvec};
use std::cmp::Ordering;

#[derive(PartialEq, Debug)]
pub(crate) struct FragmentState<R> {
    pub(crate) flags: SmallVec<[Flags; 8]>,
    records: SmallVec<[R; 8]>,
    species_nr: usize,
    bisulfite: bool,
}

type McfPair<'f> = (SmallVec<[MdCigFlags<'f>; 8]>, SmallVec<[MdCigFlags<'f>; 8]>);

impl<R: SimpleRec> FragmentState<R> {
    pub(crate) fn from_record(r: R, species_nr: usize, bisulfite: bool) -> Result<Self, Error> {
        Ok(FragmentState {
            flags: smallvec![r.flags()?],
            records: smallvec![r],
            species_nr,
            bisulfite,
        })
    }
    pub(crate) fn add_record(&mut self, r: R) -> Result<(), Error> {
        self.flags.push(r.flags()?);
        self.records.push(r);
        Ok(())
    }
    pub(crate) fn get_records(&self) -> &[R] {
        &self.records
    }
    pub(crate) fn drain_records(&mut self) -> impl Iterator<Item = R> + '_ {
        self.records.drain(..)
    }
    #[must_use]
    pub(crate) fn first_qname(&self) -> &[u8] {
        self.records
            .first()
            .and_then(|r| r.name())
            .map_or(b"", |n| n.as_ref())
    }
    #[must_use]
    pub(crate) fn get_nr(&self) -> usize {
        self.species_nr
    }
    pub(crate) fn flags(&self, i: usize) -> Option<&Flags> {
        self.flags.get(i)
    }
    /// True when all primaries in this fragment are unmapped.
    /// Used for O(1) per-stream Tier 1 classification in the N-way tournament.
    pub(crate) fn is_all_unmapped(&self) -> bool {
        let f = &self.flags[0];
        f.is_unmapped() && (!f.is_segmented() || f.is_mate_unmapped())
    }

    /// Build `MdCigFlags` for every record in the fragment, in record order.
    ///
    /// Borrows from `self` with lifetime `'f`; the returned SmallVec keeps those
    /// borrows alive.  Called before any mutable access to the fragment buffer so
    /// the borrow checker sees no aliasing.  Matches the record-order iteration
    /// expected by `score_candidate`.
    pub(crate) fn build_mcfs<'f>(&'f self) -> Result<SmallVec<[MdCigFlags<'f>; READ_CT]>, Error> {
        self.flags
            .iter()
            .zip(self.records.iter())
            .map(|(flags, rec)| MdCigFlags::try_from_record(rec, flags, self.bisulfite))
            .collect()
    }

    pub(crate) fn order_mates(&self) -> SmallVec<[usize; READ_CT]> {
        let len = self.records.len();
        let mut indices: SmallVec<[(u8, usize, usize, usize); 2]> = SmallVec::with_capacity(len);
        for i in 0..len {
            let r = &self.records[i];
            let start = match r.alignment_start() {
                Some(Ok(pos)) => pos.get(),
                _ => panic!("Mapped record has no alignment start"),
            };
            let tid = match r.ref_seq_id() {
                Some(Ok(tid)) => tid,
                _ => panic!("Mapped record has no reference sequence ID"),
            };
            let flags = &self.flags[i];
            let pos = if flags.is_reverse_complemented() {
                start + r.cigar().len()
            } else {
                start
            };
            let ord = u8::from(flags.is_last_segment()) * 2 + u8::from(flags.is_secondary());
            indices.push((ord, tid, pos, i));
        }
        indices.sort();
        indices.iter().map(|t| t.3).collect()
    }
    pub(crate) fn cmp_perfect<'f>(
        &'f self,
        other: &'f FragmentState<R>,
        ord: &mut Option<Ordering>,
    ) -> Result<McfPair<'f>, Error> {
        let mut mcfs1: SmallVec<[MdCigFlags<'f>; 8]> =
            SmallVec::with_capacity(self.get_records().len());
        let mut perfect_self = true;
        for (flags, record) in self.flags.iter().zip(self.records.iter()) {
            let mcf = MdCigFlags::try_from_record(record, flags, self.bisulfite)?;
            if flags.is_supplementary() || !mcf.is_perfect() {
                perfect_self = false;
            }
            mcfs1.push(mcf);
        }

        let mut mcfs2: SmallVec<[MdCigFlags<'f>; 8]> =
            SmallVec::with_capacity(other.get_records().len());
        let mut perfect_other = true;
        for (flags, record) in other.flags.iter().zip(other.records.iter()) {
            let mcf = MdCigFlags::try_from_record(record, flags, self.bisulfite)?;
            if flags.is_supplementary() || !mcf.is_perfect() {
                perfect_other = false;
            }
            mcfs2.push(mcf);
        }

        *ord = tie_break_bool(perfect_self, perfect_other);
        Ok((mcfs1, mcfs2))
    }
    /// Single-stream NW scoring loop, shared by LineByLine, Collated, and HashLookup.
    ///
    /// Parameters match the calling convention previously duplicated in
    /// `score_candidate`, `score_candidate_owned`, and `nw_score_fragment`.
    ///
    /// Preconditions:
    /// - `mcfs` ordered identically to `self.order_mates()` output for
    ///   primary/supplementary records; unmapped records must be absent from
    ///   `mcfs` (MdCigFlags::try_from_record rejects them).
    /// - `cancel_slot[i]` true iff mate slot i's NW contribution is provably
    ///   equal across all competing streams (same read, same quality string).
    pub(crate) fn score_state_nw(
        &self,
        mcfs: SmallVec<[MdCigFlags<'_>; READ_CT]>,
        store: Option<&dyn StoreTrait>,
        penalties: &Penalty,
        scratch: &mut Scratch,
        cancel_slot: [bool; 2],
        positive_regions: Option<(&PositiveRegions, ScoreFn)>,
    ) -> Result<f64, Error> {
        let mut segment: SmallVec<[&R; READ_CT]> = SmallVec::new();
        let mut seg_mcfs: SmallVec<[MdCigFlags<'_>; READ_CT]> = SmallVec::new();
        let mut dvnt: FragEvalVec<'_> = SmallVec::new();
        let mut supp_penalty = 0.0f64;
        let has_variants = store.is_some();

        // Positive-region bonus: sum over all primary mapped records.
        let region_bonus = if let Some((regions, score_fn)) = positive_regions {
            self.get_records()
                .iter()
                .zip(self.flags.iter())
                .filter(|(_, f)| !f.is_unmapped() && !f.is_secondary() && !f.is_supplementary())
                .map(|(rec, f)| {
                    let is_rev = f.is_reverse_complemented();
                    let ref_id = rec
                        .ref_seq_id()
                        .transpose()
                        .ok()
                        .flatten()
                        .unwrap_or(usize::MAX);
                    let r_start = rec
                        .alignment_start()
                        .transpose()
                        .ok()
                        .flatten()
                        .map(|p| p.get().saturating_sub(1))
                        .unwrap_or(0);
                    let r_end = r_start + mcfs.first().map(|m| m.get_cigar().len()).unwrap_or(0);
                    regions
                        .overlapping(ref_id, r_start, r_end, is_rev)
                        .map(|region| {
                            let ob = region.overlap_bases(r_start, r_end);
                            score_fn.apply(region.score, ob, region.len())
                        })
                        .sum::<f64>()
                })
                .sum::<f64>()
        } else {
            0.0
        };

        let mut mcfs_opt: SmallVec<[Option<MdCigFlags<'_>>; READ_CT]> =
            mcfs.into_iter().map(Some).collect();

        for idx in self.order_mates() {
            let flags = self.flags(idx).ok_or(Error::NoFlagsForRecord { idx })?;
            if flags.is_secondary() {
                mcfs_opt[idx] = None;
                continue;
            }
            let slot = mate_slot(segment_id(flags));
            if cancel_slot[slot] {
                mcfs_opt[idx] = None;
                continue;
            }
            if flags.is_supplementary() {
                supp_penalty += penalties.chimeric_junction_penalty;
            }
            let rec = &self.get_records()[idx];
            if flags.is_unmapped() {
                dvnt.push(SmallVec::new());
                // Unmapped records have no MdCigFlags; skip NW but preserve dvnt alignment.
                mcfs_opt[idx] = None;
                continue;
            }
            if has_variants {
                let tid = rec
                    .ref_seq_id()
                    .transpose()?
                    .ok_or(Error::MappedRecordNoReferenceSequenceId)?;
                let start = rec
                    .alignment_start()
                    .transpose()?
                    .ok_or(Error::MappedRecordNoAlignmentStart)?
                    .get();
                let cig_len = mcfs_opt[idx]
                    .as_ref()
                    .ok_or(Error::MdCigFlagsConsumed { idx })?
                    .get_cigar()
                    .len();
                dvnt.push(
                    store
                        .unwrap()
                        .overlapping_multi(tid, start, start + cig_len),
                );
            } else {
                dvnt.push(SmallVec::new());
            }
            segment.push(rec);
            seg_mcfs.push(
                mcfs_opt[idx]
                    .take()
                    .ok_or(Error::MdCigFlagsAlreadyConsumed { idx })?,
            );
        }

        if segment.is_empty() {
            return Ok(supp_penalty);
        }
        let base = Fragment::new(penalties, segment, seg_mcfs)?
            .score(scratch, &mut dvnt)
            .map_err(|e| Error::NeedlemanWunschError(e.to_string()))?;

        Ok(base + supp_penalty + region_bonus)
    }
}

impl<R: SimpleRec> PartialOrd for FragmentState<R> {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        tie_break_bool(other.is_all_unmapped(), self.is_all_unmapped())
    }
}

impl<R: SimpleRec> MateClassifiable for FragmentState<R> {
    /// Classifies each present mate without erroring on unmapped records:
    /// unmapped is detected from flags directly, never passed into
    /// `MdCigFlags::try_from_record` (which requires a mapped record).
    fn mate_kinds(&self) -> [Option<MateKind>; 2] {
        let mut kinds: [Option<MateKind>; 2] = [None, None];
        for (flags, rec) in self.flags.iter().zip(self.records.iter()) {
            if flags.is_secondary() || flags.is_supplementary() {
                continue; // only primaries determine the mate's classification
            }
            let slot = mate_slot(segment_id(flags));
            let kind = if flags.is_unmapped() {
                MateKind::Unmapped
            } else {
                match MdCigFlags::try_from_record(rec, flags, self.bisulfite) {
                    Ok(mcf) if mcf.is_perfect() => MateKind::Perfect,
                    Ok(_) => MateKind::Other,
                    Err(_) => MateKind::Other, // malformed MD/CIGAR -- must be scored
                }
            };
            kinds[slot] = Some(kind);
        }
        kinds
    }
}

#[cfg(test)]
mod tests;
