use super::core::LineByLine;
use crate::alignment::stringify_record;
use anyhow::{Result, anyhow};
use rust_htslib::bam::record::{Aux, Record};

impl LineByLine {
    pub(super) fn add_aux_tags(
        &mut self,
        rec: &mut Record,
        field: &[u8],
        value: Aux,
    ) -> Result<()> {
        rec.push_aux(field, value).map_err(|e| {
            anyhow!(
                "Error adding {} tag to record: {}\n{}",
                field.iter().map(|&b| b as char).collect::<String>(),
                e,
                stringify_record(rec)
            )
        })
    }
    pub(super) fn write_record(
        &mut self,
        i: usize,
        rec: Record,
        best_state: Option<bool>,
    ) -> Result<()> {
        match (i, best_state) {
            (i, Some(false)) => self.branch_counters[i << 1] += 1,
            (i, Some(true)) => self.branch_counters[1 + (i << 1)] += 1,
            (i, None) => {
                if (self.is_unmapped_skipped)(&rec) {
                    self.branch_counters[24 + i] += 1;
                    return Ok(());
                }
                self.branch_counters[16 + i] += 1;
            }
        }
        if let Some(aln) = self.aln.get_mut(i) {
            aln.write_record(rec, best_state)
        } else {
            Ok(())
        }
    }

    pub(super) fn print_counters(&self, i: usize) {
        eprintln!(
            "[{}]: Filtered from alignment {i}: {}",
            i << 1,
            self.branch_counters[i << 1]
        );
        eprintln!(
            "[{}]: Assigned to alignment {i}: {}",
            1 + (i << 1),
            self.branch_counters[1 + (i << 1)]
        );
        eprintln!(
            "[{}]: Ambiguous for alignment {i}: {}",
            16 + i,
            self.branch_counters[16 + i]
        );
        eprintln!(
            "[{}]: Unmapped filtered for alignment {i}: {}",
            24 + i,
            self.branch_counters[24 + i]
        );
    }
}
