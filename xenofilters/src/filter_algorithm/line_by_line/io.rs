use super::core::LineByLine;
use anyhow::Result;
use noodles::sam::alignment::record::data::field::Tag;
use noodles::sam::alignment::record_buf::data::field::Value;
use noodles::sam::alignment::record_buf::RecordBuf;
use crate::alignment::SimpleRec;

impl<R: SimpleRec> LineByLine<R> {
    pub(super) fn add_aux_tags(
        &mut self,
        rec: &mut RecordBuf,
        field: &[u8; 2],
        value: u8,
    ) -> Result<()> {
        let tag = Tag::new(field[0], field[1]);
        let val = Value::from(value);
        rec.data_mut().insert(tag, val);
        Ok(())
    }

    pub(super) fn write_record(
        &mut self,
        i: usize,
        rec: RecordBuf,
        best_state: Option<bool>,
    ) -> Result<()> {
        match (i, best_state) {
            (i, Some(false)) => self.branch_counters[i << 1] += 1,
            (i, Some(true)) => self.branch_counters[1 + (i << 1)] += 1,
            (i, None) => self.branch_counters[16 + i] += 1,
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
