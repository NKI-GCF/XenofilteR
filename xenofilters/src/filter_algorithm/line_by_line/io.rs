use super::core::LineByLine;
use crate::alignment::SimpleRec;
use anyhow::Result;
use noodles::sam::alignment::record::data::field::Tag;
use noodles::sam::alignment::record_buf::data::field::Value;
use noodles::sam::alignment::record_buf::RecordBuf;

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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::config::Config;
    use crate::tests::create_record;
    use smallvec::smallvec;

    #[test]
    fn test_add_aux_tags_inserts_expected_tag_and_value() -> Result<()> {
        let mut lbl: LineByLine<RecordBuf> = LineByLine::new(Config::default(), smallvec![])?;
        let mut rec = create_record(b"r", "5M", &[], &[], "5", false)?;
        lbl.add_aux_tags(&mut rec, b"XF", 42)?;

        let tag = Tag::new(b'X', b'F');
        match rec.data().get(&tag) {
            Some(Value::UInt8(v)) => assert_eq!(*v, 42),
            other => panic!("unexpected tag value: {other:?}"),
        }
        Ok(())
    }
}
