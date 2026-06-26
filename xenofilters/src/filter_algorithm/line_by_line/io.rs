use super::core::LineByLine;
use crate::alignment::SimpleRec;
use anyhow::Result;
use noodles::sam::alignment::record::data::field::Tag;
use noodles::sam::alignment::record_buf::data::field::Value;
use noodles::sam::alignment::record_buf::RecordBuf;

impl<R: SimpleRec> LineByLine<R> {
    /// Insert a single-byte aux tag into `rec`.
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

    /// Write `rec` through stream `i`.
    ///
    /// `best_state`:
    /// - `Some(true)`  → winning alignment (→ `--output`)
    /// - `Some(false)` → losing alignment  (→ `--discarded-output`)
    /// - `None`        → ambiguous         (→ `--ambiguous-output`)
    pub(super) fn write_record(
        &mut self,
        i: usize,
        rec: RecordBuf,
        best_state: Option<bool>,
    ) -> Result<()> {
        match (i, best_state) {
            (i, Some(false)) => self.routing_counters[i << 1]      += 1,
            (i, Some(true))  => self.routing_counters[1 + (i << 1)] += 1,
            (i, None)        => self.routing_counters[16 + i]        += 1,
        }
        if let Some(aln) = self.aln.get_mut(i) {
            aln.write_record(rec, best_state)
        } else {
            Ok(())
        }
    }

    /// Emit per-stream counters to the tracing backend (INFO level).
    ///
    /// Counter layout (index → meaning):
    /// ```text
    /// i*2+0  : discarded from alignment i
    /// i*2+1  : assigned to alignment i
    /// 16+i   : ambiguous for alignment i
    /// 24+i   : unmapped-discarded for alignment i
    /// ```
    pub(super) fn print_counters(&self, i: usize) {
        tracing::info!(
            stream = i,
            discarded  = self.routing_counters[i << 1],
            assigned  = self.routing_counters[1 + (i << 1)],
            ambiguous = self.routing_counters[16 + i],
            unmapped_discarded = self.routing_counters[24 + i],
            "Stream summary"
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
