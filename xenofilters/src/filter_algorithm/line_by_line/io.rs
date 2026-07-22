use super::core::{LineByLine, COUNTER_STRIDE};
use crate::alignment::SimpleRec;
use crate::Error;
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
    ) -> Result<(), Error> {
        let tag = Tag::new(field[0], field[1]);
        let val = Value::from(value);
        rec.data_mut().insert(tag, val);
        Ok(())
    }

    /// Write `rec` through stream `i`. Counter layout per stream -- stride 4:
    ///
    /// `best_state`:
    /// - `Some(false)` -> nr*4+0 -> discard    (includes unmapped-discarded when --discard-unmapped)
    /// - `Some(true)`  -> nr*4+1 -> out/winner
    /// - `None`        -> nr*4+2 -> ambiguous  (includes unmapped-ambiguous when configured)
    ///   <elsewhere>     nr*4+3 -> chimeric   (XC:Z: tagged, both streams count)
    pub(super) fn write_record(
        &mut self,
        i: usize,
        rec: RecordBuf,
        best_state: Option<bool>,
    ) -> Result<(), Error> {
        let base = i * COUNTER_STRIDE;
        match best_state {
            Some(false) => self.routing_counters[base] += 1,
            Some(true) => self.routing_counters[base + 1] += 1,
            None => self.routing_counters[base + 2] += 1,
        }
        if let Some(aln) = self.aln.get_mut(i) {
            aln.write_record(rec, best_state)
        } else {
            Ok(())
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::config::run_config::RunConfig;
    use crate::tests::create_record;
    use smallvec::smallvec;

    #[test]
    fn test_add_aux_tags_inserts_expected_tag_and_value() -> Result<(), Error> {
        let config = RunConfig::default();
        let mut lbl: LineByLine<RecordBuf> = LineByLine::new(&config, smallvec![], vec![], vec![])?;
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
