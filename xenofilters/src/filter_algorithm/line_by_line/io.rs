use super::core::{COUNTER_STRIDE, LineByLine};
use crate::Error;
use crate::alignment::SimpleRec;
use noodles::sam::alignment::record_buf::RecordBuf;

impl<R: SimpleRec> LineByLine<R> {
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
