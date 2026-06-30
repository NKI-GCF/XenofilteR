use super::chimeric::ChimericDecision;
use super::core::FragmentBuffer;
use super::core::LineByLine;
use crate::Error;
use crate::alignment::SimpleRec;
use noodles::sam::alignment::record::data::field::Tag;
use noodles::sam::alignment::record_buf::RecordBuf;
use noodles::sam::alignment::record_buf::data::field::Value;

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
    ) -> Result<(), Error> {
        match (i, best_state) {
            (i, Some(false)) => self.routing_counters[i * 4] += 1,
            (i, Some(true)) => self.routing_counters[1 + (i * 4)] += 1,
            (i, None) => self.routing_counters[2 + (i * 4)] += 1,
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
            discarded = self.routing_counters[i * 4],
            assigned = self.routing_counters[1 + (i * 4)],
            ambiguous = self.routing_counters[2 + (i * 4)],
            chimeric = self.routing_counters[3 + (i * 4)],
            "Stream summary"
        );
    }
    /// Write all records from a chimeric fragment.
    ///
    /// For streams in the chimeric pair: records go to assigned output with
    /// `XC:Z:<other_stream_label>` tag.  For every other stream present in
    /// `best`: records go to filtered (discarded) output.
    ///
    /// The `XC:Z:` tag value is the human-readable label of the *other* stream
    /// (configured via `--stream-labels`), making it easy to filter chimeric
    /// reads with `samtools view -d XC:hpv`.
    pub(super) fn emit_chimeric(
        &mut self,
        best: &mut FragmentBuffer<R>,
        decision: ChimericDecision,
    ) -> Result<(), Error> {
        let (chimeric_a, chimeric_b, kind) = match decision {
            ChimericDecision::Chimeric {
                stream_a,
                stream_b,
                kind,
            } => (stream_a, stream_b, kind),
            ChimericDecision::Normal => {
                unreachable!("emit_chimeric called on non-chimeric decision")
            }
        };

        let label_a = self.chimeric_label(chimeric_a);
        let label_b = self.chimeric_label(chimeric_b);
        let tag_xc = Tag::new(b'X', b'C');

        tracing::debug!(
            kind = ?kind,
            stream_a = chimeric_a,
            stream_b = chimeric_b,
            label_a  = %label_a,
            label_b  = %label_b,
            "Emitting chimeric fragment"
        );

        best.drain(..)
            .try_for_each(|mut state| -> Result<(), Error> {
                let nr = state.get_nr();
                let is_chimeric_stream = nr == chimeric_a || nr == chimeric_b;

                if is_chimeric_stream {
                    let other_label = if nr == chimeric_a { &label_b } else { &label_a };
                    let xc_value = Value::String(other_label.as_bytes().into());

                    // For ReadSplit chimerism, stream A may contain a supplementary
                    // alignment that is a false-positive mapping of the split read's
                    // complementary portion to the wrong reference.  It is written
                    // here with the XC:Z tag so that its provenance is clear;
                    // downstream tools can identify it via the supplementary flag
                    // (SAM 0x800 / `samtools view -F 2048`).
                    //
                    // A future `--chimeric-suppress-supplementary` option could
                    // silently drop these records.
                    state
                        .drain_records()
                        .try_for_each(|r| -> Result<(), Error> {
                            let header = self.aln[nr].header();
                            let mut rb = r.as_record_buf(header)?;
                            rb.data_mut().insert(tag_xc, xc_value.clone());
                            self.routing_counters[3 + (nr * 4)] += 1;
                            self.aln[nr].write_record(rb, Some(true))
                        })
                } else {
                    // Streams outside the chimeric pair are discarded normally.
                    state
                        .drain_records()
                        .try_for_each(|r| -> Result<(), Error> {
                            let header = self.aln[nr].header();
                            let rb = r.as_record_buf(header)?;
                            self.routing_counters[nr << 1] += 1;
                            self.aln[nr].write_record(rb, Some(false))
                        })
                }
            })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::config::Config;
    use crate::tests::create_record;
    use smallvec::smallvec;

    #[test]
    fn test_add_aux_tags_inserts_expected_tag_and_value() -> Result<(), Error> {
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
