// src/filter_algorithm/hash_lookup/name_codec.rs
use crate::aln_stream::AlignmentStream;
use crate::config::{run_config::RunConfig, NameEncoderKind};
use crate::{ensure, Error};
use noodles::sam::alignment::record_buf::RecordBuf;
use smallvec::SmallVec;

pub(crate) trait NameEncoder: Send + Sync {
    /// Encode a stripped (suffix-removed) read name to a heap key.
    fn encode(&self, name: &[u8]) -> Box<[u8]>;
    /// Called on stream 0's first name. Returns prefix bytes for stream 1 verification.
    fn detect(&mut self, _name: &[u8]) -> &[u8] {
        b""
    }
    /// Called on stream 1's first name. Implementations may disable compression on mismatch.
    fn verify(&mut self, _name: &[u8], _expected: &[u8]) -> bool {
        true
    }
}

/// Passthrough: full name -> heap key. Correct for any platform.
pub(crate) struct PassthroughEncoder;
impl NameEncoder for PassthroughEncoder {
    fn encode(&self, name: &[u8]) -> Box<[u8]> {
        name.into()
    }
}

/// Illumina: strip machine:run:flowcell prefix (3 colon-delimited fields).
/// Falls back to PassthroughEncoder on format mismatch or cross-stream prefix divergence.
pub(crate) struct IlluminaEncoder {
    prefix_len: usize,
    detected: bool,
    prefix_buf: Vec<u8>,
}

impl IlluminaEncoder {
    pub(crate) fn new() -> Self {
        Self {
            prefix_len: 0,
            detected: false,
            prefix_buf: Vec::new(),
        }
    }
    const STRIP_FIELDS: usize = 3;
}

impl NameEncoder for IlluminaEncoder {
    fn detect(&mut self, name: &[u8]) -> &[u8] {
        let mut colons = 0;
        for (i, &b) in name.iter().enumerate() {
            if b == b':' {
                colons += 1;
                if colons == Self::STRIP_FIELDS {
                    self.prefix_len = i + 1;
                    self.prefix_buf = name[..self.prefix_len].to_vec();
                    self.detected = true;
                    return &self.prefix_buf;
                }
            }
        }
        self.prefix_len = 0;
        b""
    }

    fn verify(&mut self, name: &[u8], expected: &[u8]) -> bool {
        if self.prefix_len > 0 && !name.starts_with(expected) {
            tracing::warn!("Read name prefix mismatch; disabling key compression");
            self.prefix_len = 0;
            false
        } else {
            true
        }
    }

    fn encode(&self, name: &[u8]) -> Box<[u8]> {
        if self.prefix_len > 0 && name.len() > self.prefix_len {
            name[self.prefix_len..].into()
        } else {
            name.into()
        }
    }
}

/// Select encoder from config. Unknown platform -> PassthroughEncoder.
pub(crate) fn make_encoder(
    cfg: &RunConfig,
    aln: &SmallVec<[Box<dyn AlignmentStream<RecordBuf>>; 2]>,
) -> Result<Box<dyn NameEncoder>, Error> {
    if let Some(encoder) = cfg.name_encoder {
        Ok(match encoder {
            NameEncoderKind::Illumina => Box::new(IlluminaEncoder::new()),
            NameEncoderKind::Passthrough => Box::new(PassthroughEncoder),
        })
    } else {
        // auto-detect based on first record of stream 0
        let name = aln[0].next_qname();
        let mut encoder = IlluminaEncoder::new();
        let first_prefix = encoder.detect(name);
        let n = aln.len();
        if !first_prefix.is_empty() {
            tracing::info!(
                "Detected Illumina read name format; prefix length {}",
                first_prefix.len()
            );
            for i in 1..n {
                let name_i = aln[i].next_qname();
                ensure!(
                    name_i.starts_with(first_prefix),
                    Error::InvalidInput(format!(
                        "Read name prefix mismatch between streams 0 and {}",
                        i
                    ))
                );
            }
            Ok(Box::new(encoder))
        } else {
            tracing::info!("Using PassthroughEncoder for read names");
            for i in 1..n {
                let name_i = aln[i].next_qname();
                let prefix = encoder.detect(name_i);
                ensure!(
                    prefix.is_empty(),
                    Error::InvalidInput(format!(
                        "Read name prefix mismatch between streams 0 and {}",
                        i
                    ))
                );
            }
            Ok(Box::new(PassthroughEncoder))
        }
    }
}
