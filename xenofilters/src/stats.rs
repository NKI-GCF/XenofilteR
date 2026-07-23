//! Machine-readable JSON summary. MultiQC custom-content compatible.
//!
//! MultiQC format reference:
//! <https://multiqc.info/docs/development/custom_content/#json-data>

use crate::Error;
use crate::filter_algorithm::COUNTER_STRIDE;
use serde::Serialize;
use std::collections::HashMap;
use std::path::Path;

#[derive(Serialize)]
struct MultiQcData<'a> {
    id: &'static str,
    plot_type: &'static str,
    pconfig: PConfig,
    data: HashMap<&'a str, StreamStats>,
}

#[derive(Serialize)]
struct PConfig {
    namespace: &'static str,
}

#[derive(Serialize)]
pub(crate) struct StreamStats {
    pub(crate) fragments_total: u64,
    pub(crate) assigned: u64,
    pub(crate) discarded: u64,
    pub(crate) ambiguous: u64,
    pub(crate) chimeric: u64,
    pub(crate) assigned_pct: f64,
    pub(crate) ambiguous_pct: f64,
    pub(crate) chimeric_pct: f64,
}

impl StreamStats {
    pub(crate) fn from_counters(counters: &[u64], nr: usize) -> Self {
        let b = nr * COUNTER_STRIDE;
        let disc = counters[b];
        let out = counters[b + 1];
        let ambig = counters[b + 2];
        let chimeric = counters[b + 3];
        let total = disc + out + ambig + chimeric;
        let pct = |n: u64| {
            if total == 0 {
                0.0
            } else {
                n as f64 * 100.0 / total as f64
            }
        };
        StreamStats {
            fragments_total: total,
            assigned: out,
            discarded: disc,
            ambiguous: ambig,
            chimeric,
            assigned_pct: pct(out),
            ambiguous_pct: pct(ambig),
            chimeric_pct: pct(chimeric),
        }
    }
}

pub(crate) fn write_stats(
    path: &Path,
    counters: &[u64],
    stream_count: usize,
    labels: &[String],
    sample_name: &str,
) -> Result<(), Error> {
    let _data: HashMap<String, StreamStats> = HashMap::new();
    for nr in 0..stream_count {
        let label = labels.get(nr).map(|s| s.as_str()).unwrap_or("unknown");
        let key = format!("{sample_name}:stream_{nr}:{label}");
        // We need owned keys; use a Vec to back them.
        let _ = (key, StreamStats::from_counters(counters, nr));
    }
    // Rebuild with owned strings to satisfy lifetime:
    let owned: Vec<(String, StreamStats)> = (0..stream_count)
        .map(|nr| {
            let label = labels.get(nr).map(|s| s.as_str()).unwrap_or("unknown");
            let key = format!("{sample_name}:stream_{nr}:{label}");
            (key, StreamStats::from_counters(counters, nr))
        })
        .collect();
    let _map: HashMap<&str, &StreamStats> = owned.iter().map(|(k, v)| (k.as_str(), v)).collect();

    let _doc = MultiQcData {
        id: "xenofilters",
        plot_type: "generalstats",
        pconfig: PConfig {
            namespace: "xenofilters",
        },
        data: owned
            .iter()
            .map(|(k, _v)| {
                (k.as_str(), StreamStats::from_counters(counters, 0)) // placeholder
            })
            .collect(),
    };
    // Re-derive correctly:
    let data_map: HashMap<String, StreamStats> = owned.into_iter().collect();
    let final_doc = serde_json::json!({
        "id": "xenofilters",
        "plot_type": "generalstats",
        "pconfig": { "namespace": "xenofilters" },
        "data": data_map,
    });
    let f = std::fs::File::create(path)?;
    serde_json::to_writer_pretty(f, &final_doc)?;
    Ok(())
}
