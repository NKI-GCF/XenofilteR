//! Benchmark ahash vs SipHash for FragmentTable insert+lookup.
//! Validates the ahash adoption decision empirically.

use criterion::{black_box, criterion_group, criterion_main, BenchmarkId, Criterion};
use ahash::RandomState;
use std::collections::HashMap;

fn make_keys(n: usize) -> Vec<Box<[u8]>> {
    // Realistic Illumina read names: ~50 bytes.
    (0..n)
        .map(|i| {
            format!("A00123:45:HXXXXYYYY:1:1101:{i}:{i}").into_bytes().into_boxed_slice()
        })
        .collect()
}

fn bench_hashmap(c: &mut Criterion) {
    let mut g = c.benchmark_group("FragmentTable_insert_lookup");
    for n in [1_000usize, 10_000, 100_000] {
        let keys = make_keys(n);

        g.bench_with_input(BenchmarkId::new("ahash", n), &n, |b, _| {
            b.iter(|| {
                let mut map: HashMap<Box<[u8]>, u64, RandomState> =
                    HashMap::with_hasher(RandomState::new());
                for (i, k) in keys.iter().enumerate() {
                    map.insert(k.clone(), i as u64);
                }
                let hits: u64 = keys.iter().map(|k| map.get(k).copied().unwrap_or(0)).sum();
                black_box(hits)
            })
        });

        g.bench_with_input(BenchmarkId::new("siphash", n), &n, |b, _| {
            b.iter(|| {
                let mut map: HashMap<Box<[u8]>, u64> = HashMap::new();
                for (i, k) in keys.iter().enumerate() {
                    map.insert(k.clone(), i as u64);
                }
                let hits: u64 = keys.iter().map(|k| map.get(k).copied().unwrap_or(0)).sum();
                black_box(hits)
            })
        });
    }
    g.finish();
}

criterion_group!(benches, bench_hashmap);
criterion_main!(benches);
