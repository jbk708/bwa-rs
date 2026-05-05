//! Benchmarks for wavelet tree occ queries vs sampling.
//!
//! Run with: cargo bench

use bwa_mem::{FMIndex, Reference};
use criterion::{black_box, criterion_group, criterion_main, BenchmarkId, Criterion};

fn create_test_index(size: usize) -> FMIndex {
    let sequence: Vec<u8> = (0..size)
        .map(|i| match i % 4 {
            0 => b'A',
            1 => b'C',
            2 => b'G',
            _ => b'T',
        })
        .collect();
    let fasta = format!(">test\n{}", String::from_utf8(sequence).unwrap());
    let reference = Reference::parse_fasta(&fasta).unwrap();
    FMIndex::build(&reference)
}

fn occ_query_benchmark(c: &mut Criterion) {
    let sizes = [1_000, 10_000, 100_000, 1_000_000];

    let mut group = c.benchmark_group("wavelet_occ_query");

    for size in sizes {
        let index = create_test_index(size);
        let pattern = [0u8, 1, 2, 3]; // ACGT pattern

        group.bench_with_input(BenchmarkId::from_parameter(size), &size, |b, _| {
            b.iter(|| {
                // Perform multiple occ queries like FM-index search does
                let (mut left, mut right) = (0usize, index.len);
                for &c in pattern.iter().rev() {
                    let occ_left = index.occ.occ(c, left);
                    let occ_right = index.occ.occ(c, right);
                    let f_c = (index.f_column[c as usize] as usize).saturating_sub(1);
                    left = f_c + occ_left as usize;
                    right = f_c + occ_right as usize;
                    if left >= right {
                        break;
                    }
                }
                black_box((left, right))
            });
        });
    }

    group.finish();
}

fn naive_occ_query(seq: &[u8], c: u8, idx: usize) -> usize {
    seq[..idx.min(seq.len())].iter().filter(|&&x| x == c).count()
}

fn occ_query_comparison(c: &mut Criterion) {
    let sizes = [1_000, 10_000, 100_000];

    let mut group = c.benchmark_group("occ_comparison");

    for size in sizes {
        let index = create_test_index(size);
        let bwt = index.bwt.as_slice().to_vec();
        let chars: [u8; 4] = [0, 1, 2, 3];

        group.bench_with_input(BenchmarkId::from_parameter(size), &size, |b, _| {
            b.iter(|| {
                for &c in &chars {
                    for idx in (0..size).step_by(size / 10).take(10) {
                        let wavelet_result = index.occ.occ(c, idx);
                        let naive_result = naive_occ_query(&bwt, c, idx);
                        assert_eq!(wavelet_result, naive_result as u32);
                    }
                }
            });
        });
    }

    group.finish();
}

criterion_group! {
    name = benches;
    config = Criterion::default();
    targets = occ_query_benchmark, occ_query_comparison
}
criterion_main!(benches);