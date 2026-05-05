//! Benchmarks for alignment parameter tuning.
//!
//! Run with: cargo bench

use bwa_mem::{alignment::optimal_bandwidth, Aligner, Alignment, FMIndex, Reference};
use criterion::{black_box, criterion_group, criterion_main, Criterion};

fn create_test_aligner(size: usize) -> (Aligner, Vec<u8>) {
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
    let index = FMIndex::build(&reference);
    let aligner = Aligner::new(index, reference.as_slice().to_vec()).min_seed_len(10);
    (aligner, sequence)
}

fn bandwidth_calculation_benchmark(c: &mut Criterion) {
    let sizes = [50, 100, 150, 200, 500, 1000];

    c.bench_function_over_inputs(
        "optimal_bandwidth",
        |b, &size| {
            b.iter(|| black_box(optimal_bandwidth(size)));
        },
        sizes,
    );
}

fn alignment_benchmark(c: &mut Criterion) {
    let sizes = [100, 500, 1000];

    let mut group = c.benchmark_group("alignment");

    for size in sizes {
        let (aligner, _sequence) = create_test_aligner(size * 10);
        let read = vec![0u8, 1, 2, 3].repeat(size / 4);

        group.bench_function(format!("align_{}", size), |b| {
            b.iter(|| black_box(aligner.align_read(&read, None)));
        });
    }

    group.finish();
}

fn bandwidth_scaling_benchmark(c: &mut Criterion) {
    let query_lens = [50, 100, 200, 500];

    c.bench_function_over_inputs(
        "bandwidth_vs_query_len",
        |b, &query_len| {
            let bandwidth = optimal_bandwidth(query_len);
            let expected = (16 + query_len / 2).min(256);
            b.iter(|| {
                assert_eq!(black_box(bandwidth), expected);
            });
        },
        query_lens,
    );
}

fn alignment_with_different_bw(c: &mut Criterion) {
    // Test that different bandwidths don't break alignment
    let (aligner, _sequence) = create_test_aligner(10000);
    let read = vec![0u8, 1, 2, 3].repeat(25); // 100bp read

    c.bench_function("alignment_consistency", |b| {
        b.iter(|| {
            let result = aligner.align_read(&read, None);
            black_box(result.is_ok())
        });
    });
}

criterion_group!(
    benches,
    bandwidth_calculation_benchmark,
    alignment_benchmark,
    bandwidth_scaling_benchmark,
    alignment_with_different_bw
);
criterion_main!(benches);
