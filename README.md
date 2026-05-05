# bwa-rs

[![CI](https://github.com/jbk708/bwa-rs/actions/workflows/ci.yml/badge.svg)](https://github.com/jbk708/bwa-rs/actions)
[![version](https://img.shields.io/badge/version-1.0.0-blue.svg)](Cargo.toml)

**Pure Rust implementation of BWA-MEM** for genomic alignment with paired-end support.

## Features

- **Pure Rust** - No C dependencies, compiles to a single binary
- **FM-Index** - O(n) suffix array, O(1) rank queries via wavelet tree
- **SIMD Alignment** - AVX2/AVX-512 Smith-Waterman (x86) / scalar fallback (ARM)
- **Memory Mapped** - Zero-copy index access for 3GB+ genomes
- **Parallel** - Multi-threaded alignment and seeding via Rayon
- **2-bit Encoding** - Memory-efficient storage
- **SAM Output** - Standard format with BGZF compression support
- **Streaming** - Process large FASTQ files with low memory

## Installation

```bash
cargo build --release
```

Binary: `target/release/bwa-rs`

## Usage

```bash
# Build index
bwa-rs index -r reference.fa -p ref

# Align single-end
bwa-rs mem -r reference.fa reads.fq -o output.sam

# Align paired-end
bwa-rs mem -r reference.fa -1 reads_1.fq -2 reads_2.fq -o output.sam
```

## Options

| Flag | Description | Default |
|------|-------------|---------|
| `-k` | Minimum seed length | 19 |
| `-t` | Thread count | auto |
| `-o` | Output file | stdout |

> **Tip:** Use `-k 10` for small references (<10KB). BWA-MEM uses 19 for large genomes.

## Architecture

```
src/
├── lib.rs              # Library exports
├── main.rs             # CLI entry point
├── error.rs            # Error types
├── types.rs            # Core types
├── reference.rs        # FASTA parsing, 2-bit encoding
├── sa.rs               # Suffix array construction
├── fm_index.rs         # FM-Index, BWT
├── mmap_index.rs       # Memory-mapped FM-Index
├── compact.rs          # Compact bit-packed structures
├── occ/                # Occurrence table
│   ├── mod.rs
│   └── wavelet_tree.rs # Wavelet tree for O(log σ) queries
├── mem_finder.rs       # Supermaximal MEM finding
├── seed.rs             # MEM discovery
├── chaining.rs         # Seed chaining
├── alignment.rs        # Smith-Waterman, affine gaps
├── simd_sw.rs          # SIMD Smith-Waterman (AVX2/AVX-512)
├── simd_affine.rs      # SIMD affine alignment
├── parallel.rs         # Multi-threaded alignment
├── parallel_seed.rs    # Chunked parallel seeding
├── sam.rs              # SAM records
├── bam.rs              # BAM output
├── paired.rs           # Pairing, rescue
├── fastq.rs            # FASTQ streaming
└── utils.rs            # Utilities
```

## Testing

See [testing.md](testing.md) for:
- Unit and integration tests
- BWA-MEM comparison tests
- Performance benchmarking

```bash
# Run tests
cargo test

# Compare against bwa mem (requires bwa installed)
brew install bwa
BWA_PATH=bwa cargo test test_compare_against_bwa_mem
```

## Documentation

- [docs/testing.md](docs/testing.md) - Testing guide
- [docs/PERFORMANCE.md](docs/PERFORMANCE.md) - Benchmark results
- [docs/RELEASE-1.0.0.md](docs/RELEASE-1.0.0.md) - Release checklist

## Benchmarks

Run performance benchmarks:

```bash
cargo bench
```

See [docs/PERFORMANCE.md](docs/PERFORMANCE.md) for detailed results.

## License

MIT