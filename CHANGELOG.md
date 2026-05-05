# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added
- Initial release preparation

## [0.1.0] - 2026-05-05

### Added

#### Core Algorithm
- **FM-Index** with backward search (O(m) per pattern)
- **Suffix Array construction** using libsais-rs (O(n) algorithm)
- **Wavelet Tree** for O(log σ) occurrence queries
- **MEM finding** via supermaximal algorithm
- **Smith-Waterman alignment** with banded dynamic programming
- **Affine gap penalties** with 3-matrix dynamic programming
- **Paired-end alignment** with FR orientation and rescue

#### Performance
- **SIMD vectorization** for Smith-Waterman (AVX2/AVX-512 on x86)
- **SIMD affine alignment** for gap extension
- **Multi-threaded alignment** via Rayon thread pool
- **Parallel MEM seeding** with chunked query processing
- **Memory-mapped FM-Index** for zero-copy access to 3GB+ genomes
- **Compact bit-packed encoding** (2-bit BWT, wavelet tree occ table)

#### I/O and Format Support
- **SAM output** with full format compliance
- **BAM output** with BGZF compression
- **FASTQ streaming** for low-memory processing of large files
- **FASTA reference parsing** with 2-bit encoding

#### CLI Interface
- `index` command for building FM-index
- `mem` command for read alignment
- Configurable seed length (`-k`), threads (`-t`), and output (`-o`)

#### Testing
- 232 tests including unit, integration, and BWA-MEM comparison tests
- Benchmark suite with criterion
- FM-index correctness verification

#### Documentation
- README with features, architecture, and usage examples
- Testing guide (`docs/testing.md`)
- Performance benchmarks (`docs/PERFORMANCE.md`)
- Project tracking (`docs/tickets.md`)

### Fixed
- FM-index sentinel handling for long references
- CIGAR query length calculation in affine alignment
- MEM boundary bug (filtering positions that extend past reference)
- SA-IS crash on large sequences
- Unused variable warnings
- Code formatting issues

### Changed
- Minimum Rust version: 1.89
- Default minimum seed length: 10 (from 19)
- FM-Index now stores only reference data, not entire index
- CIGAR operators normalized (M/= for matches)

### Performance Improvements
- 7x faster indexing on chr1 compared to reference implementation
- Optimized bandwidth formula: `min(256, 16 + query_len / 2)`

[Unreleased]: https://github.com/jbk708/bwa-rs/compare/v0.1.0...HEAD
[0.1.0]: https://github.com/jbk708/bwa-rs/releases/tag/v0.1.0