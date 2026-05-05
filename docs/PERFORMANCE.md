# Performance Benchmarking

Benchmark results comparing `bwa-rs` (Rust) against `bwa-mem` (C reference).

## Indexing Performance

### vs bwa (C reference)

| Size | bwa (C) | bwa-rs (Rust) | Ratio |
|------|---------|---------------|-------|
| 5KB | 3ms | 259ms | 0.01x |
| 50KB | 8ms | 255ms | 0.03x |
| 500KB | 50ms | 280ms | 0.18x |
| 5MB | 494ms | 470ms | 1.05x |
| **248MB (chr1)** | **236s** | **33s** | **7.0x** ✓ |

**Notes:**
- bwa-rs overhead dominates for small sequences due to O(n) SA construction initialization
- At real genome sizes (chr1+), bwa-rs **outperforms bwa** by 7x
- Both use multithreading internally (libsais/OpenMP)

### bwa-rs Indexing Rate

| Size | Time | Rate |
|------|------|------|
| 5KB | 0.27s | 19 KB/s |
| 50KB | 0.26s | 195 KB/s |
| 500KB | 0.28s | 1.8 MB/s |
| 5MB | 0.49s | 10 MB/s |
| 50MB | 2.49s | 20 MB/s |

O(n) SA-IS construction scales linearly with sequence size.

## Alignment Parameters

| Parameter | bwa-rs | bwa-mem | Notes |
|-----------|--------|---------|-------|
| match_score | 1 | 1 | Match bonus |
| mismatch_penalty | 4 | 4 | Mismatch penalty |
| gap_open | 6 | 5 | Gap opening penalty |
| gap_extend | 1 | 1 | Gap extension penalty |
| bandwidth | `min(256, 16 + len/2)` | dynamic | Banded DP width |
| min_seed_len | 10 | 19 | Minimum MEM seed |

**Bandwidth formula:** `min(256, 16 + query_len / 2)`
- Short reads (<500bp): bandwidth 26-266
- Long reads (>500bp): capped at 256

## SIMD Vectorization

Smith-Waterman alignment with runtime dispatch on x86:

| Hardware | Lanes | Vector Width |
|----------|-------|--------------|
| AVX-512 | 16 | 512-bit |
| AVX-2 | 8 | 256-bit |
| Scalar | 1 | 64-bit (ARM/fallback) |

**Detection:**
```rust
let config = get_simd_config();
// config.use_avx512 -> 16 lanes
// config.use_avx2 -> 8 lanes
// config.enabled -> false (scalar)
```

SIMD provides ~4-8x speedup for alignment vs scalar on x86.

## System Info

```
CPU: x86_64 Linux
Rust version: 1.89
bwa version: 0.7.19-r1273
```

## Running Benchmarks

```bash
# Install bwa for comparison
brew install bwa  # macOS
# or: apt install bwa  # Linux

# Build release
cargo build --release

# Run all benchmarks
cargo bench

# Compare against bwa
BWA_PATH=bwa cargo test test_compare_against_bwa_mem

# Time indexing
time ./target/release/bwa-rs index -r ref.fa -p ref

# Time alignment
time ./target/release/bwa-rs mem -R ref.fa -1 reads.fq -o out.sam
```

## Accuracy Verification

```bash
# Compare SAM outputs
./target/release/bwa-rs mem -R ref.fa -1 reads.fq -o rust.sam -k 10
bwa mem ref.fa reads.fq > c.sam

# Run integration tests
cargo test --all-features
```

**Note:** CIGAR normalization - bwa uses `M`, bwa-rs uses `=` for matches. Both are equivalent.

## Known Limitations

- Occ queries use O(log σ) wavelet tree (bwa uses O(1) sampling)
- MEM finding uses O(n log m) supermaximal (bwa uses O(n))
- MAPQ calculation may differ from bwa