# BWA-MEM Performance Optimization

**Project Status: ✅ Near Complete**

Pure Rust implementation targeting C BWA-MEM performance.

## Summary

| Status | Count |
|--------|-------|
| ✅ Done | 47 |
| 🔄 In Progress | 0 |
| ⬜ Pending | 1 |
| **Total** | **48** |

---

## Completed Phases

| Phase | Description | Tickets |
|-------|-------------|---------|
| Phase 1 | Core Implementation | T1-T24 |
| Phase 2 | Fast Suffix Array | T25, T26 |
| Phase 3 | Succinct Occ Table | T27, T28 |
| Phase 4 | SIMD Alignment | T29, T30 |
| Phase 5 | Memory Optimization | T31, T32 |
| Phase 6 | Parallelization | T33, T34 |
| Phase 7 | Code Cleanup | T35-T41 |
| Phase 8 | Optimization & Fixes | T42-T50 |

---

## Future Optimization Tickets

### T47: Benchmark and validate implementation ⬜

**Description:** Comprehensive benchmarking and correctness verification for the Rust BWA-MEM implementation on x86 hardware.

**Deliverables:**
- [x] Verify T43 supermaximal MEM correctness against legacy algorithm
  - ✅ Fixed MEM boundary bug (filtering positions that extend past reference)
  - ✅ Integration test `test_compare_against_bwa_mem` passes with C BWA-MEM
- [x] Compare alignment accuracy against C BWA-MEM reference
  - ✅ CIGAR matching works (M/= normalized)
  - ✅ Position matching works for non-repetitive sequences
  - ⚠️  Multiple alignment handling differs (MAPQ calculation)
- [x] Fill in PERFORMANCE.md with benchmark results
- [x] Benchmark wavelet tree occ queries vs sampling on large genomes
  - Benchmarks added in `benches/occ_benchmark.rs` and `benches/alignment_benchmark.rs`
- [x] Tune `optimal_bandwidth()` based on empirical results
- [ ] Document SIMD speedups on x86 (AVX2/AVX-512)

**Dependencies:** None (T44 and T45 completed)

**Known Issues:**
1. **min_seed_len default:** CLI defaults to 10 (updated from 19)

**Note:** Running on x86 compute node with AVX2/AVX-512 support.

---

### T44: Integrate wavelet tree for O(log σ) occ queries ✅

**Description:** Replace sampling-based OccTable with wavelet tree for guaranteed O(log σ) queries.

**Deliverables:**
- [x] Implement wavelet tree occurrence table in `src/occ/wavelet_tree.rs`
- [x] Add `CompactOccTable` in `src/compact.rs` using wavelet tree
- [x] Integrate wavelet tree into main `FMIndex` (see T46)
- [x] Remove or deprecate sampling-based `OccTable` (no longer exists - FMIndex uses CompactOccTable)
- [x] Benchmark improvement on large genomes (benchmarks added in `benches/`)

**Dependencies:** T42

**Note:** T46 handles the FM-Index integration work.

**Verification:**
- [x] FMIndex uses CompactOccTable with wavelet tree
- [x] All 14 integration tests pass
- [x] Benchmark infrastructure added for occ query comparison

---

### T45: Tune alignment parameters ✅

**Description:** Benchmark and tune bandwidth, scoring parameters for optimal alignment quality.

**Deliverables:**
- [x] Fill in PERFORMANCE.md with benchmarks
- [x] Compare alignment accuracy against C BWA-MEM (integration test available)
- [x] Tune `optimal_bandwidth()` based on empirical results

**Dependencies:** T42, T43, T44

**Verification:**
- [x] Added alignment parameter table to PERFORMANCE.md
- [x] Added alignment benchmarks (benches/alignment_benchmark.rs)
- [x] Current bandwidth formula: `min(256, 16 + query_len / 2)`
- [x] Integration test `test_compare_against_bwa_mem` available for accuracy verification

---

## Testing

The implementation uses **inline test data** in `tests/integration.rs`:

- `TEST_REFERENCE`: 144bp synthetic reference (`ACGT` repeated)
- `TEST_FASTQ`: 3 synthetic reads for testing

**Comparison tests** against real BWA-MEM:
- `test_compare_against_bwa_mem()` - runs `bwa mem` and compares results
- Uses `$BWA_PATH` env var (defaults to `bwa`)
- Creates temp reference/reads files, aligns both, compares SAM output

### Running BWA comparison tests

```bash
# Install bwa
brew install bwa  # macOS
# or: apt install bwa  # linux

# Run tests with comparison
BWA_PATH=bwa cargo test test_compare_against_bwa_mem

# Run all tests
cargo test
```

---

## Big-O Complexity

| Component | Algorithm | Complexity | C BWA-MEM |
|-----------|-----------|------------|-----------|
| Index Build | SA-IS | O(n) | O(n) |
| FM-Index Search | Backward search | O(m) | O(m) |
| MEM Finding | Supermaximal | O(n log m) | O(n) |
| Chaining | Greedy DP | O(k log k) | O(k log k) |
| Alignment | Banded SW | O(m × w) | O(m × w) |
| Occ Queries | Wavelet Tree | O(log σ) | O(1) |

**Goal:** Achieve O(n) MEM finding and O(1) occ queries for full parity.