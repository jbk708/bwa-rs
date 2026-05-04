# BWA-MEM Performance Optimization

**Project Status: 🔄 In Progress**

Pure Rust implementation targeting C BWA-MEM performance.

## Summary

| Status | Count |
|--------|-------|
| ✅ Done | 38 |
| 🔄 In Progress | 1 |
| ⬜ Pending | 1 |
| **Total** | **40** |

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

---

## Future Optimization Tickets

### T47: Benchmark and validate implementation 🔄

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
- [ ] Benchmark wavelet tree occ queries vs sampling on large genomes
  - ⚠️  SA-IS crashes on sequences > ~2000bp (sa-is crate bug)
- [ ] Tune `optimal_bandwidth()` based on empirical results
- [ ] Document SIMD speedups on x86 (AVX2/AVX-512)

**Dependencies:** T43, T44, T45, T46

**Known Issues:**
1. **SA-IS scaling bug:** `sa-is` crate panics on large sequences with "index out of bounds: the len is 256 but the index is NNNN"
2. **min_seed_len default:** CLI defaults to 10 (updated from 19)

**Note:** Running on x86 compute node with AVX2/AVX-512 support.

---

### T48: Investigate and fix SA-IS crash on large sequences ⬜

**Description:** The `sa-is` crate panics when building suffix arrays for sequences larger than ~2000bp.

**Error:**
```
thread 'main' panicked at sa-is-0.1.0/src/lib.rs:342:16:
index out of bounds: the len is 256 but the index is NNNN
```

**Reproducer:**
```rust
let ref = Reference::parse_fasta(">test\nACGT...
"); // ~5000bp
let index = FMIndex::build(&ref); // crashes
```

**Deliverables:**
- [ ] Identify root cause (alphabet size handling? buffer overflow?)
- [ ] Fix sa-is integration or switch to alternative SA construction
- [ ] Verify with E. coli genome (~4.7MB)
- [ ] Add regression test for large sequences

**Dependencies:** None

**Potential solutions:**
1. Fix sa-is crate configuration (alphabet size parameter)
2. Use alternative SA construction library (sais, sais-lite, libsais)
3. Implement fallback to O(n²) construction for small inputs

---

### T42: Implement true SIMD Smith-Waterman ✅

**Description:** Vectorize the Smith-Waterman DP using AVX2/AVX-512 instead of scalar fallback.

**Deliverables:**
- [x] Implement 8-lane AVX2 SW using intrinsics
- [x] Implement 16-lane AVX-512 SW for supported hardware
- [x] Verify correctness against scalar baseline
- [x] Add CPU feature detection with runtime dispatch

**Dependencies:** None

---

### T43: Implement supermaximal MEM finding ✅

**Description:** Replace recursive binary search with BWA's supermaximal MEM algorithm for O(n) expected time.

**Deliverables:**
- [x] Implement supermaximal MEM discovery algorithm
- [x] Remove `find_mems_recursive` in favor of binary search
- [x] Add comprehensive tests for MEM finding
- [x] Integrate `find_supermaximal_mems` into `Aligner` pipeline
- [ ] Verify correctness against legacy algorithm

**Dependencies:** None

**Implementation:** Added `src/mem_finder.rs` with binary search based MEM finding. The `find_mems()` function in `src/seed.rs` delegates to `find_supermaximal_mems()`.

---

### T44: Integrate wavelet tree for O(log σ) occ queries

**Description:** Replace sampling-based OccTable with wavelet tree for guaranteed O(log σ) queries.

**Deliverables:**
- [x] Implement wavelet tree occurrence table in `src/occ/wavelet_tree.rs`
- [x] Add `CompactOccTable` in `src/compact.rs` using wavelet tree
- [x] Integrate wavelet tree into main `FMIndex` (see T46)
- [ ] Remove or deprecate sampling-based `OccTable`
- [ ] Benchmark improvement on large genomes

**Dependencies:** T42

**Note:** T46 handles the FM-Index integration work.

---

### T45: Tune alignment parameters

**Description:** Benchmark and tune bandwidth, scoring parameters for optimal alignment quality.

**Deliverables:**
- [ ] Fill in PERFORMANCE.md with benchmarks
- [ ] Compare alignment accuracy against C BWA-MEM
- [ ] Tune `optimal_bandwidth()` based on empirical results

**Dependencies:** T42, T43, T44

---

### T46: Integrate wavelet tree into FMIndex ✅

**Description:** Replace the sampling-based `OccTable` in `src/fm_index.rs` with `CompactOccTable` that uses wavelet tree for O(log σ) occ queries.

**Deliverables:**
- [x] Replace `OccTable` with `CompactOccTable` in `FMIndex` struct
- [x] Update `FMIndex::build()` to use `CompactOccTable::from_bwt()`
- [x] Update serialization (`save`/`load`) to be compatible with wavelet tree
- [x] Ensure all existing tests pass
- [ ] Add benchmark comparing sampling vs wavelet tree occ queries

**Dependencies:** T44

**Implementation:**
- `src/occ/wavelet_tree.rs` - added `Clone` and `Debug` implementations with `original_data` field
- `src/compact.rs` - `CompactOccTable` derives `Clone` and `Debug`
- `src/fm_index.rs` - uses `CompactOccTable` instead of `OccTable`, updated to version 3

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

| Phase | Algorithm | Complexity | C BWA-MEM |
|-------|-----------|------------|-----------|
| Index Build | SA-IS | O(n) | O(n) |
| FM-Index Search | Backward search | O(m) | O(m) |
| MEM Finding | Supermaximal | O(n log m) | O(n) |
| Chaining | Greedy DP | O(k log k) | O(k log k) |
| Alignment | Banded SW | O(m × w) | O(m × w) |
| Occ Queries | Wavelet Tree | O(log σ) | O(1) |

**Goal:** Achieve O(n) MEM finding and O(1) occ queries for full parity.