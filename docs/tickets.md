# BWA-MEM Performance Optimization

**Project Status: ✅ Complete**

Pure Rust implementation targeting C BWA-MEM performance.

## Summary

| Status | Count |
|--------|-------|
| ✅ Done | 35 |
| 🔄 In Progress | 0 |
| ⬜ Pending | 4 |
| **Total** | **39** |

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

### T42: Implement true SIMD Smith-Waterman

**Description:** Vectorize the Smith-Waterman DP using AVX2/AVX-512 instead of scalar fallback.

**Deliverables:**
- [ ] Implement 8-lane AVX2 SW using `wide` crate
- [ ] Implement 16-lane AVX-512 SW for supported hardware
- [ ] Verify correctness against scalar baseline
- [ ] Add CPU feature detection with runtime dispatch

**Dependencies:** None

---

### T43: Implement supermaximal MEM finding

**Description:** Replace recursive binary search with BWA's supermaximal MEM algorithm for O(n) expected time.

**Deliverables:**
- [ ] Implement Z-array based supermaximal MEM discovery
- [ ] Remove `find_mems_recursive` in favor of linear algorithm
- [ ] Verify correctness on test sequences

**Dependencies:** None

---

### T44: Integrate wavelet tree for O(1) occ queries

**Description:** Replace sampling-based OccTable with wavelet tree for guaranteed constant-time queries.

**Deliverables:**
- [ ] Integrate `src/occ/wavelet_tree.rs` into FM-Index
- [ ] Remove sampling-based occ table
- [ ] Benchmark improvement on large genomes

**Dependencies:** T42

---

### T45: Tune alignment parameters

**Description:** Benchmark and tune bandwidth, scoring parameters for optimal alignment quality.

**Deliverables:**
- [ ] Fill in PERFORMANCE.md with benchmarks
- [ ] Compare alignment accuracy against C BWA-MEM
- [ ] Tune `optimal_bandwidth()` based on empirical results

**Dependencies:** T42, T43, T44

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
| MEM Finding | Current: O(n log n) | O(n log n) | O(n) |
| Chaining | Greedy DP | O(k log k) | O(k log k) |
| Alignment | Banded SW | O(m × w) | O(m × w) |
| Occ Queries | Sampling | O(w) avg | O(1) |

**Goal:** Achieve O(n) MEM finding and O(1) occ queries for full parity.