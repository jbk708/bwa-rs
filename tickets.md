# BWA-MEM Performance Optimization Tickets

Pure Rust implementation targeting C BWA-MEM performance.

## Summary

| Status | Count |
|--------|-------|
| 🟡 In Progress | 0 |
| ⬜ Pending | 7 |
| ✅ Done | 3 |
| **Total** | **10** |

---

## Phase 1: Fast Suffix Array

### T25: SA-IS Algorithm
**Status:** ✅ Done
**Description:** O(n) suffix array construction using sa-is crate.
**Deliverables:**
- [x] Integrated sa-is crate for O(n) SA construction
- [x] Streaming support for memory-bounded processing
- [x] All tests pass including large sequences

**Notes:** Using the `sa-is` crate for O(n) suffix array construction.
This replaces the sort-based O(n² log n) approach with true O(n) performance.

### T26: Integer Alphabet SA
**Status:** ✅ Done
**Description:** Radix sort on integer-encoded sequences.
**Deliverables:**
- [x] Encode sequences as u8/u16 integers
- [x] Radix sort suffixes using integer comparison
- [x] Eliminate string slice overhead

**Notes:** Using libsais-rs crate for O(n) integer alphabet SA construction. Provides `encode_sequence()`, `build_sa_integer()`, and `build_sa_i32()` functions. Radix-style sorting on compact integers eliminates string comparison overhead.

---

## Phase 2: Succinct Occ Table

### T27: Wavelet Tree
**Status:** ✅ Done
**Description:** O(1) rank queries with compressed storage.
**Deliverables:**
- [x] Wavelet tree implementation
- [x] O(1) occ(c, k) queries (using wavelet-matrix crate)
- [x] Build from BWT in O(n)

**Notes:** Using the `wavelet-matrix` crate for efficient rank queries. The wrapper handles edge cases (empty, single-element, all-same sequences).

### T28: RRR / SDArray
**Status:** ✅ Done
**Description:** Succinct bitvectors for rank queries.
**Deliverables:**
- [x] RRR bitvector implementation
- [x] popcount-based occ queries
- [x] 2-3x compression vs raw counts

**Notes:** Using the `succinct` crate's `Rank9` structure. Each character in the alphabet (A, C, G, T, N) has its own bitvector with O(1) rank queries. Compression ratio is ~10x for sparse character distributions, better for denser ones.

---

## Phase 3: SIMD Alignment

### T29: SIMD Smith-Waterman
**Status:** Pending
**Description:** Vectorized SW using AVX2/AVX-512.
**Deliverables:**
- [ ] `packed_simd` or `wide` integration
- [ ] Process 16-32 bases per cycle
- [ ] CPU feature detection at runtime

### T30: SIMD Affine DP
**Status:** Pending
**Description:** Vectorized 3-matrix affine gap alignment.
**Deliverables:**
- [ ] SIMD M/X/G matrices
- [ ] Banded traversal with vectorization
- [ ] Verify correctness vs scalar version

---

## Phase 4: Memory Optimization

### T31: Full Memory Mapping
**Status:** Pending
**Description:** Memory map entire index for zero-copy access.
**Deliverables:**
- [ ] `memmap2` for index files
- [ ] Never load index fully into RAM
- [ ] Lazy loading for BWT/SA/Occ

### T32: Compact Encoding
**Status:** Pending
**Description:** Bit-packed structures for minimal memory.
**Deliverables:**
- [ ] 2-bit BWT (done)
- [ ] Bit-packed occurrence tables
- [ ] Streaming index construction

---

## Phase 5: Parallelization

### T33: Multi-threaded Alignment
**Status:** Pending
**Description:** Parallel read processing with Rayon.
**Deliverables:**
- [ ] Thread pool for batch alignment
- [ ] Lock-free FM-Index queries
- [ ] Configurable thread count

### T34: Parallel Seeding
**Status:** Pending
**Description:** Partition query for multi-threaded MEM finding.
**Deliverables:**
- [ ] Divide long reads into chunks
- [ ] Parallel MEM discovery
- [ ] Merge and dedupe results

---

## Priority Order

1. **T25** - SA-IS (biggest bottleneck, 100x potential)
2. **T27** - Wavelet tree (32x occ improvement)
3. **T29** - SIMD SW (4-8x alignment speedup)
4. **T31** - Memory mapping (enables 3GB support)
5. **T33** - Parallelization (multi-core throughput)
6. T26, T28, T30, T32, T34 (refinements)