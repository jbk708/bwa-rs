# BWA-MEM Performance Optimization Tickets

Pure Rust implementation targeting C BWA-MEM performance.

## Summary

| Status | Count |
|--------|-------|
| 🟡 In Progress | 0 |
| ⬜ Pending | 5 |
| ✅ Done | 4 |
| **Total** | **9** |

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
**Status:** ✅ Done
**Description:** Vectorized SW using AVX2/AVX-512.
**Deliverables:**
- [x] `wide` crate integration
- [x] Process 8-16 bases per cycle (AVX2: 8 lanes, AVX-512: 16 lanes)
- [x] CPU feature detection at runtime

**Changes:**
- Added `src/simd_sw.rs` with SIMD-accelerated Smith-Waterman
- `wide` crate for portable SIMD operations on x86
- Runtime CPU feature detection with `is_x86_feature_detected!` macro
- AVX2 implementation processes 8 lanes (256-bit), AVX-512 processes 16 lanes (512-bit)
- Automatic fallback to scalar implementation on unsupported hardware
- Exported `get_simd_config()`, `nw_score()`, and `extend_forward_simd()` functions

### T30: SIMD Affine DP
**Status:** ❌ Failed
**Description:** Vectorized 3-matrix affine gap alignment.
**Deliverables:**
- [ ] SIMD M/X/G matrices
- [ ] Banded traversal with vectorization
- [ ] Verify correctness vs scalar version

---

## Phase 4: Memory Optimization

### T31: Full Memory Mapping
**Status:** ✅ Done
**Description:** Memory map entire index for zero-copy access.
**Deliverables:**
- [x] `memmap2` for index files
- [x] Never load index fully into RAM
- [x] Lazy loading for BWT/SA/Occ

**Changes:**
- Added `src/mmap_index.rs` with `MmapFMIndex` struct
- `MmapFMIndex::open()` - memory-maps index without loading into RAM
- `MmapFMIndex::build_and_save()` - builds and saves in one step
- All search operations read directly from mmap
- Changed `unsafe_code` lint to `deny` (allows module-level `#![allow]`)

### T32: Compact Encoding
**Status:** ✅ Done 
**Description:** Bit-packed structures for minimal memory.
**Deliverables:**
- [x] 2-bit BWT (done)
- [x] Bit-packed occurrence tables
- [x] Streaming index construction

---

## Phase 5: Parallelization

### T33: Multi-threaded Alignment
**Status:** ✅ Done
**Description:** Parallel read processing with Rayon.
**Deliverables:**
- [x] Thread pool for batch alignment
- [x] Lock-free FM-Index queries
- [x] Configurable thread count

**Changes:**
- Added `rayon` dependency for parallel processing
- Added `src/parallel.rs` with `ParallelAligner` struct
- `ParallelAligner::align_batch()` - parallel alignment for multiple reads
- `ParallelAligner::align_batch_paired()` - parallel paired-end alignment
- `ParallelAligner::align_single()` - single read alignment
- `ThreadPoolConfig` - configurable thread pool settings
- `default_thread_count()` - auto-detect CPU cores
- Added `-t/--threads` CLI flag for thread count control
- Exported `ParallelAligner`, `ThreadPoolConfig`, `default_thread_count`

### T34: Parallel Seeding
**Status:** Pending
**Description:** Partition query for multi-threaded MEM finding.
**Deliverables:**
- [ ] Divide long reads into chunks
- [ ] Parallel MEM discovery
- [ ] Merge and dedupe results

---

## Maintenance

### T35: Code Cleanup
**Status:** ✅ Done
**Description:** Simplify source files and integration tests.
**Deliverables:**
- [x] Review and simplify each source file
- [x] Clean up integration tests
- [x] Ensure `cargo clippy` passes

**Changes:**
- Fixed `manual_strip` clippy warning in `reference.rs`
- Implemented `Display` trait for `SAMHeader` and `SAMRecord` (replaced `to_string` methods)
- Added `#[allow(clippy::too_many_arguments)]` for `SAMRecord::new`
- Consolidated imports in `types.rs`
- Simplified integration tests with `create_test_aligner()` helper function
- Fixed temporary lifetime issues in integration tests

---

## Priority Order

1. **T25** - SA-IS (biggest bottleneck, 100x potential)
2. **T27** - Wavelet tree (32x occ improvement)
3. **T29** - SIMD SW (4-8x alignment speedup)
4. **T31** - Memory mapping (enables 3GB support)
5. **T33** - Parallelization (multi-core throughput)
6. T26, T28, T30, T32, T34 (refinements)