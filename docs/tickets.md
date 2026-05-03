# BWA-MEM Performance Optimization

**Project Status: 🟡 In Progress**

Pure Rust implementation targeting C BWA-MEM performance.

## Summary

| Status | Count |
|--------|-------|
| ✅ Done | 16 |
| 🔄 In Progress | 2 |
| **Total** | **18** |

---

## Completed Phases

| Phase | Description | Tickets |
|-------|-------------|---------|
| Phase 1 | Fast Suffix Array | T25, T26 |
| Phase 2 | Succinct Occ Table | T27, T28 |
| Phase 3 | SIMD Alignment | T29, T30 |
| Phase 4 | Memory Optimization | T31, T32 |
| Phase 5 | Parallelization | T33, T34 |
| Phase 6 | Code Cleanup | T35 |

---

## Current Tickets

### T1: Consolidate duplicate encoding functions

**Status:** ✅ Done

**Description:** encode_sequence() is duplicated in sa.rs and reference.rs. Extract to a shared location.

**Deliverables:**
- [ ] Move encoding to utils.rs or types.rs
- [ ] Remove duplicates from sa.rs and reference.rs

---

### T2: Merge duplicate MD tag generation

**Status:** ✅ Done

**Description:** mdz_string() exists in both alignment.rs and types.rs producing identical output.

**Deliverables:**
- [ ] Keep one implementation (prefer types.rs since it lives on AlignmentResult)
- [ ] Remove duplicate from alignment.rs

---

### T3: Deduplicate Cigar compression

**Status:** ✅ Done

**Description:** compress_cigar() is copy-pasted in alignment.rs and simd_sw.rs.

**Deliverables:**
- [ ] Move to shared location (types.rs or utils.rs)
- [ ] Remove duplicate implementations

---

### T4: Remove dead SIMD code paths

**Status:** ✅ Done

**Description:** simd_affine.rs AVX2 function falls back to scalar at end. AVX512 is stub-only.

**Deliverables:**
- [ ] Remove incomplete affine_extend_forward_avx2 or complete traceback
- [ ] Remove empty affine_extend_forward_avx512 stub
- [ ] Keep #[allow(dead_code)] only if intentional for future work

---

### T5: Remove unused dependencies

**Status:** ✅ Done

**Description:** succinct crate imported but RrrBitvec never used by main code.

**Deliverables:**
- [ ] Remove succinct from Cargo.toml if not needed
- [ ] Or add #[allow(dead_code)] with justification comment

---

### T6: Optimize FM-index serialization

**Status:** 🟡 In Progress

**Description:** OccTable::read_from() reconstructs counts from BWT instead of reading pre-computed values.

**Deliverables:**
- [ ] Store counts directly in index file format
- [ ] Skip reconstruction on load
- [ ] Add test for large index load performance

---

### T7: Remove stored reference from FMIndex

**Status:** ⬜ Pending

**Description:** reference: Vec<u8> field is stored but never used (method is #[allow(dead_code)]).

**Deliverables:**
- [ ] Remove reference field from FMIndex struct
- [ ] Remove reference() accessor
- [ ] Update FMIndex::build() to not store reference

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