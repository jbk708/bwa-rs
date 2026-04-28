# BWA-MEM Implementation Tickets

Pure Rust implementation of BWA-MEM for human-scale genomic alignment with paired-end support.

## Phase 1: Foundation

### ✅ T1: Project Setup
**Status:** Complete
**Description:** Initialize Cargo project with minimal dependencies, define core types.
**Deliverables:**
- ✅ `Cargo.toml` with `clap`, `thiserror`, `memmap2`
- ✅ `lib.rs` with module exports
- ✅ Core structs: `Sequence`, `Quality`, `AlignmentRecord`

### ✅ T2: Reference Storage
**Status:** Complete
**Description:** 2-bit encoded genome storage, FASTA parsing helpers.
**Deliverables:**
- ✅ `Reference` struct with 2-bit encoding (A=0, C=1, G=2, T=3, N=4)
- ✅ Memory-efficient storage for 3GB genome
- ✅ `Reference::from_fasta()` parser

---

## Phase 2: FM-Index

### ✅ T3: Suffix Array
**Status:** Complete
**Description:** SA construction using O(n log n) sort-based method.
**Deliverables:**
- ✅ `SuffixArray` struct
- ✅ `SuffixArray::build()` implementation

### ✅ T4: BWT Construction
**Status:** Complete
**Description:** Build BWT from suffix array.
**Deliverables:**
- ✅ `BWT` struct
- ✅ `BWT::from_sa()`

### ✅ T5: FM-Index Queries
**Status:** Complete
**Description:** LF-mapping, backward search, OCC counts.
**Deliverables:**
- ✅ `FMIndex` struct
- ✅ `FMIndex::search()` - backward search
- ✅ `FMIndex::find_all()` - get all positions
- ✅ `OccTable::occ()` - efficient rank queries

### ⬜ T6: Index I/O
**Status:** Pending
**Description:** Save/load FM-index to disk with memory mapping.
**Deliverables:**
- [ ] `FMIndex::save(path)`
- [ ] `FMIndex::load(path)` with memory-mapped files
- [ ] Index file format with magic header

---

## Phase 3: Seeding

### ✅ T7: Backward Search
**Status:** Complete
**Description:** Extend FM-index search base-by-base.
**Deliverables:**
- ✅ Recursive backward search implementation

### ✅ T8: MEM Discovery
**Status:** Complete
**Description:** Find all Maximal Exact Matches recursively.
**Deliverables:**
- ✅ `find_mems()` with configurable `min_seed_len`
- ✅ Recursive backtracking implementation

### ✅ T9: Seed Filtering
**Status:** Complete
**Description:** Filter overlapping and weak seeds.
**Deliverables:**
- ✅ `filter_mems()` removing redundant seeds

---

## Phase 4: Alignment Extension

### ✅ T10: SW Matrix
**Status:** Complete
**Description:** Initialize scoring matrix (match/mismatch/scoring).
**Deliverables:**
- ✅ `Scoring` struct with default BWA-MEM values

### ⬜ T11: Banded SW
**Status:** Pending
**Description:** Smith-Waterman with banded DP for speed.
**Deliverables:**
- [ ] `extend_seed_forward()`
- [ ] `extend_seed_backward()`
- [ ] Band width optimization

### ⬜ T12: Affine Gaps
**Status:** Pending
**Description:** Support gap open/extension penalties.
**Deliverables:**
- [ ] `AffineGapAlignment` implementation
- [ ] Separate open vs extension penalties

---

## Phase 5: Chaining & Scoring

### ✅ T13: Seed Chaining
**Status:** Complete
**Description:** Chain seeds into coherent alignment path.
**Deliverables:**
- ✅ `chain_seeds()` DP implementation
- ✅ `ChainedSeed` struct with scores

### 🟡 T14: Best Alignment
**Status:** Partial
**Description:** Select best chain, compute position and CIGAR.
**Deliverables:**
- [x] `AlignmentResult` struct
- [ ] Full CIGAR generation from alignment path
- [ ] MAPQ calculation from alignment score

### ⬜ T15: Mismatch Annotation
**Status:** Pending
**Description:** Generate MD:Z tag for mismatches.
**Deliverables:**
- [ ] `mdz_string()` function
- [ ] Proper SAM MD:Z tag format

---

## Phase 6: Paired-End

### ✅ T16: Insert Size Model
**Status:** Complete
**Description:** Track insert size distribution.
**Deliverables:**
- ✅ `InsertSizeDistribution` struct
- ✅ Online mean/variance estimation

### ✅ T17: Pairing Logic
**Status:** Complete
**Description:** Match read1/read2 by expected distance.
**Deliverables:**
- ✅ `pair_reads()` with FR orientation handling
- ✅ Proper pair flag setting

### ✅ T18: Rescue
**Status:** Complete
**Description:** Handle anomalous pairs, orphan rescue.
**Deliverables:**
- ✅ `rescue_unpaired()` for orphan reads
- ✅ Anomalous pair detection

---

## Phase 7: SAM Output

### ✅ T19: SAM Record
**Status:** Complete
**Description:** Build all SAM fields.
**Deliverables:**
- ✅ `SAMRecord` struct with all 11 required fields
- ✅ FLAG calculation for paired/mapped

### ✅ T20: SAM Writer
**Status:** Complete
**Description:** Stream alignments to SAM file.
**Deliverables:**
- ✅ `SAMWriter` struct
- ✅ Header generation (@HD, @SQ, @PG)

### ⬜ T21: BAM Support
**Status:** Pending
**Description:** Binary BAM output.
**Deliverables:**
- [ ] `BAMWriter` struct
- [ ] BGZF compression via hts-sys
- [ ] BAM header handling

---

## Phase 8: CLI & Integration

### 🟡 T22: CLI Interface
**Status:** Partial
**Description:** Implement `bwa mem` equivalent CLI.
**Deliverables:**
- [x] Basic CLI structure with `index` and `mem` subcommands
- [ ] FASTQ parsing for reads
- [ ] Full BWA-MEM option support (-k, -w, -d, -T, etc.)
- [ ] Progress reporting for large files

### ⬜ T23: Integration Tests
**Status:** Pending
**Description:** End-to-end alignment tests.
**Deliverables:**
- [ ] Test on chr1 reference
- [ ] Compare output to reference BWA-MEM
- [ ] Verify SAM format correctness

### ✅ T24: README & Docs
**Status:** Complete
**Description:** Usage documentation.
**Deliverables:**
- ✅ README.md with examples
- ✅ Module documentation

---

## Summary

| Status | Count |
|--------|-------|
| ✅ Complete | 14 |
| 🟡 Partial | 3 |
| ⬜ Pending | 7 |
| **Total** | **24** |

---

## Priority Order for Remaining Work

1. **T6** - Index I/O (enables saving/loading indexes)
2. **T11, T12** - Proper Smith-Waterman alignment
3. **T14, T15** - CIGAR generation, MD:Z tag
4. **T22** - Full CLI with FASTQ parsing
5. **T21** - BAM output
6. **T23** - Integration tests