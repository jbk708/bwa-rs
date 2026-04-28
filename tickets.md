# BWA-MEM Implementation Tickets

Pure Rust implementation of BWA-MEM for human-scale genomic alignment with paired-end support.

## Phase 1: Foundation

### T1: Project Setup
**Description:** Initialize Cargo project with minimal dependencies, define core types.
**Deliverables:**
- `Cargo.toml` with `clap`, `thiserror`, `memmap2`
- `lib.rs` with module exports
- Core structs: `Sequence`, `Quality`, `AlignmentRecord`

### T2: Reference Storage
**Description:** Implement 2-bit encoded genome storage, FASTA parsing helpers.
**Deliverables:**
- `Reference` struct with 2-bit encoding (A=00, C=01, G=10, T=11)
- Memory-efficient storage for 3GB genome
- `Reference::from_fasta()` parser

---

## Phase 2: FM-Index (Core Indexing)

### T3: Suffix Array
**Description:** Implement SA construction (SA-IS algorithm).
**Deliverables:**
- `SuffixArray` struct
- `SuffixArray::construct()` using induced copying
- Support for memory-mapped output

### T4: BWT Construction
**Description:** Build BWT from suffix array.
**Deliverables:**
- `BWT` struct
- `BWT::from_suffix_array()`
- FMD-index variant for FASTQ compatibility

### T5: FM-Index Queries
**Description:** Implement LF-mapping, backward search, OCC counts.
**Deliverables:**
- `FMIndex` struct
- `FMIndex::search(pattern) -> Vec<SAInterval>`
- Efficient OCC sampling

### T6: Index I/O
**Description:** Save/load FM-index to disk with memory mapping.
**Deliverables:**
- `FMIndex::save(path)`
- `FMIndex::load(path)` with memory-mapped files

---

## Phase 3: Seeding (MEM Finding)

### T7: Backward Search
**Description:** Extend FM-index search base-by-base.
**Deliverables:**
- `backward_search()` function
- Handle ambiguous bases (N)

### T8: MEM Discovery
**Description:** Find all Maximal Exact Matches recursively.
**Deliverables:**
- `find_mems(pattern, min_len) -> Vec<MEM>`
- Recursive backtracking implementation

### T9: Seed Filtering
**Description:** Filter overlapping and weak seeds.
**Deliverables:**
- `filter_seeds()` removing redundant seeds
- Seed scoring by uniqueness

---

## Phase 4: Alignment Extension

### T10: SW Matrix
**Description:** Initialize scoring matrix (match/mismatch/scoring).
**Deliverables:**
- `ScoringMatrix` config
- Default BWA-MEM scoring values

### T11: Banded SW
**Description:** Smith-Waterman with banded DP for speed.
**Deliverables:**
- `extend_seed_forward()`
- `extend_seed_backward()`
- Band width optimization

### T12: Affine Gaps
**Description:** Support gap open/extension penalties.
**Deliverables:**
- `AffineGapAlignment`
- Separate open vs extension penalties

---

## Phase 5: Chaining & Scoring

### T13: Seed Chaining
**Description:** DP chain seeds, score colinear paths.
**Deliverables:**
- `chain_seeds()`
- `ChainedSeed` struct with position and score

### T14: Best Alignment
**Description:** Select best chain, compute final position and CIGAR.
**Deliverables:**
- `AlignmentResult` struct
- CIGAR generation from alignment path

### T15: Mismatch Annotation
**Description:** Generate MD:Z tag for mismatches.
**Deliverables:**
- `mdz_string()` function
- Proper SAM MD:Z tag format

---

## Phase 6: Paired-End

### T16: Insert Size Model
**Description:** Track insert size distribution.
**Deliverables:**
- `InsertSizeDistribution` struct
- Mean and variance estimation

### T17: Pairing Logic
**Description:** Match read1/read2 by expected distance.
**Deliverables:**
- `pair_reads()`
- Proper FR orientation handling

### T18: Rescue
**Description:** Handle anomalous pairs, orphan rescue.
**Deliverables:**
- `rescue_unpaired()`
- Anomalous pair flagging

---

## Phase 7: SAM Output

### T19: SAM Record
**Description:** Build all SAM fields.
**Deliverables:**
- `SAMRecord::new()` with all 11 required fields
- FLAG calculation

### T20: SAM Writer
**Description:** Stream alignments to SAM file.
**Deliverables:**
- `SAMWriter` struct
- Header generation with @HD, @SQ

### T21: BAM Support
**Description:** Binary BAM output via hts-sys.
**Deliverables:**
- `BAMWriter` struct (optional)
- BGZF compression support

---

## Phase 8: CLI & Integration

### T22: CLI Interface
**Description:** Implement `bwa mem` equivalent CLI.
**Deliverables:**
- Index subcommand
- Align subcommand with all BWA-MEM options

### T23: Integration Tests
**Description:** End-to-end alignment tests.
**Deliverables:**
- Test on chr1 reference
- Verify SAM output correctness

### T24: README & Docs
**Description:** Usage documentation.
**Deliverables:**
- README.md with examples
- Cargo documentation

---

## Ticket Dependencies

```
T1 ──► T2 ──► T3 ──► T4 ──► T5 ──► T6 ──► T7 ──► T8 ──► T9 ──► T10 ──► T11 ──► T12 ──► T13 ──► T14 ──► T15 ──► T19 ──► T20 ──► T21 ──► T22 ──► T23 ──► T24
                                                                                                                                              │
                                                                                                                                              └──► T16 ──► T17 ──► T18 ──► (back to T22)
```

## Priority Order

1. T1, T2 - Foundation
2. T3, T4, T5, T6 - FM-Index (critical path)
3. T7, T8, T9 - Seeding
4. T10, T11, T12 - Alignment
5. T13, T14, T15 - Chaining
6. T16, T17, T18 - Paired-end
7. T19, T20, T21 - SAM output
8. T22, T23, T24 - CLI & docs