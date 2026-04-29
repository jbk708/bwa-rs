# BWA-MEM Performance Optimization Plan

Target: Match or exceed C BWA-MEM performance for human-scale genomes.

---

## Performance Bottlenecks

| Component | Current | Target | Improvement |
|-----------|---------|--------|-------------|
| Suffix Array | O(n² log n) | O(n) | **100x** |
| Occ Table | O(32) scan | O(1) | **32x** |
| SW Extension | Scalar | SIMD | **4-8x** |
| Memory | 25GB | 5GB | **5x** |
| Index Time | Hours | ~5 min | **100x** |

---

## Optimization Roadmap

### Phase 1: Fast Suffix Array

**T25: Replace SA with SA-IS / libsais**
- Use induced sorting algorithm
- O(n) time, O(n) space
- Integrate via FFI or pure Rust port (libsais-rs)

**T26: Integer alphabet SA**
- Encode sequences as integers, not bytes
- Enable radix sort on suffix prefixes
- Remove string comparison overhead

---

### Phase 2: Succinct Occ Table

**T27: Wavelet Tree implementation**
- O(1) or O(log σ) rank queries
- Compressed storage
- Build from BWT in O(n)

**T28: RRR or SDArray**
- Succinct bitvectors for rank
- 2-3x compression over raw counts
- popcount-based occ queries

---

### Phase 3: SIMD Alignment

**T29: SIMD Smith-Waterman**
- Use `std::simd` or `wide` crate
- Process 16-32 bases per cycle
- AVX2/AVX-512 support

**T30: SIMD affine DP**
- Vectorized 3-matrix DP
- Banded traversal with SIMD

---

### Phase 4: Memory Optimization

**T31: Full memory mapping**
- Memory map entire index for 3GB genomes
- Never load fully into RAM
- Use `memmap2` for zero-copy access

**T32: Compact BWT encoding**
- 2 bits per character (already done)
- Bit-packed occurrence tables
- Streaming index construction

---

### Phase 5: Parallelization

**T33: Multi-threaded alignment**
- Rayon for parallel read processing
- Thread pool for batch alignment
- Lock-free FM-Index queries

**T34: Parallel seed finding**
- Partition query across threads
- Reduce per-read latency

---

## Implementation Order

```
T25 (SA-IS)          → 4-6 weeks
T27 (Wavelet Tree)   → 2-3 weeks  
T29 (SIMD SW)        → 2-4 weeks
T31 (Memory map)     → 1 week
T33 (Parallel)       → 2 weeks
T26, T28, T30, T32   → 4-6 weeks
T34                  → 1 week
```

---

## Verification

| Metric | C BWA-MEM | Target |
|--------|-----------|--------|
| Index 3GB | ~5 min | <10 min |
| Memory (index) | ~5GB | <10GB |
| Alignment throughput | ~10M reads/hr | >5M reads/hr |
| occ(c, k) query | O(1) | O(1) |

---

## Dependencies

| Crate | Purpose |
|-------|---------|
| `packed_simd` or `wide` | SIMD operations |
| `rayon` | Parallelization |
| `memmap2` | Memory mapping |
| `libsais-sys` (FFI) | Fast SA construction |

---

## Risk Mitigation

| Risk | Mitigation |
|------|------------|
| FFI complexity | Use pure Rust SA-IS if available |
| SIMD portability | Detect CPU features at runtime |
| Memory pressure | Aggressive streaming, no full loads |