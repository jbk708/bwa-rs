# BWA-MEM 1.0.0 Release Checklist

**Target:** Pure Rust BWA-MEM implementation ready for production use

---

## Pre-Release Requirements

### Documentation
- [x] README.md - Complete with features, installation, usage
- [x] docs/testing.md - Testing guide
- [x] docs/PERFORMANCE.md - Benchmark results
- [x] CHANGELOG.md - Release notes
- [x] Cargo.toml - Complete metadata

### Code Quality
- [x] `cargo test` passes (231+ tests)
- [x] `cargo clippy` clean - no warnings
- [x] `cargo fmt` compliant
- [x] No `unimplemented!()` or `todo!()` in code
- [x] No TODO/FIXME/XXX/HACK comments
- [x] No panics in library code - all fallibles return `Result<BwaError>`

### CI/CD
- [x] GitHub Actions CI workflow
- [x] Tests on Rust 1.89 and stable
- [x] Clippy checks
- [x] Format checks

### Core Functionality
- [x] FM-Index with backward search (O(m) per pattern)
- [x] Suffix array construction (libsais-rs, O(n))
- [x] Wavelet tree for occ queries (O(log σ))
- [x] MEM finding (supermaximal algorithm)
- [x] Smith-Waterman alignment (banded DP)
- [x] Affine gap penalties (3-matrix DP)
- [x] SIMD vectorization (AVX2/AVX-512)
- [x] Memory-mapped index for 3GB+ genomes
- [x] Multi-threaded alignment (Rayon)
- [x] SAM/BAM output with BGZF compression
- [x] Paired-end alignment (FR orientation)
- [x] MD:Z tag generation

---

## Pre-Release Tasks

### Version & Metadata
- [ ] Update `Cargo.toml` version from `0.1.0` to `1.0.0`
- [x] Update repository URL from `USER/bwa-mem` to actual repo
- [ ] Update rust-version if needed (currently 1.89)

### API Stability
- [ ] Review public API for any breaking changes needed
- [ ] Ensure error types are comprehensive
- [ ] Document any deprecated items

### Testing
- [ ] Run full integration test suite with bwa comparison
- [ ] Test on large genomes (chr1, human reference)
- [ ] Verify SAM output format compliance
- [ ] Test paired-end scenarios

### Performance
- [ ] Verify PERFORMANCE.md benchmarks are current
- [ ] Confirm 7x speedup on chr1 indexing is reproducible
- [ ] Test memory usage on 3GB+ genomes

### Documentation
- [ ] Add usage examples to README
- [ ] Document CLI options completely
- [ ] Add API documentation with `cargo doc`
- [x] CHANGELOG.md - Created with release notes

### Security
- [ ] Review unsafe code blocks
- [ ] Audit dependencies for vulnerabilities
- [ ] Verify no sensitive data in tests

---

## Release Steps

1. **Branch:** Create `release/1.0.0` branch
2. **Update:** Bump version, update docs
3. **Test:** Run full test suite including bwa comparison
4. **Review:** Code review focused on API stability
5. **Tag:** Create `v1.0.0` tag
6. **Publish:** `cargo publish --dry-run` then `cargo publish`
7. **Announce:** Update GitHub releases, notify users

---

## Post-Release Considerations

### Future Enhancements (After 1.0.0)
- T44: Remove/deprecate sampling-based OccTable
- T45: Tune alignment parameters with benchmarking
- T47: Complete benchmarking validation

### Known Limitations
- Occ queries use O(log σ) wavelet tree (vs O(1) in bwa)
- MEM finding uses O(n log m) supermaximal (vs O(n) in bwa)
- MAPQ calculation may differ from bwa

---

## Verification Commands

```bash
# Full test suite
cargo test --all-features

# Compare against bwa
BWA_PATH=bwa cargo test test_compare_against_bwa_mem

# Clippy
cargo clippy -- -D warnings

# Format
cargo fmt --check

# Documentation
cargo doc --no-deps

# Build release
cargo build --release
```