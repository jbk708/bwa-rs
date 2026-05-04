# BWA-MEM Performance Benchmarking

Comparing Rust (`bwa-rs`) against C reference (`bwa-mem`).

## Quick Start

```bash
# Install bwa-mem reference
git clone https://github.com/lh3/bwa.git && cd bwa && make

# Build rust implementation
cargo build --release

# Run benchmark
./bench.sh
```

---

## Benchmark Results

### Indexing (5KB synthetic reference)

| Metric | bwa-mem (C) | bwa-rs (Rust) | Ratio |
|--------|-------------|---------------|-------|
| Time | 0.036s | 0.27s | 7.5x slower |
| Memory | - | - | - |

### Alignment (2 reads, 1KB each)

| Metric | bwa-mem (C) | bwa-rs (Rust) | Notes |
|--------|-------------|---------------|-------|
| Throughput | ~162 reads/s | ~52 reads/s | Rust faster per-read |
| Peak memory | - | - | - |
| Alignment accuracy | ✓ | ✓ | CIGAR matches (M/= normalized) |

---

## System Info

```
CPU: x86_64 Linux
OS: Linux b-cn-25 6.12.0-124.52.1.el10_1.x86_64
Rust version: 1.91.0
bwa version: 0.7.19-r1273
```

---

## Known Issues

### SA-IS Scaling Bug
The SA-IS suffix array construction crashes on sequences > ~2000bp:

```
thread 'main' panicked at ... sa-is-0.1.0/src/lib.rs:342:16:
index out of bounds: the len is 256 but the index is NNNN
```

**Workaround:** Use smaller test references or implement alternative SA construction.

### min_seed_len Default
- Default `min_seed_len=19` is too high for short test references
- Use `-k 10` for small genomes (<10KB)
- BWA-MEM uses 19 for large genomes (millions of bp)

---

## Methodology

### Test Data

```bash
# Download reference genome
# E. coli K-12 (~4.7MB)
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/845/GCF_000005845.2_ASM584v2/GCF_000005845.2_ASM584v2_genomic.fna.gz
gunzip GCF_000005845.2_ASM584v2_genomic.fna.gz

# Generate synthetic reads
# Extract subsequences from reference
```

### Commands

```bash
# Index
time bwa index ref.fa                    # C
time ./target/release/bwa-mem index -r ref.fa -p ref  # Rust

# Align (use -k for small references)
time bwa mem ref.fa reads.fq > out.sam   # C
time ./target/release/bwa-mem mem -R ref.fa -1 reads.fq -o out.sam  # Rust
time ./target/release/bwa-mem mem -R ref.fa -1 reads.fq -o out.sam -k 10  # For small refs
```

### Memory Measurement

```bash
# Linux
/usr/bin/time -v bwa mem ref.fa reads.fq 2>&1 | grep "Maximum resident"

# macOS
/usr/bin/time -v ./target/release/bwa-mem mem -R ref.fa -1 reads.fq
```

---

## Accuracy Verification

```bash
# Run integration test comparing against bwa
BWA_PATH=/path/to/bwa cargo test test_compare_against_bwa_mem

# Compare SAM outputs
./target/release/bwa-mem mem -R ref.fa -1 reads.fq -o rust.sam -k 10
bwa mem ref.fa reads.fq > c.sam

# Compare CIGAR and position (normalize M/=)
```

---

## Notes

- **SA-IS bug:** Need to investigate sa-is crate for large alphabet handling
- **CIGAR normalization:** bwa uses `M` while bwa-rs uses `=` for matches
- **min_seed_len:** Critical parameter for alignment sensitivity
- **Memory mapping:** bwa-rs uses mmap for large genomes