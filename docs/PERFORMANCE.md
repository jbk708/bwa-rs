# BWA-MEM Performance Benchmarking

Comparing Rust (`bwa-rs`) against C reference (`bwa-mem`).

## Quick Start

```bash
# Install bwa-mem reference
brew install bwa  # macOS
# or: sudo apt install bwa  # Linux

# Build rust implementation
cargo build --release

# Run benchmark
./bench.sh
```

---

## Benchmark Results

Fill in results from your benchmarks:

### Indexing

| Metric | bwa-mem (C) | bwa-rs (Rust) | Ratio |
|--------|-------------|---------------|-------|
| Time (3GB genome) | | | |
| Memory (3GB genome) | | | |
| Time (E. coli ~5MB) | | | |
| Memory (E. coli) | | | |

### Alignment

| Metric | bwa-mem (C) | bwa-rs (Rust) | Ratio |
|--------|-------------|---------------|-------|
| Throughput (reads/sec) | | | |
| Peak memory | | | |
| Alignment accuracy | | | |

---

## Methodology

### Test Data

```bash
# Download reference genome
# E. coli K-12
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/845/GCF_000005845.2_ASM584v2/GCF_000005845.2_ASM584v2_genomic.fna.gz
gunzip GCF_000005845.2_ASM584v2_genomic.fna.gz

# Generate synthetic reads (requires seqtk)
seqtk seq -a ref.fna -s 42 - 100000 > reads.fq
```

### Commands

```bash
# Index
time bwa index ref.fa                    # C
time ./target/release/bwa-rs index -r ref.fa -p ref  # Rust

# Align
time bwa mem ref.fa reads.fq > out.sam   # C
time ./target/release/bwa-rs mem -r ref.fa reads.fq -o out.sam  # Rust
```

### Memory Measurement

```bash
# Linux
/usr/bin/time -v bwa mem ref.fa reads.fq 2>&1 | grep "Maximum resident"

# macOS
/usr/bin/time -v ./target/release/bwa-rs mem -r ref.fa reads.fq
```

---

## Accuracy Verification

Verify alignment correctness:

```bash
cargo test test_compare_against_bwa_mem
```

Compare SAM outputs directly:

```bash
bwa mem ref.fa reads.fq > c_output.sam
./target/release/bwa-rs mem -r ref.fa reads.fq -o rust_output.sam

# Diff alignments (excluding timing/metadata)
diff <(grep -v "^@" c_output.sam | cut -f1-6 | sort) \
     <(grep -v "^@" rust_output.sam | cut -f1-6 | sort)
```

---

## System Info

Fill in your system specs:

```
CPU: 
RAM: 
OS: 
Rust version: 
bwa version: 
```

---

## Notes

- Hardware: Apple Silicon falls back to scalar SIMD (slower than x86 AVX)
- Memory mapping: bwa-rs uses mmap, may show different RSS vs bwa-mem
- Thread count: Use `-t` flag to control parallelism