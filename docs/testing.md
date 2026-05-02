# Testing Guide

## Running Tests

```bash
# All tests
cargo test

# With output
cargo test -- --nocapture

# Specific test
cargo test test_compare_against_bwa_mem
```

## Test Structure

### Unit Tests

Inline in each source file covering:
- Sequence encoding (2-bit)
- CIGAR parsing/validation
- SAM record formatting
- Error handling edge cases

### Integration Tests

Located in `tests/integration.rs`:

| Test | Description |
|------|-------------|
| `test_complete_alignment_pipeline` | Full workflow: index → align → SAM |
| `test_sam_header_format` | @HD, @SQ, @PG header generation |
| `test_perfect_match` | Exact match detection |
| `test_mismatch_position` | Mismatch handling |
| `test_multiple_alignments_same_read` | Multi-MEM handling |
| `test_reference_sequence_encoding` | 2-bit encoding verification |
| `test_mapq_calculation` | Mapping quality scoring |
| `test_cigar_validation` | CIGAR string validation |
| `test_unmapped_read` | Unmapped read handling |
| `test_sam_record_fields` | Field completeness |
| `test_chr1_scaled_reference` | Large reference scaling |
| `test_compare_against_bwa_mem` | **BWA-MEM comparison** |
| `test_sam_format_complete` | Full SAM output validation |

## Compare Against BWA-MEM

The `test_compare_against_bwa_mem` test aligns the same reads with both implementations and compares:
- QNAME (read name)
- RNAME (reference name)
- CIGAR (alignment)
- Position (1-indexed)

### Prerequisites

```bash
# Install bwa
brew install bwa  # macOS
# or
sudo apt install bwa  # Ubuntu/Debian

# Verify
bwa 2>&1 | head -3
```

### Run Comparison

```bash
# Default bwa path
BWA_PATH=bwa cargo test test_compare_against_bwa_mem

# Custom path
BWA_PATH=/usr/local/bin/bwa cargo test test_compare_against_bwa_mem
```

If BWA is not installed, the test prints `SKIP: bwa not available` and passes.

## Test Data

Tests use inline synthetic data:

```rust
// 144bp reference: ACGT repeated
const TEST_REFERENCE: &str = ">test_ref\nACGTACGT...";

// 3 synthetic reads with varying complexity
const TEST_FASTQ: &str = "@read1\nACGTACGT...\n...";
```

## Performance Testing

For meaningful benchmarks, use real genomic data:

```bash
# Download E. coli genome (~5MB)
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/845/GCF_000005845.2_ASM584v2/GCF_000005845.2_ASM584v2_genomic.fna.gz
gunzip GCF_000005845.2_ASM584v2_genomic.fna.gz

# Build index
cargo run --release -- index -r GCF_000005845.2_ASM584v2_genomic.fna.gz -p ecoli

# Generate synthetic reads (requires seqtk or similar)
# seqtk seq -a ecoli.fna | head -10000 > reads.fa

# Align
time cargo run --release -- mem -r ecoli.fna reads.fa -o output.sam
```

## CI/Development

```bash
# Full check
cargo fmt --check
cargo clippy
cargo test
cargo build --release
```

## Coverage

Run with cargo-llvm-cov for coverage reports:

```bash
cargo install cargo-llvm-cov
cargo llvm-cov --lcov --output-path lcov.info
```