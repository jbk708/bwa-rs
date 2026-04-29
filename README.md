# bwa-rs

Pure Rust implementation of BWA-MEM for genomic alignment with paired-end support.

## Features

- **Pure Rust** - No C dependencies, compiles to a single binary
- **2-bit encoded** - Memory-efficient storage for 3GB genomes
- **FM-Index** - Fast substring search with BWT and suffix array
- **MEM seeding** - Maximal Exact Matches for sensitive alignment
- **Smith-Waterman** - Banded alignment with affine gap penalties
- **Paired-end** - FR orientation, insert size estimation, orphan rescue
- **SAM/BAM** - Standard format output with BGZF compression
- **Streaming** - Process large FASTQ files with low memory footprint

## Installation

```bash
cargo build --release
```

The binary will be at `target/release/bwa-rs`.

## Usage

### Build index

```bash
bwa-rs index -r reference.fa -p ref
```

### Align single-end reads

```bash
bwa-rs mem -R reference.fa -1 reads.fq -o output.sam
```

### Align paired-end reads

```bash
bwa-rs mem -R reference.fa -1 reads_1.fq -2 reads_2.fq -o output.sam
```

### Options

| Flag | Description | Default |
|------|-------------|---------|
| `-k` | Minimum seed length | 19 |
| `-o` | Output file (SAM) | stdout |

## Architecture

```
src/
├── lib.rs           # Library exports
├── main.rs          # CLI entry point
├── error.rs         # Error types (BwaError)
├── types.rs         # Core types (Sequence, Cigar, AlignmentResult)
├── reference.rs     # FASTA parsing, 2-bit encoding
├── fm_index.rs      # FM-Index, BWT, suffix array, save/load
├── fastq.rs         # FASTQ streaming parser
├── seed.rs          # MEM finding, filtering
├── alignment.rs     # Smith-Waterman, affine gaps, chaining
├── chaining.rs      # Seed chaining DP
├── sam.rs           # SAM record construction
├── bam.rs           # BAM/BCF output with BGZF
└── paired.rs        # Insert size model, pairing, rescue
```

## Testing

```bash
cargo test
```

Integration tests verify alignment correctness against reference BWA-MEM.

## License

MIT