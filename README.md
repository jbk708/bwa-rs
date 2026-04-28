# bwa-rs

Pure Rust implementation of BWA-MEM for genomic alignment.

## Status

Early development. Currently implements:
- 2-bit encoded reference storage
- FM-Index with suffix array and BWT
- MEM (Maximal Exact Match) finding
- Basic alignment with CIGAR generation
- SAM format output
- Paired-end read handling

## Installation

```bash
cargo build --release
```

## Usage

### Index a reference

```bash
cargo run --release -- index -r genome.fa -p genome
```

### Align reads

```bash
cargo run --release -- mem -R genome.fa -1 reads_1.fq -2 reads_2.fq -o output.sam
```

## Architecture

```
src/
├── lib.rs           # Library exports
├── main.rs          # CLI entry point
├── error.rs         # Error types
├── types.rs         # Core types (Sequence, Cigar, MEM)
├── reference.rs     # Reference genome handling
├── fm_index.rs      # FM-Index, BWT, suffix array
├── seed.rs          # MEM finding
├── alignment.rs     # Smith-Waterman alignment
├── chaining.rs      # Seed chaining
├── sam.rs           # SAM format output
└── paired.rs        # Paired-end logic
```

## License

MIT