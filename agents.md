# Coding Standards for BWA-MEM Rust Implementation

## General Principles

1. **Zero-cost abstractions** - Prefer Rust idioms over C patterns; no unnecessary heap allocation
2. **Memory efficiency** - Human genome is 3GB; use 2-bit encoding, memory mapping, streaming
3. **No unnecessary dependencies** - Minimal `Cargo.toml` deps; implement from scratch where feasible
4. **Explicit error handling** - No panics in library code; use `Result<T, BwaError>` consistently
5. **Minimal inline comments** - Code should be self-documenting; comments explain *why*, not *what*

---

## Code Organization

### Module Structure

```
src/
├── lib.rs              # Public exports only
├── error.rs            # BwaError enum (thiserror)
├── types.rs            # Public types: Sequence, Quality, etc.
├── reference.rs        # Reference genome handling
├── fm_index/           # FM-Index submodules
│   ├── mod.rs
│   ├── suffix_array.rs
│   ├── bwt.rs
│   └── query.rs
├── seed.rs             # MEM finding
├── alignment.rs        # Smith-Waterman
├── chaining.rs         # Seed chaining
├── sam.rs              # SAM/BAM output
├── paired.rs           # Paired-end logic
└── cli.rs              # CLI (separate binary crate)
```

### File Naming
- One module per file, named after module
- Private helpers in `mod.rs` or `_inner.rs`

---

## Type Conventions

### Naming
- `CamelCase` for types: `Sequence`, `FMIndex`, `AlignmentResult`
- `snake_case` for functions and variables: `find_mems`, `insert_size`
- `SCREAMING_SNAKE_CASE` for constants: `DEFAULT_MIN_SEED_LEN`

### Structs
```rust
// Immutable structs prefer Clone over mutability
pub struct Sequence {
    name: String,
    bases: Vec<u8>,  // 2-bit encoded: A=0, C=1, G=2, T=3, N=4
}

pub struct FMIndex {
    bwt: BWT,
    sa: SuffixArray,
    occ: OccTable,
}
```

### Enums
```rust
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Strand {
    Forward,
    Reverse,
}

#[derive(Debug, Clone)]
pub enum BwaError {
    Io(std::io::Error),
    Parse(String),
    Index(String),
    Alignment(String),
}
```

---

## Function Design

### Signatures
```rust
// Prefer Result for fallible operations
pub fn find_mems(sequence: &Sequence, min_len: usize) -> Result<Vec<MEM>, BwaError>

// Use slices over Vec when borrowing
pub fn extend_seed(ref_seq: &[u8], query: &[u8], start: usize) -> AlignmentExt

// Use Iterator for streaming
pub fn iter_alignments(&self) -> impl Iterator<Item = Result<SAMRecord, BwaError>>
```

### No Unnecessary Cloning
```rust
// Bad
let seq_copy = sequence.clone();
process(seq_copy);

// Good - borrow when possible
process(&sequence);

// Or use indices for large data
fn find_in_range(data: &[u8], range: std::ops::Range<usize>) -> ...
```

---

## Error Handling

### Error Type
```rust
use thiserror::Error;

#[derive(Error, Debug)]
pub enum BwaError {
    #[error("I/O error: {0}")]
    Io(#[from] std::io::Error),

    #[error("Parse error: {0}")]
    Parse(String),

    #[error("Index error: {0}")]
    Index(String),

    #[error("Alignment error: {0}")]
    Alignment(String),
}
```

### Propagation
```rust
// Use ? operator
fn load_index(path: &Path) -> Result<FMIndex, BwaError> {
    let file = File::open(path)?;
    let bwt = BWT::read(&file)?;
    Ok(FMIndex { bwt, ... })
}
```

---

## Performance Guidelines

### Memory Mapping
```rust
// Use memmap2 for large files
use memmap2::Mmap;

pub struct MemoryMappedIndex {
    _mmap: Mmap,
    // Access via slices into mmap
}
```

### Bounded Loops
```rust
// Avoid O(n²) in hot paths
for seed in seeds.iter() {
    // Linear operations only
}
```

### SIMD (Future)
- Consider `std::simd` for batch operations when stable
- Profile before optimizing

---

## Testing Standards

### Unit Tests
```rust
#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_bwt_construction() {
        let seq = b"ACGTACGT";
        let bwt = BWT::from_sequence(seq);
        assert_eq!(bwt[0], b'T'); // Last char + sentinel
    }
}
```

### Property Tests (with `proptest`)
```rust
#[cfg(test)]
mod property_tests {
    use proptest::prelude::*;

    proptest! {
        #[test]
        fn test_mem_finding_never_exceeds_bounds(seq: Vec<u8>) {
            // seq is 2-bit encoded
            let mems = find_mems(&seq, 10)?;
            for mem in mems {
                assert!(mem.ref_end <= seq.len());
            }
        }
    }
}
```

### Integration Tests
```rust
// tests/integration.rs
#[test]
fn test_chr1_alignment() {
    let ref = Reference::from_fasta("test data/chr1.fa").unwrap();
    let index = FMIndex::build(&ref);
    let reads = load_fastq("test data/reads.fq");
    // Verify output SAM
}
```

---

## Documentation

### Module Docs
```rust
//! FM-Index implementation for fast substring search.
//!
//! Uses the Ferragina-Manzini index with O(m) search time
//! where m is the pattern length.
//!
//! # Example
//! ```
//! let index = FMIndex::build(&reference);
!```

### Function Docs
```rust
/// Find all Maximal Exact Matches (MEMs) in the reference.
///
/// A MEM is a maximal substring match between query and reference.
/// Returns positions sorted by start coordinate.
///
/// # Arguments
/// * `query` - The query sequence to search
/// * `min_len` - Minimum MEM length (BWA default: 19)
///
/// # Errors
/// Returns `BwaError::Index` if FM-index is not built.
```

---

## CLI Design

### Use Clap Derive
```rust
use clap::{Parser, Subcommand};

#[derive(Parser)]
#[command(name = "bwa-rs")]
#[command(about = "Pure Rust BWA-MEM implementation")]
pub enum Cli {
    /// Index a reference genome
    Index(IndexArgs),
    /// Align reads to reference
    Mem(MemArgs),
}

#[derive(Args)]
pub struct MemArgs {
    /// Reference database
    #[arg(short = 'R')]
    pub reference: PathBuf,

    /// FASTQ file for read 1
    #[arg(short = '1')]
    pub read1: PathBuf,

    /// FASTQ file for read 2
    #[arg(short = '2')]
    pub read2: Option<PathBuf>,

    /// Output file (SAM format)
    #[arg(short = 'o', default_value = "-")]
    pub output: PathBuf,

    /// Minimum seed length
    #[arg(short = 'k', default_value = "19")]
    pub min_seed_len: u32,
}
```

---

## Naming Conventions for Biology

| Concept | Variable/Type | Notes |
|---------|---------------|-------|
| Reference sequence | `ref_seq`, `reference` | 0-indexed |
| Read sequence | `read`, `query` | From FASTQ |
| Position | `pos`, `start`, `end` | 0-indexed in code, 1-indexed in SAM |
| CIGAR | `cigar` | String like "10M2I5M" |
| Quality | `qual`, `quality` | Phred scaled |
| Mate read | `mate`, `read2` | Paired-end partner |

---

## Commit Messages

Format: `T{N}: {Brief description}`

Examples:
- `T5: Implement FM-index backward search`
- `T12: Add affine gap penalty support`
- `T20: Stream SAM output to file`

---

## Code Review Checklist

- [ ] No panics in library code
- [ ] All fallible operations return `Result`
- [ ] Memory-mapped large structures
- [ ] No unnecessary clones (profile if unsure)
- [ ] Tests pass: `cargo test`
- [ ] Tests cover edge cases: empty input, N bases, reverse complement
- [ ] Documentation compiles: `cargo doc`
- [ ] Clippy clean: `cargo clippy`