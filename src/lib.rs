//! Pure Rust BWA-MEM implementation for genomic alignment.
//!
//! This crate provides a complete implementation of the BWA-MEM algorithm
//! for aligning sequencing reads to a reference genome.
//!
//! # Features
//!
//! - Pure Rust implementation (no C dependencies)
//! - FM-Index for fast substring search
//! - MEM (Maximal Exact Match) seeding
//! - Smith-Waterman alignment extension
//! - Paired-end read alignment
//! - SAM format output
//!
//! # Example
//!
//! ```ignore
//! use bwa_mem::{Reference, FMIndex, Aligner};
//!
//! // Load reference and build index
//! let reference = Reference::from_fasta("genome.fa")?;
//! let index = FMIndex::build(&reference);
//!
//! // Align reads
//! let aligner = Aligner::new(index);
//! // Note: queries must be 2-bit encoded sequences
//! let results = aligner.align_read(&[0, 1, 2, 3, 0, 1, 2, 3], None)?;
//! ```

pub mod alignment;
pub mod bam;
pub mod chaining;
pub mod compact;
pub mod error;
pub mod fastq;
pub mod fm_index;
pub mod mmap_index;
pub mod occ;
pub mod paired;
pub mod reference;
pub mod sa;
pub mod sam;
pub mod seed;
pub mod simd_sw;
pub mod types;

pub use alignment::Aligner;
pub use compact::{BitPackedBWT, CompactOccTable, StreamingFMIndex, StreamingFMIndexBuilder};
pub use error::BwaError;
pub use fm_index::FMIndex;
pub use reference::Reference;
pub use sa::SuffixArray;
pub use simd_sw::{extend_forward_simd, get_simd_config, nw_score, SimdConfig};
pub use types::*;
