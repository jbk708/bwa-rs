//! Suffix Array construction module using SA-IS algorithm.
//!
//! This module provides O(n) suffix array construction via the sa-is crate
//! and integer alphabet support via libsais-rs for efficient radix sorting.

use std::io::{Read, Write};

pub use crate::utils::{encode_sequence, encode_sequence_u16};

/// Build suffix array using integer alphabet via libsais-rs.
/// Returns indices into the original sequence.
///
/// Uses the 64-bit libsais entry point so references larger than the i32 SA
/// limit (~2.147 Gbp) can be indexed — e.g. the full human genome (3.1 Gbp).
pub fn build_sa_integer(seq: &[u8]) -> Vec<u64> {
    if seq.is_empty() {
        return vec![];
    }
    let n = seq.len();
    let mut sa = vec![0i64; n];
    // fs=0 for standard SA construction, freq=None for no frequency output.
    // Use the multithreaded (rayon-backed) entry point; libsais only parallelizes
    // for large inputs (n >= 65536) and falls back to serial otherwise.
    let threads = rayon::current_num_threads().max(1) as i64;
    libsais_rs::libsais64::libsais64_omp(seq, &mut sa, 0, None, threads);
    sa.into_iter().map(|x| x as u64).collect()
}

/// Build suffix array from pre-encoded integer sequence.
pub fn build_sa_from_encoded(encoded: &[u8]) -> Vec<u64> {
    if encoded.is_empty() {
        return vec![];
    }
    let n = encoded.len();
    let mut sa = vec![0i64; n];
    libsais_rs::libsais64(encoded, &mut sa, 0, None);
    sa.into_iter().map(|x| x as u64).collect()
}

/// Build suffix array using integer-encoded alphabet (64-bit).
pub fn build_sa_i32(encoded: &[i32], alphabet_size: i32) -> Vec<u64> {
    if encoded.is_empty() {
        return vec![];
    }
    let n = encoded.len();
    let mut sa = vec![0i64; n];
    let mut encoded_i64: Vec<i64> = encoded.iter().map(|&x| x as i64).collect();
    // k = alphabet size, fs = 0
    libsais_rs::libsais64::libsais64_int(&mut encoded_i64, &mut sa, alphabet_size as i64, 0);
    sa.into_iter().map(|x| x as u64).collect()
}

#[derive(Clone, Debug)]
pub struct SuffixArray {
    pub sa: Vec<u64>,
    len: usize,
}

impl SuffixArray {
    pub fn build(sequence: &[u8]) -> Self {
        if sequence.is_empty() {
            return Self { sa: vec![], len: 0 };
        }

        // Use libsais-rs for fast suffix array construction
        // Fixes T48: sa-is crate crashes on sequences > ~2000bp
        let sa = build_sa_integer(sequence);

        Self {
            sa,
            len: sequence.len(),
        }
    }

    /// Create a new SuffixArray with pre-computed SA values and explicit length
    /// Used by FM-index which needs n+1 entries including sentinel
    pub fn with_len(sa: Vec<u64>, len: usize) -> Self {
        Self { sa, len }
    }

    pub fn get(&self, idx: usize) -> Option<u64> {
        self.sa.get(idx).copied()
    }

    pub fn len(&self) -> usize {
        self.len
    }

    pub fn is_empty(&self) -> bool {
        self.len == 0
    }

    pub fn as_slice(&self) -> &[u64] {
        &self.sa
    }

    pub fn write_to(&self, writer: &mut impl Write) -> std::io::Result<()> {
        for &idx in &self.sa {
            writer.write_all(&idx.to_le_bytes())?;
        }
        Ok(())
    }

    pub fn read_from(reader: &mut impl Read, count: usize) -> std::io::Result<Self> {
        let mut sa = Vec::with_capacity(count);
        for _ in 0..count {
            let mut bytes = [0u8; 8];
            reader.read_exact(&mut bytes)?;
            sa.push(u64::from_le_bytes(bytes));
        }
        Ok(Self { sa, len: count })
    }
}

pub fn build_sa_streaming<'a>(
    sequence: &'a [u8],
    chunk_size: usize,
) -> impl Iterator<Item = Vec<u64>> + 'a
where
    &'a [u8]: 'a,
{
    sequence
        .chunks(chunk_size)
        .map(|chunk| SuffixArray::build(chunk).sa)
}

pub fn build_sa_with_sentinel(sequence: &[u8]) -> Vec<u64> {
    let n = sequence.len();
    let mut padded = Vec::with_capacity(n + 1);
    padded.extend_from_slice(sequence);
    padded.push(4);

    let sa = SuffixArray::build(&padded);
    let mut result: Vec<u64> = sa.into_iter().filter(|&x| x != n as u64).collect();
    result.sort();
    result
}

impl IntoIterator for SuffixArray {
    type Item = u64;
    type IntoIter = std::vec::IntoIter<u64>;

    fn into_iter(self) -> Self::IntoIter {
        self.sa.into_iter()
    }
}

impl<'a> IntoIterator for &'a SuffixArray {
    type Item = &'a u64;
    type IntoIter = std::slice::Iter<'a, u64>;

    fn into_iter(self) -> Self::IntoIter {
        self.sa.iter()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn verify_sa(seq: &[u8], sa: &[u64]) -> bool {
        if sa.len() != seq.len() {
            return false;
        }

        let mut seen = vec![false; seq.len()];
        for &pos in sa {
            if pos as usize >= seq.len() {
                return false;
            }
            seen[pos as usize] = true;
        }

        for seen_pos in &seen {
            if !seen_pos {
                return false;
            }
        }

        for i in 0..sa.len().saturating_sub(1) {
            let pos1 = sa[i] as usize;
            let pos2 = sa[i + 1] as usize;
            if seq[pos1..].cmp(&seq[pos2..]) != std::cmp::Ordering::Less {
                return false;
            }
        }

        true
    }

    #[test]
    fn test_sa_small() {
        let seq = b"A";
        let sa = SuffixArray::build(seq);
        let vals: Vec<u64> = sa.into_iter().collect();
        assert!(verify_sa(seq, &vals));
    }

    #[test]
    fn test_sa_two_chars() {
        let seq = b"AC";
        let sa = SuffixArray::build(seq);
        let vals: Vec<u64> = sa.into_iter().collect();
        assert!(verify_sa(seq, &vals), "SA for AC: {:?}", vals);
    }

    #[test]
    fn test_sa_acgt() {
        let seq = b"ACGT";
        let sa = SuffixArray::build(seq);
        let vals: Vec<u64> = sa.into_iter().collect();
        assert!(verify_sa(seq, &vals), "SA for ACGT: {:?}", vals);
    }

    #[test]
    fn test_sa_repeated() {
        let seq = b"AAAA";
        let sa = SuffixArray::build(seq);
        let vals: Vec<u64> = sa.into_iter().collect();
        assert!(verify_sa(seq, &vals));
    }

    #[test]
    fn test_sa_medium() {
        let seq = b"AACGAACGG";
        let sa = SuffixArray::build(seq);
        let vals: Vec<u64> = sa.into_iter().collect();
        assert!(verify_sa(seq, &vals), "SA for AACGAACGG: {:?}", vals);
    }

    #[test]
    fn test_sa_acgtacgt() {
        let seq = b"ACGTACGT";
        let sa = SuffixArray::build(seq);
        let vals: Vec<u64> = sa.into_iter().collect();
        assert!(verify_sa(seq, &vals), "SA for ACGTACGT: {:?}", vals);
    }

    #[test]
    fn test_sa_random() {
        let seq = b"GATCGATCGA";
        let sa = SuffixArray::build(seq);
        let vals: Vec<u64> = sa.into_iter().collect();
        assert!(verify_sa(seq, &vals));
    }

    #[test]
    fn test_sa_longer() {
        let seq = b"ABABABABA";
        let sa = SuffixArray::build(seq);
        let vals: Vec<u64> = sa.into_iter().collect();
        assert!(verify_sa(seq, &vals));
    }

    #[test]
    fn test_sa_empty() {
        let seq: &[u8] = &[];
        let sa = SuffixArray::build(seq);
        assert!(sa.is_empty());
    }

    #[test]
    fn test_sa_with_n() {
        let seq = b"ACGNACGT";
        let sa = SuffixArray::build(seq);
        let vals: Vec<u64> = sa.into_iter().collect();
        assert!(verify_sa(seq, &vals));
    }

    #[test]
    fn test_build_sa_with_sentinel() {
        let seq = b"ACGT";
        let sa = build_sa_with_sentinel(seq);
        assert_eq!(sa.len(), 4);
        assert!(verify_sa(seq, &sa));
    }

    #[test]
    fn test_streaming() {
        let seq = b"ACGTACGTACGT";
        let chunks: Vec<_> = build_sa_streaming(seq, 4).collect();
        assert_eq!(chunks.len(), 3);

        for chunk in &chunks {
            assert!(!chunk.is_empty());
        }
    }

    #[test]
    fn test_large_sequence() {
        let seq: Vec<u8> = b"ACGT".repeat(100);
        let sa = SuffixArray::build(&seq);
        let vals: Vec<u64> = sa.into_iter().collect();
        assert!(verify_sa(&seq, &vals), "SA for large sequence failed");
    }

    #[test]
    fn test_large_sequence_t48() {
        // T48: SA-IS crash on sequences > ~2000bp
        let seq: Vec<u8> = b"ACGT".repeat(1250); // 5000bp
        let sa = SuffixArray::build(&seq);
        let vals: Vec<u64> = sa.into_iter().collect();
        assert_eq!(vals.len(), 5000, "SA length should match sequence length");
        assert!(
            verify_sa(&seq, &vals),
            "SA for 5000bp sequence should be valid"
        );
    }

    #[test]
    fn test_all_same_char() {
        let seq = vec![0u8; 100];
        let sa = SuffixArray::build(&seq);
        let vals: Vec<u64> = sa.into_iter().collect();
        assert!(verify_sa(&seq, &vals), "SA for all same char failed");
    }

    #[test]
    fn test_sa_consistency() {
        for len in [1, 2, 5, 10, 20, 50, 100, 500, 1000] {
            let seq: Vec<u8> = (0..len).map(|i| (i % 5) as u8).collect();
            let sa = SuffixArray::build(&seq);
            let vals: Vec<u64> = sa.into_iter().collect();
            assert!(verify_sa(&seq, &vals), "SA for len={} should be valid", len);
        }
    }

    // Integer alphabet encoding tests
    #[test]
    fn test_encode_sequence_basic() {
        let seq = b"ACGT";
        let encoded = encode_sequence(seq);
        assert_eq!(encoded, &[0, 1, 2, 3]);
    }

    #[test]
    fn test_encode_sequence_with_n() {
        let seq = b"ACGNACGT";
        let encoded = encode_sequence(seq);
        assert_eq!(encoded, &[0, 1, 2, 4, 0, 1, 2, 3]);
    }

    #[test]
    fn test_encode_sequence_lowercase() {
        let seq = b"acgt";
        let encoded = encode_sequence(seq);
        assert_eq!(encoded, &[0, 1, 2, 3]);
    }

    #[test]
    fn test_encode_sequence_empty() {
        let seq: &[u8] = b"";
        let encoded = encode_sequence(seq);
        assert!(encoded.is_empty());
    }

    #[test]
    fn test_encode_sequence_unknown_char() {
        let seq = b"ACXGT";
        let encoded = encode_sequence(seq);
        assert_eq!(encoded, &[0, 1, 4, 2, 3]);
    }

    #[test]
    fn test_build_sa_integer_small() {
        let seq = b"ACGT";
        let sa = build_sa_integer(seq);
        assert_eq!(sa.len(), 4);
        assert!(verify_sa(seq, &sa));
    }

    #[test]
    fn test_build_sa_integer_medium() {
        let seq = b"AACGAACGG";
        let sa = build_sa_integer(seq);
        assert_eq!(sa.len(), 9);
        assert!(verify_sa(seq, &sa));
    }

    #[test]
    fn test_build_sa_integer_repeated() {
        let seq = b"AAAA";
        let sa = build_sa_integer(seq);
        assert_eq!(sa.len(), 4);
        assert!(verify_sa(seq, &sa));
    }

    #[test]
    fn test_build_sa_from_encoded() {
        let encoded = &[0u8, 1, 2, 3]; // ACGT
        let sa = build_sa_from_encoded(encoded);
        assert_eq!(sa.len(), 4);
        // Verify by checking that encoded[sa[i]..] < encoded[sa[i+1]..]
        for i in 0..sa.len().saturating_sub(1) {
            let idx1 = sa[i] as usize;
            let idx2 = sa[i + 1] as usize;
            assert!(encoded[idx1..].cmp(&encoded[idx2..]) == std::cmp::Ordering::Less);
        }
    }

    #[test]
    fn test_build_sa_integer_consistency() {
        for len in [1, 2, 5, 10, 20, 50, 100] {
            let seq: Vec<u8> = (0..len).map(|i| [b'A', b'C', b'G', b'T'][i % 4]).collect();
            let sa = build_sa_integer(&seq);
            assert_eq!(sa.len(), len, "SA length mismatch for len={}", len);
            assert!(
                verify_sa(&seq, &sa),
                "SA verification failed for len={}",
                len
            );
        }
    }

    #[test]
    fn test_build_sa_integer_empty() {
        let seq: &[u8] = b"";
        let sa = build_sa_integer(seq);
        assert!(sa.is_empty());
    }

    #[test]
    fn test_encode_sequence_u16() {
        let seq = b"ACGT";
        let encoded = encode_sequence_u16(seq);
        assert_eq!(encoded, &[0u16, 1, 2, 3]);
    }

    #[test]
    fn test_build_sa_i32() {
        let encoded = vec![0i32, 1, 2, 3];
        let sa = build_sa_i32(&encoded, 5);
        assert_eq!(sa.len(), 4);
        // Verify suffix array is correct for integer sequence
        for i in 0..sa.len().saturating_sub(1) {
            let idx1 = sa[i] as usize;
            let idx2 = sa[i + 1] as usize;
            assert!(encoded[idx1..].cmp(&encoded[idx2..]) == std::cmp::Ordering::Less);
        }
    }

    #[test]
    fn test_integer_vs_byte_sa_equivalence() {
        // Both methods should produce the same suffix array
        let seq = b"ACGTACGT";
        let sa_byte = SuffixArray::build(seq);
        let sa_int = build_sa_integer(seq);

        // Both should have the same length
        assert_eq!(sa_byte.len(), sa_int.len());

        // Both should be valid suffix arrays
        let vals_byte: Vec<u64> = sa_byte.into_iter().collect();
        assert!(verify_sa(seq, &vals_byte));
        assert!(verify_sa(seq, &sa_int));
    }

    #[test]
    fn test_large_sequence_integer() {
        let seq: Vec<u8> = b"ACGT".repeat(100);
        let sa = build_sa_integer(&seq);
        assert_eq!(sa.len(), 400);
        assert!(verify_sa(&seq, &sa));
    }
}
