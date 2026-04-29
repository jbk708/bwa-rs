//! Suffix Array construction module using SA-IS algorithm.
//!
//! This module provides O(n) suffix array construction via the sa-is crate
//! and integer alphabet support via libsais-rs for efficient radix sorting.

use std::io::{Read, Write};

/// Encode a DNA sequence as integers: A=0, C=1, G=2, T=3, N=4
/// This enables efficient integer-based suffix array construction.
pub fn encode_sequence(seq: &[u8]) -> Vec<u8> {
    seq.iter()
        .map(|&b| match b.to_ascii_uppercase() {
            b'A' => 0,
            b'C' => 1,
            b'G' => 2,
            b'T' => 3,
            b'N' => 4,
            _ => 4,
        })
        .collect()
}

/// Encode a DNA sequence as u16 integers for larger alphabets.
pub fn encode_sequence_u16(seq: &[u8]) -> Vec<u16> {
    seq.iter()
        .map(|&b| match b.to_ascii_uppercase() {
            b'A' => 0,
            b'C' => 1,
            b'G' => 2,
            b'T' => 3,
            b'N' => 4,
            _ => 4,
        })
        .collect()
}

/// Build suffix array using integer alphabet via libsais-rs.
/// Returns indices into the original sequence.
pub fn build_sa_integer(seq: &[u8]) -> Vec<u32> {
    if seq.is_empty() {
        return vec![];
    }
    let encoded = encode_sequence(seq);
    let n = encoded.len();
    let mut sa = vec![0i32; n];
    // k=5 for 5-symbol alphabet (A,C,G,T,N), fs=0, freq=None
    libsais_rs::libsais(&encoded, &mut sa, 0, None);
    sa.into_iter().map(|x| x as u32).collect()
}

/// Build suffix array from pre-encoded integer sequence.
pub fn build_sa_from_encoded(encoded: &[u8]) -> Vec<u32> {
    if encoded.is_empty() {
        return vec![];
    }
    let n = encoded.len();
    let mut sa = vec![0i32; n];
    libsais_rs::libsais(encoded, &mut sa, 0, None);
    sa.into_iter().map(|x| x as u32).collect()
}

/// Build suffix array using i32-encoded integer alphabet.
pub fn build_sa_i32(encoded: &[i32], alphabet_size: i32) -> Vec<u32> {
    if encoded.is_empty() {
        return vec![];
    }
    let n = encoded.len();
    let mut sa = vec![0i32; n];
    let mut encoded_i32 = encoded.to_vec();
    // k = alphabet size, fs = 0
    libsais_rs::libsais_int(&mut encoded_i32, &mut sa, alphabet_size, 0);
    sa.into_iter().map(|x| x as u32).collect()
}

#[derive(Clone, Debug)]
pub struct SuffixArray {
    pub sa: Vec<u32>,
    len: usize,
}

impl SuffixArray {
    pub fn build(sequence: &[u8]) -> Self {
        if sequence.is_empty() {
            return Self { sa: vec![], len: 0 };
        }

        // Use sa-is for O(n) suffix array construction
        // key_bound = 256 for byte alphabet
        let sa_usize: Vec<usize> = sa_is::make_suffix_array(sequence, 256);
        let sa: Vec<u32> = sa_usize.into_iter().map(|x| x as u32).collect();

        Self {
            sa,
            len: sequence.len(),
        }
    }

    pub fn get(&self, idx: usize) -> Option<u32> {
        self.sa.get(idx).copied()
    }

    pub fn len(&self) -> usize {
        self.len
    }

    pub fn is_empty(&self) -> bool {
        self.len == 0
    }

    pub fn as_slice(&self) -> &[u32] {
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
            let mut bytes = [0u8; 4];
            reader.read_exact(&mut bytes)?;
            sa.push(u32::from_le_bytes(bytes));
        }
        Ok(Self { sa, len: count })
    }
}

pub fn build_sa_streaming<'a>(
    sequence: &'a [u8],
    chunk_size: usize,
) -> impl Iterator<Item = Vec<u32>> + 'a
where
    &'a [u8]: 'a,
{
    sequence
        .chunks(chunk_size)
        .map(|chunk| SuffixArray::build(chunk).sa)
}

pub fn build_sa_with_sentinel(sequence: &[u8]) -> Vec<u32> {
    let n = sequence.len();
    let mut padded = Vec::with_capacity(n + 1);
    padded.extend_from_slice(sequence);
    padded.push(4);

    let sa = SuffixArray::build(&padded);
    let mut result: Vec<u32> = sa.into_iter().filter(|&x| x != n as u32).collect();
    result.sort();
    result
}

impl IntoIterator for SuffixArray {
    type Item = u32;
    type IntoIter = std::vec::IntoIter<u32>;

    fn into_iter(self) -> Self::IntoIter {
        self.sa.into_iter()
    }
}

impl<'a> IntoIterator for &'a SuffixArray {
    type Item = &'a u32;
    type IntoIter = std::slice::Iter<'a, u32>;

    fn into_iter(self) -> Self::IntoIter {
        self.sa.iter()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn verify_sa(seq: &[u8], sa: &[u32]) -> bool {
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
        let vals: Vec<u32> = sa.into_iter().collect();
        assert!(verify_sa(seq, &vals));
    }

    #[test]
    fn test_sa_two_chars() {
        let seq = b"AC";
        let sa = SuffixArray::build(seq);
        let vals: Vec<u32> = sa.into_iter().collect();
        assert!(verify_sa(seq, &vals), "SA for AC: {:?}", vals);
    }

    #[test]
    fn test_sa_acgt() {
        let seq = b"ACGT";
        let sa = SuffixArray::build(seq);
        let vals: Vec<u32> = sa.into_iter().collect();
        assert!(verify_sa(seq, &vals), "SA for ACGT: {:?}", vals);
    }

    #[test]
    fn test_sa_repeated() {
        let seq = b"AAAA";
        let sa = SuffixArray::build(seq);
        let vals: Vec<u32> = sa.into_iter().collect();
        assert!(verify_sa(seq, &vals));
    }

    #[test]
    fn test_sa_medium() {
        let seq = b"AACGAACGG";
        let sa = SuffixArray::build(seq);
        let vals: Vec<u32> = sa.into_iter().collect();
        assert!(verify_sa(seq, &vals), "SA for AACGAACGG: {:?}", vals);
    }

    #[test]
    fn test_sa_acgtacgt() {
        let seq = b"ACGTACGT";
        let sa = SuffixArray::build(seq);
        let vals: Vec<u32> = sa.into_iter().collect();
        assert!(verify_sa(seq, &vals), "SA for ACGTACGT: {:?}", vals);
    }

    #[test]
    fn test_sa_random() {
        let seq = b"GATCGATCGA";
        let sa = SuffixArray::build(seq);
        let vals: Vec<u32> = sa.into_iter().collect();
        assert!(verify_sa(seq, &vals));
    }

    #[test]
    fn test_sa_longer() {
        let seq = b"ABABABABA";
        let sa = SuffixArray::build(seq);
        let vals: Vec<u32> = sa.into_iter().collect();
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
        let vals: Vec<u32> = sa.into_iter().collect();
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
        let vals: Vec<u32> = sa.into_iter().collect();
        assert!(verify_sa(&seq, &vals), "SA for large sequence failed");
    }

    #[test]
    fn test_all_same_char() {
        let seq = vec![0u8; 100];
        let sa = SuffixArray::build(&seq);
        let vals: Vec<u32> = sa.into_iter().collect();
        assert!(verify_sa(&seq, &vals), "SA for all same char failed");
    }

    #[test]
    fn test_sa_consistency() {
        for len in [1, 2, 5, 10, 20, 50, 100, 500, 1000] {
            let seq: Vec<u8> = (0..len).map(|i| (i % 5) as u8).collect();
            let sa = SuffixArray::build(&seq);
            let vals: Vec<u32> = sa.into_iter().collect();
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
        let vals_byte: Vec<u32> = sa_byte.into_iter().collect();
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
