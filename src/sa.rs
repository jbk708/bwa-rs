//! Suffix Array construction module using SA-IS algorithm.
//!
//! This module provides O(n) suffix array construction via the sa-is crate.

use std::io::{Read, Write};

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
        
        Self { sa, len: sequence.len() }
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
    sequence.chunks(chunk_size).map(|chunk| {
        SuffixArray::build(chunk).sa
    })
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
}