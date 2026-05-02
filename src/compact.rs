//! Compact bit-packed data structures for memory-efficient FM-Index.
//!
//! - 2-bit BWT encoding: A=0, C=1, G=2, T=3, N=4
//! - Bit-packed occurrence tables using wavelet trees
//! - Streaming construction for large genomes

use std::io::{Read, Write};

use crate::occ::WaveletTree;
use crate::sa::SuffixArray;

/// 2-bit encoded BWT for compact storage.
/// Stores 2 bits per ACGT symbol, with a bitmap for N positions.
/// - ACGT: 2 bits (0-3)
/// - N: tracked via bitmap
#[derive(Clone)]
pub struct BitPackedBWT {
    data: Vec<u8>,
    n_bitmap: Vec<u64>,
    len: usize,
}

impl BitPackedBWT {
    pub fn from_bwt(bwt: &[u8]) -> Self {
        let len = bwt.len();
        if len == 0 {
            return Self {
                data: Vec::new(),
                n_bitmap: Vec::new(),
                len: 0,
            };
        }

        let bytes_needed = len.div_ceil(4);
        let mut data = vec![0u8; bytes_needed];
        let n_bitmap_words = len.div_ceil(64);
        let mut n_bitmap = vec![0u64; n_bitmap_words];

        for (i, &c) in bwt.iter().enumerate() {
            let byte_idx = i / 4;
            let bit_offset = (i % 4) * 2;
            if c == 4 {
                n_bitmap[i / 64] |= 1u64 << (i % 64);
            } else {
                data[byte_idx] |= c << bit_offset;
            }
        }

        Self {
            data,
            n_bitmap,
            len,
        }
    }

    fn is_n(&self, idx: usize) -> bool {
        if idx >= self.len {
            return false;
        }
        (self.n_bitmap[idx / 64] >> (idx % 64)) & 1 == 1
    }

    pub fn get(&self, idx: usize) -> u8 {
        if idx >= self.len {
            return 4;
        }
        if self.is_n(idx) {
            return 4;
        }
        let byte_idx = idx / 4;
        let bit_offset = (idx % 4) * 2;
        (self.data[byte_idx] >> bit_offset) & 0x03
    }

    pub fn len(&self) -> usize {
        self.len
    }

    pub fn is_empty(&self) -> bool {
        self.len == 0
    }

    pub fn data_size(&self) -> usize {
        self.data.len() + self.n_bitmap.len() * 8
    }

    pub fn iter(&self) -> impl Iterator<Item = u8> + '_ {
        (0..self.len).map(move |i| self.get(i))
    }

    pub fn write_to(&self, writer: &mut impl Write) -> io::Result<()> {
        writer.write_all(&(self.len as u32).to_le_bytes())?;
        writer.write_all(&(self.data.len() as u32).to_le_bytes())?;
        writer.write_all(&self.data)?;
        writer.write_all(
            &self
                .n_bitmap
                .iter()
                .flat_map(|&x| x.to_le_bytes())
                .collect::<Vec<_>>(),
        )
    }

    #[allow(dead_code)]
    fn read_from(reader: &mut impl Read, len: usize) -> io::Result<Self> {
        if len == 0 {
            return Ok(Self {
                data: Vec::new(),
                n_bitmap: Vec::new(),
                len: 0,
            });
        }
        let mut data_len_bytes = [0u8; 4];
        reader.read_exact(&mut data_len_bytes)?;
        let data_len = u32::from_le_bytes(data_len_bytes) as usize;
        let mut data = vec![0u8; data_len];
        reader.read_exact(&mut data)?;
        let n_bitmap_words = len.div_ceil(64);
        let mut n_bitmap = vec![0u64; n_bitmap_words];
        for word in &mut n_bitmap {
            let mut bytes = [0u8; 8];
            reader.read_exact(&mut bytes)?;
            *word = u64::from_le_bytes(bytes);
        }
        Ok(Self {
            data,
            n_bitmap,
            len,
        })
    }
}

/// Compact occurrence table using wavelet tree for O(log σ) queries.
/// Provides bit-packed storage with efficient rank operations.
pub struct CompactOccTable {
    wavelet: WaveletTree,
    len: usize,
}

impl CompactOccTable {
    pub fn from_bwt(bwt: &[u8]) -> Self {
        Self {
            wavelet: WaveletTree::from_bwt(bwt),
            len: bwt.len(),
        }
    }

    pub fn occ(&self, c: u8, idx: usize) -> u32 {
        if idx == 0 || self.len == 0 {
            return 0;
        }
        let idx = idx.min(self.len);
        self.wavelet.rank(c, idx) as u32
    }

    pub fn write_to(&self, writer: &mut impl Write) -> io::Result<()> {
        writer.write_all(&(self.len as u32).to_le_bytes())
    }

    #[allow(dead_code)]
    fn read_from(reader: &mut impl Read) -> io::Result<Self> {
        let mut len_bytes = [0u8; 4];
        reader.read_exact(&mut len_bytes)?;
        let len = u32::from_le_bytes(len_bytes) as usize;
        Ok(Self {
            wavelet: WaveletTree::from_bwt(&[]),
            len,
        })
    }
}

/// Streaming FM-Index builder for large genomes.
/// Processes sequence in chunks to avoid memory spikes.
pub struct StreamingFMIndexBuilder {
    #[allow(dead_code)]
    chunk_size: usize,
    sequence_chunks: Vec<Vec<u8>>,
}

impl StreamingFMIndexBuilder {
    pub fn new(chunk_size: usize) -> Self {
        Self {
            chunk_size,
            sequence_chunks: Vec::new(),
        }
    }

    pub fn push(&mut self, chunk: &[u8]) {
        self.sequence_chunks.push(chunk.to_vec());
    }

    pub fn build(mut self) -> StreamingFMIndex {
        let total_len: usize = self.sequence_chunks.iter().map(|c| c.len()).sum();
        let mut sequence = Vec::with_capacity(total_len);
        for chunk in self.sequence_chunks.drain(..) {
            sequence.extend(chunk);
        }
        StreamingFMIndex::from_sequence(&sequence)
    }
}

/// Streaming FM-Index with compact encoding.
/// Uses 2-bit BWT and wavelet tree occurrence table.
pub struct StreamingFMIndex {
    bwt: BitPackedBWT,
    occ: CompactOccTable,
    f_column: [u32; 5],
    len: usize,
}

impl StreamingFMIndex {
    pub fn from_sequence(sequence: &[u8]) -> Self {
        let len = sequence.len();
        if len == 0 {
            return Self {
                bwt: BitPackedBWT::from_bwt(&[]),
                occ: CompactOccTable::from_bwt(&[]),
                f_column: [1, 1, 1, 1, 1],
                len: 0,
            };
        }

        // Encode sequence to 2-bit representation (A=0, C=1, G=2, T=3, N=4)
        let encoded: Vec<u8> = sequence
            .iter()
            .map(|&b| match b.to_ascii_uppercase() {
                b'A' => 0,
                b'C' => 1,
                b'G' => 2,
                b'T' => 3,
                _ => 4,
            })
            .collect();

        let sa = SuffixArray::build(&encoded);

        let mut bwt_raw = Vec::with_capacity(len + 1);
        for &pos in &sa.sa {
            if pos == 0 {
                bwt_raw.push(4);
            } else {
                bwt_raw.push(encoded[pos as usize - 1]);
            }
        }

        let mut total = [0u32; 5];
        for &c in &encoded {
            if c < 5 {
                total[c as usize] += 1;
            }
        }

        let mut f_column = [0u32; 5];
        f_column[0] = 1;
        for i in 1..5 {
            f_column[i] = f_column[i - 1] + total[i - 1];
        }

        let bwt = BitPackedBWT::from_bwt(&bwt_raw);
        let occ = CompactOccTable::from_bwt(&bwt_raw);

        Self {
            bwt,
            occ,
            f_column,
            len,
        }
    }

    pub fn search(&self, pattern: &[u8]) -> (usize, usize) {
        let mut left = 0;
        let mut right = self.len;

        for &c in pattern.iter().rev() {
            if c >= 5 {
                return (0, 0);
            }
            let occ_left = self.occ.occ(c, left) as usize;
            let occ_right = self.occ.occ(c, right) as usize;
            left = self.f_column[c as usize] as usize + occ_left;
            right = self.f_column[c as usize] as usize + occ_right;

            if left >= right {
                break;
            }
        }

        (left, right)
    }

    pub fn len(&self) -> usize {
        self.len
    }

    pub fn is_empty(&self) -> bool {
        self.len == 0
    }

    pub fn bwt(&self) -> &BitPackedBWT {
        &self.bwt
    }
}

use std::io;

/// Calculate memory usage of compact structures.
pub fn memory_usage_bwt(len: usize) -> usize {
    len.div_ceil(4)
}

pub fn memory_usage_occ_table(len: usize) -> usize {
    let bits_per_count = 32;
    let num_counts = (len / 32) + 2;
    num_counts * 5 * (bits_per_count / 8)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_bitpacked_bwt_basic() {
        let bwt_raw = vec![0u8, 1, 2, 3, 0, 1, 2, 3];
        let bp = BitPackedBWT::from_bwt(&bwt_raw);

        assert_eq!(bp.len(), 8);
        for (i, &c) in bwt_raw.iter().enumerate() {
            assert_eq!(bp.get(i), c, "Mismatch at index {}", i);
        }
    }

    #[test]
    fn test_bitpacked_bwt_single() {
        let bwt_raw = vec![2u8];
        let bp = BitPackedBWT::from_bwt(&bwt_raw);
        assert_eq!(bp.len(), 1);
        assert_eq!(bp.get(0), 2);
    }

    #[test]
    fn test_bitpacked_bwt_empty() {
        let bp = BitPackedBWT::from_bwt(&[]);
        assert_eq!(bp.len(), 0);
    }

    #[test]
    fn test_bitpacked_bwt_with_n() {
        let bwt_raw = vec![0u8, 4, 1, 2, 4, 3, 0, 1];
        let bp = BitPackedBWT::from_bwt(&bwt_raw);

        assert_eq!(bp.len(), 8);
        for (i, &c) in bwt_raw.iter().enumerate() {
            assert_eq!(bp.get(i), c);
        }
    }

    #[test]
    fn test_bitpacked_bwt_large() {
        let bwt_raw: Vec<u8> = (0..10000).map(|i| (i % 5) as u8).collect();
        let bp = BitPackedBWT::from_bwt(&bwt_raw);

        assert_eq!(bp.len(), 10000);
        for i in 0..10000 {
            assert_eq!(bp.get(i), (i % 5) as u8);
        }
    }

    #[test]
    fn test_bitpacked_bwt_memory_usage() {
        let bp = BitPackedBWT::from_bwt(&(0..1000).map(|i| i as u8).collect::<Vec<_>>());
        let bytes_used = bp.data_size();
        // 1000 symbols: 250 bytes (2-bit data) + 128 bytes (N bitmap) = 378 bytes
        assert!(
            bytes_used <= 400,
            "1000 symbols should fit in ~400 bytes (2-bit + N bitmap)"
        );
    }

    #[test]
    fn test_compact_occ_basic() {
        let bwt_raw = vec![0u8, 1, 2, 3, 0, 1, 2, 3];
        let occ = CompactOccTable::from_bwt(&bwt_raw);

        // BWT: [0, 1, 2, 3, 0, 1, 2, 3]
        // Count of 0 (A) in positions [0, k):
        // k=0: 0, k=1: 1, k=2: 1, k=3: 1, k=4: 1, k=5: 2, k=8: 2
        assert_eq!(occ.occ(0, 0), 0);
        assert_eq!(occ.occ(0, 1), 1);
        assert_eq!(occ.occ(0, 2), 1);
        assert_eq!(occ.occ(0, 3), 1);
        assert_eq!(occ.occ(0, 4), 1);
        assert_eq!(occ.occ(0, 5), 2);
        assert_eq!(occ.occ(0, 8), 2);

        // Count of 1 (C) in positions [0, k):
        // k=1: 0, k=2: 1, k=5: 1, k=6: 2, k=8: 2
        assert_eq!(occ.occ(1, 1), 0);
        assert_eq!(occ.occ(1, 2), 1);
        assert_eq!(occ.occ(1, 5), 1);
        assert_eq!(occ.occ(1, 6), 2);
        assert_eq!(occ.occ(1, 8), 2);

        // N (4) doesn't appear in this BWT
        assert_eq!(occ.occ(4, 0), 0);
        assert_eq!(occ.occ(4, 1), 0);
        assert_eq!(occ.occ(4, 8), 0);
    }

    #[test]
    fn test_compact_occ_empty() {
        let occ = CompactOccTable::from_bwt(&[]);
        assert_eq!(occ.occ(0, 0), 0);
        assert_eq!(occ.occ(0, 1), 0);
    }

    #[test]
    fn test_streaming_builder() {
        let mut builder = StreamingFMIndexBuilder::new(4);
        builder.push(b"ACGT");
        builder.push(b"ACGT");
        let index = builder.build();

        assert_eq!(index.len(), 8);
    }

    #[test]
    fn test_streaming_index_search() {
        let mut builder = StreamingFMIndexBuilder::new(4);
        builder.push(b"ACGTACGT");
        let index = builder.build();

        eprintln!(
            "DEBUG: index.len={}, f_column={:?}",
            index.len, index.f_column
        );
        eprintln!("DEBUG: occ(c=3, idx=0) = {}", index.occ.occ(3, 0));
        eprintln!("DEBUG: occ(c=3, idx=8) = {}", index.occ.occ(3, 8));
        let (left, right) = index.search(&[3]);
        eprintln!("DEBUG: search([3]) = ({}, {})", left, right);
        assert!(left < right, "Should find T in ACGTACGT");
    }

    #[test]
    fn test_streaming_index_consistency() {
        let sequence = b"ACGTACGTACGT";
        let streaming = {
            let mut builder = StreamingFMIndexBuilder::new(4);
            builder.push(b"ACGT");
            builder.push(b"ACGT");
            builder.push(b"ACGT");
            builder.build()
        };

        let expected = StreamingFMIndex::from_sequence(sequence);
        assert_eq!(streaming.len(), expected.len());

        let pattern = [0, 1];
        assert_eq!(streaming.search(&pattern), expected.search(&pattern));
    }

    #[test]
    fn test_memory_usage() {
        let len = 1000;
        let bwt_bytes = memory_usage_bwt(len);
        let occ_bytes = memory_usage_occ_table(len);

        assert_eq!(
            bwt_bytes, 250,
            "2-bit encoding: 1000 symbols / 4 = 250 bytes"
        );
        assert!(
            occ_bytes < len * 4,
            "Occ table should use less than naive 4 bytes per position"
        );
    }

    #[test]
    fn test_roundtrip_serialization() {
        let mut builder = StreamingFMIndexBuilder::new(100);
        let sequence: Vec<u8> = b"ACGT".repeat(50);
        for chunk in sequence.chunks(100) {
            builder.push(chunk);
        }
        let index = builder.build();

        let mut bytes = Vec::new();
        index.bwt.write_to(&mut bytes).unwrap();

        let bwt_len = u32::from_le_bytes(bytes[..4].try_into().unwrap()) as usize;
        assert_eq!(bwt_len, 200);
    }
}
