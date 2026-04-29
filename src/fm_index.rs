//! FM-Index implementation for fast substring search.

use std::path::Path;
use std::io::{Read, Write};
use std::fs::File;
use std::io;

use crate::error::BwaError;
use crate::reference::Reference;
use crate::sa::SuffixArray;

const MAGIC: &[u8; 6] = b"BWAIDX";
const VERSION: u8 = 2;

#[derive(Clone, Debug)]
pub struct BWT {
    bwt: Vec<u8>,
    len: usize,
}

impl BWT {
    pub fn from_sa(sequence: &[u8], sa: &SuffixArray) -> Self {
        let len = sequence.len();
        let mut bwt = Vec::with_capacity(len + 1);

        for &pos in &sa.sa {
            if pos == 0 {
                bwt.push(4);
            } else {
                bwt.push(sequence[pos as usize - 1]);
            }
        }

        Self { bwt, len }
    }

    pub fn as_slice(&self) -> &[u8] {
        &self.bwt
    }

    pub fn len(&self) -> usize {
        self.len
    }

    pub fn is_empty(&self) -> bool {
        self.len == 0
    }

    pub fn get(&self, idx: usize) -> u8 {
        self.bwt[idx]
    }

    fn write_to(&self, writer: &mut impl Write) -> io::Result<()> {
        writer.write_all(&self.bwt)
    }

    fn read_from(reader: &mut impl Read, len: usize) -> io::Result<Self> {
        let mut bwt = vec![0u8; len];
        reader.read_exact(&mut bwt)?;
        Ok(Self { bwt, len })
    }
}

#[derive(Clone, Debug)]
pub struct OccTable {
    counts: Vec<[u32; 5]>,
    sample_rate: usize,
    bwt: Vec<u8>,
}

impl OccTable {
    pub fn new(bwt: &BWT, sample_rate: usize) -> Self {
        let mut counts = Vec::new();
        let mut total = [0u32; 5];
        let mut last_sample = [0u32; 5];
        let bwt_vec = bwt.as_slice().to_vec();

        for (i, &c) in bwt_vec.iter().enumerate() {
            // Don't count sentinel (value 4) in occurrence table
            if c < 4 {
                total[c as usize] += 1;
            }
            if i % sample_rate == 0 {
                counts.push(last_sample);
                last_sample = total;
            }
        }
        
        // Push the final sample if there were any characters processed
        if !bwt_vec.is_empty() {
            counts.push(total);
        }

        Self { counts, sample_rate, bwt: bwt_vec }
    }

    pub fn occ(&self, c: u8, idx: usize) -> u32 {
        if idx == 0 {
            return 0;
        }

        let sample_idx = (idx - 1) / self.sample_rate;
        let mut total = if sample_idx < self.counts.len() {
            self.counts[sample_idx]
        } else {
            [0; 5]
        };

        let start = if sample_idx == 0 { 0 } else { sample_idx * self.sample_rate + 1 };
        for i in start..idx {
            if i < self.bwt.len() && self.bwt[i] == c {
                total[c as usize] += 1;
            }
        }

        total[c as usize]
    }

    fn write_to(&self, writer: &mut impl Write) -> io::Result<()> {
        writer.write_all(&self.bwt)
    }

    fn read_from(reader: &mut impl Read, sample_rate: usize, bwt_len: usize) -> io::Result<Self> {
        let mut bwt = vec![0u8; bwt_len];
        reader.read_exact(&mut bwt)?;
        
        // Recreate the occurrence table from the BWT
        let mut counts = Vec::new();
        let mut total = [0u32; 5];
        let mut last_sample = [0u32; 5];

        for (i, &c) in bwt.iter().enumerate() {
            if c < 4 {
                total[c as usize] += 1;
            }
            if i % sample_rate == 0 {
                counts.push(last_sample);
                last_sample = total;
            }
        }
        
        if !bwt.is_empty() {
            counts.push(total);
        }

        Ok(Self { counts, sample_rate, bwt })
    }
}

#[derive(Clone, Debug)]
pub struct FMIndex {
    bwt: BWT,
    sa: SuffixArray,
    occ: OccTable,
    f_column: [u32; 5],
    len: usize,
    #[allow(dead_code)]
    reference: Vec<u8>,
}

impl FMIndex {
    pub fn build(reference: &Reference) -> Self {
        let sequence = reference.as_slice();
        let len = sequence.len();

        let sa = SuffixArray::build(&sequence);
        let bwt = BWT::from_sa(&sequence, &sa);
        let occ = OccTable::new(&bwt, 32);

        let mut total = [0u32; 5];
        for &c in &sequence {
            total[c as usize] += 1;
        }

        let mut f_column = [0u32; 5];
        f_column[0] = 1;
        for i in 1..5 {
            f_column[i] = f_column[i - 1] + total[i - 1];
        }

        Self {
            bwt,
            sa,
            occ,
            f_column,
            len,
            reference: sequence.to_vec(),
        }
    }

    #[allow(dead_code)]
    pub fn reference(&self) -> &[u8] {
        &self.reference
    }

    pub fn search(&self, pattern: &[u8]) -> (usize, usize) {
        let mut left = 0;
        let mut right = self.len;

        for &c in pattern.iter().rev() {
            let occ_left = self.occ.occ(c, left);
            let occ_right = self.occ.occ(c, right);
            left = self.f_column[c as usize] as usize + occ_left as usize;
            right = self.f_column[c as usize] as usize + occ_right as usize;

            if left >= right {
                break;
            }
        }

        (left, right)
    }

    pub fn get_position(&self, sa_idx: usize) -> Option<u32> {
        self.sa.get(sa_idx)
    }

    pub fn count(&self, pattern: &[u8]) -> usize {
        let (left, right) = self.search(pattern);
        right - left
    }

    pub fn find_all(&self, pattern: &[u8]) -> Vec<u32> {
        let (left, right) = self.search(pattern);
        (left..right)
            .filter_map(|i| self.sa.get(i))
            .collect()
    }

    pub fn save(&self, path: &Path) -> Result<(), BwaError> {
        let mut file = File::create(path)?;

        file.write_all(MAGIC)?;
        file.write_all(&[VERSION])?;
        file.write_all(&(self.len as u32).to_le_bytes())?;
        file.write_all(&(self.occ.sample_rate as u32).to_le_bytes())?;

        for &fc in &self.f_column {
            file.write_all(&fc.to_le_bytes())?;
        }

        self.bwt.write_to(&mut file)?;
        self.sa.write_to(&mut file)?;
        self.occ.write_to(&mut file)?;

        Ok(())
    }

    pub fn load(path: &Path) -> Result<Self, BwaError> {
        let mut file = File::open(path)?;

        let mut magic = [0u8; 6];
        file.read_exact(&mut magic)?;
        if &magic != MAGIC {
            return Err(BwaError::Index("Invalid index file: bad magic".into()));
        }

        let mut version = [0u8; 1];
        file.read_exact(&mut version)?;
        if version[0] != VERSION {
            return Err(BwaError::Index(format!(
                "Unsupported index version: expected {}, got {}",
                VERSION, version[0]
            )));
        }

        let mut len_bytes = [0u8; 4];
        file.read_exact(&mut len_bytes)?;
        let len = u32::from_le_bytes(len_bytes) as usize;

        let mut rate_bytes = [0u8; 4];
        file.read_exact(&mut rate_bytes)?;
        let sample_rate = u32::from_le_bytes(rate_bytes) as usize;

        let mut f_column = [0u32; 5];
        for fc in &mut f_column {
            let mut bytes = [0u8; 4];
            file.read_exact(&mut bytes)?;
            *fc = u32::from_le_bytes(bytes);
        }

        let bwt = BWT::read_from(&mut file, len)?;
        let sa = SuffixArray::read_from(&mut file, len)?;
        let occ = OccTable::read_from(&mut file, sample_rate, len)?;

        Ok(Self {
            bwt,
            sa,
            occ,
            f_column,
            len,
            reference: Vec::new(),
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;


    #[test]
    fn test_suffix_array() {
        let seq = b"AACGAACGG";
        let sa = SuffixArray::build(seq);
        assert_eq!(sa.len(), 9);
        
        // Check SA contains expected positions
        let positions: Vec<u32> = (0..9).filter_map(|i| sa.get(i)).collect();
        assert_eq!(positions.len(), 9);
    }
    
    #[test]
    fn test_suffix_array_order() {
        let seq = vec![0u8, 1, 2, 3, 0, 1, 2, 3]; // ACGTACGT
        let sa = SuffixArray::build(&seq);
        
        // Check that SA is sorted correctly
        let sa_vals: Vec<u32> = (0..8).filter_map(|i| sa.get(i)).collect();
        
        // The suffixes should be sorted:
        // Position 4: ACGT (shorter, comes first)
        // Position 0: ACGTACGT
        // Position 5: CGT
        // Position 1: CGTACGT
        // Position 6: GT
        // Position 2: GTACGT
        // Position 7: T
        // Position 3: TACGT
        assert_eq!(sa_vals, vec![4, 0, 5, 1, 6, 2, 7, 3]);
    }

    #[test]
    fn test_bwt() {
        let seq = vec![0u8, 1, 2, 3, 0, 1, 2, 3]; // ACGTACGT
        let sa = SuffixArray::build(&seq);
        let bwt = BWT::from_sa(&seq, &sa);
        assert_eq!(bwt.len(), 8);
        
        // SA = [4, 0, 5, 1, 6, 2, 7, 3]
        // For pos=4: sequence[3] = 3 (T), so BWT[0] = 3
        // For pos=0: push 4 (sentinel), so BWT[1] = 4
        // For pos=5: sequence[4] = 0 (A), so BWT[2] = 0
        // etc.
        // BWT should be [3, 4, 0, 1, 2, 3, 0, 1]
        assert_eq!(bwt.get(0), 3, "BWT[0] should be 3 (T from sequence[3])");
        assert_eq!(bwt.get(1), 4, "BWT[1] should be 4 (sentinel from SA[1]=0)");
    }

    #[test]
    fn test_fm_index_search() {
        let ref_seq = Reference::parse_fasta(">test\nACGTACGT").unwrap();
        let index = FMIndex::build(&ref_seq);
        assert_eq!(index.len, 8);
        
        // Check search range for single character T
        let (left, right) = index.search(&[3]); // T
        assert!(left < right, "Search for T should find matches: ({}, {})", left, right);
        
        let positions = index.find_all(&[3]);
        assert!(!positions.is_empty(), "FM-index should find T, got: {:?}", positions);
        
        // Now search for G
        let (left_g, right_g) = index.search(&[2]); // G
        assert!(left_g < right_g, "Search for G should find matches: ({}, {})", left_g, right_g);
        
        // Search for GT (last two chars)
        let (left_tg, right_tg) = index.search(&[2, 3]); // GT
        assert!(left_tg < right_tg, "Search for GT should find matches: ({}, {})", left_tg, right_tg);
        
        // Search for AC (first two chars)
        let (left_ac, right_ac) = index.search(&[0, 1]); // AC
        assert!(left_ac < right_ac, "Search for AC should find matches: ({}, {})", left_ac, right_ac);
        
        // Search for CG (middle two chars)
        let (left_cg, right_cg) = index.search(&[1, 2]); // CG
        assert!(left_cg < right_cg, "Search for CG should find matches: ({}, {})", left_cg, right_cg);
    }
    
    #[test]
    fn test_simple_search() {
        // Test with a simple sequence ACAC
        let ref_seq = Reference::parse_fasta(">test\nACAC").unwrap();
        let index = FMIndex::build(&ref_seq);
        
        // Manually trace through the search for AC
        // First search for C (pattern[1] = 1)
        let (c_left, c_right) = index.search(&[1]);
        assert!(c_left < c_right, "C search failed: ({}, {})", c_left, c_right);
        
        // Now search for AC - should use C result as starting point
        let pattern = vec![0, 1]; // AC
        let (left, right) = index.search(&pattern);
        assert!(left < right, "AC search failed: ({}, {})", left, right);
        
        let positions = index.find_all(&pattern);
        assert!(!positions.is_empty(), "FM-index should find AC in ACAC, got: {:?}", positions);
    }
    
    #[test]
    fn test_occ_directly() {
        let ref_seq = Reference::parse_fasta(">test\nACAC").unwrap();
        let index = FMIndex::build(&ref_seq);
        
        // Check BWT directly - it should be [1, 4, 0, 0]
        // For ACAC: SA=[2, 0, 3, 1], BWT[0]=sequence[1]=C=1, BWT[1]=4, BWT[2]=sequence[2]=A=0, BWT[3]=sequence[0]=A=0
        let bwt_slice = index.bwt.as_slice();
        assert_eq!(bwt_slice.len(), 4, "BWT should have 4 entries");
        assert_eq!(bwt_slice[0], 1, "BWT[0] should be 1 (C)");
        assert_eq!(bwt_slice[1], 4, "BWT[1] should be 4 (sentinel)");
        assert_eq!(bwt_slice[2], 0, "BWT[2] should be 0 (A)");
        assert_eq!(bwt_slice[3], 0, "BWT[3] should be 0 (A)");
        
        // Debug: check the counts array
        assert_eq!(index.occ.counts.len(), 2, "Should have 2 samples: initial and final");
        assert_eq!(index.occ.counts[0], [0, 0, 0, 0, 0], "First sample should be all zeros");
        assert_eq!(index.occ.counts[1], [2, 1, 0, 0, 0], "Second sample should be final counts");
        
        // Check occ at intermediate positions
        assert_eq!(index.occ.occ(0, 1), 0, "occ(A, 1) should be 0");
        assert_eq!(index.occ.occ(0, 2), 0, "occ(A, 2) should be 0");
        assert_eq!(index.occ.occ(0, 3), 1, "occ(A, 3) should be 1");
        assert_eq!(index.occ.occ(0, 4), 2, "occ(A, 4) should be 2");
        
        // Check occ for C
        assert_eq!(index.occ.occ(1, 1), 1, "occ(C, 1) should be 1");
        assert_eq!(index.occ.occ(1, 4), 1, "occ(C, 4) should be 1");
        assert_eq!(index.occ.occ(1, 0), 0, "occ(C, 0) should be 0");
    }

    #[test]
    fn test_reference_length() {
        let ref_seq = Reference::parse_fasta(">test\nACGT").unwrap();
        let index = FMIndex::build(&ref_seq);
        assert_eq!(index.len, 4);
    }

    #[test]
    fn test_save_load_roundtrip() {
        let reference = Reference::parse_fasta(">test\nACGTACGT").unwrap();
        let original = FMIndex::build(&reference);

        let pattern = [0, 1]; // AC
        let original_positions = original.find_all(&pattern);
        assert!(!original_positions.is_empty(), "Original should find AC pattern");

        let temp_path = std::env::temp_dir().join("test_index.idx");
        original.save(&temp_path).unwrap();

        let loaded = FMIndex::load(&temp_path).unwrap();
        let loaded_positions = loaded.find_all(&pattern);

        assert_eq!(original_positions, loaded_positions);
        assert_eq!(original.len, loaded.len);

        std::fs::remove_file(temp_path).ok();
    }

    #[test]
    fn test_save_load_large() {
        let seq = b"ACGT".repeat(1000);
        let reference = Reference::parse_fasta(&(">test\n".to_string() + &String::from_utf8(seq).unwrap())).unwrap();
        let original = FMIndex::build(&reference);

        let pattern = [0, 1, 2];
        let original_count = original.count(&pattern);

        let temp_path = std::env::temp_dir().join("test_large_index.idx");
        original.save(&temp_path).unwrap();

        let loaded = FMIndex::load(&temp_path).unwrap();

        assert_eq!(original_count, loaded.count(&pattern));
        assert_eq!(original.len, loaded.len);

        std::fs::remove_file(temp_path).ok();
    }

    #[test]
    fn test_load_invalid_magic() {
        let temp_path = std::env::temp_dir().join("test_invalid.idx");
        std::fs::write(&temp_path, b"INVALID").ok();

        let result = FMIndex::load(&temp_path);
        assert!(result.is_err());

        std::fs::remove_file(temp_path).ok();
    }

    #[test]
    fn test_load_invalid_version() {
        let temp_path = std::env::temp_dir().join("test_bad_version.idx");
        let mut file = std::fs::File::create(&temp_path).unwrap();
        file.write_all(b"BWAIDX").unwrap();
        file.write_all(&[99u8]).unwrap();

        let result = FMIndex::load(&temp_path);
        assert!(result.is_err());

        std::fs::remove_file(temp_path).ok();
    }

    #[test]
    fn test_load_corrupted() {
        let temp_path = std::env::temp_dir().join("test_corrupt.idx");
        std::fs::write(&temp_path, b"BWAIDX\x01\x00\x00\x00\x00\x00\x00\x00\x00\x00").ok();

        let result = FMIndex::load(&temp_path);
        assert!(result.is_err());

        std::fs::remove_file(temp_path).ok();
    }
}
