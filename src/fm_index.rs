//! FM-Index implementation for fast substring search.

use std::fs::File;
use std::io;
use std::io::{Read, Write};
use std::path::Path;

use crate::compact::CompactOccTable;
use crate::error::BwaError;
use crate::reference::Reference;
use crate::sa::SuffixArray;

const MAGIC: &[u8; 6] = b"BWAIDX";
const VERSION: u8 = 3;

#[derive(Clone, Debug)]
pub struct BWT {
    bwt: Vec<u8>,
    len: usize,
}

impl BWT {
    /// Build BWT from SA with n+1 entries (including sentinel suffix)
    pub fn from_sa_with_sentinel(sequence: &[u8], sa: &SuffixArray) -> Self {
        let len = sa.len(); // n+1
        let n = sequence.len();
        let mut bwt = Vec::with_capacity(len);

        for &pos in &sa.sa {
            if pos == 0 {
                // Suffix starting at position 0: BWT entry is sentinel (last char is at position n-1)
                bwt.push(4);
            } else {
                // All other suffixes: BWT entry is the character before the suffix start
                // This includes the sentinel suffix at position n (pos = n)
                // For the sentinel suffix, pos-1 = n-1, which is the last character of the original sequence
                bwt.push(sequence[(pos - 1) as usize]);
            }
        }

        Self { bwt, len }
    }

    /// Build BWT from SA (legacy method for compatibility)
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
pub struct FMIndex {
    pub bwt: BWT,
    pub sa: SuffixArray,
    pub occ: CompactOccTable,
    pub f_column: [u32; 5],
    pub len: usize,
}

impl FMIndex {
    pub fn build(reference: &Reference) -> Self {
        let sequence = reference.as_slice();
        let len = sequence.len();

        // Build suffix array with n+1 entries (including sentinel)
        // We need to construct suffixes with sentinel appended, then sort
        let mut padded = Vec::with_capacity(len + 1);
        padded.extend_from_slice(&sequence);
        padded.push(4); // Sentinel character at position n

        let sa = SuffixArray::build(&padded);

        // The SA from padded sequence has len+1 entries.
        // The sentinel suffix (position n) is already in the correct sorted position.
        // We just need to use it directly as our SA.
        // The SA values still point to positions in the padded sequence,
        // so we need to remap them to original positions (position n becomes sentinel).

        let mut final_sa = Vec::with_capacity(len + 1);
        for &pos in sa.as_slice() {
            if pos as usize == len {
                // This is the sentinel suffix - keep it as sentinel
                final_sa.push(len as u32); // Use n to represent sentinel
            } else {
                final_sa.push(pos); // Keep original position
            }
        }

        let sa = SuffixArray::with_len(final_sa, len + 1);

        // Build BWT with n+1 entries (including sentinel)
        let bwt = BWT::from_sa_with_sentinel(&sequence, &sa);
        let occ = CompactOccTable::from_bwt(bwt.as_slice());

        let mut total = [0u32; 5];
        for &c in &sequence {
            total[c as usize] += 1;
        }
        // Add sentinel to total count
        total[4] = 1; // Sentinel is encoded as 4

        let mut f_column = [0u32; 5];
        // F[i] = position of first suffix starting with character i in sorted SA
        // Formula: F[i] = 1 + sum of counts of all characters before i
        // This gives 1-indexed positions: F[0]=1, F[1]=1+A, F[2]=1+A+C, etc.
        // F[4] (sentinel) = n + 1 (points past the last character)
        f_column[0] = 1; // Sentinel (or A) is first
        for i in 1..5 {
            f_column[i] = 1 + total[..i].iter().sum::<u32>();
        }
        // Override F[4] to be n + 1 (pointing past the last suffix)
        // The sentinel is implicitly at the end of the BWT, after all other characters
        f_column[4] = (len + 1) as u32;

        Self {
            bwt,
            sa,
            occ,
            f_column,
            len: len + 1, // Include sentinel in length
        }
    }

    pub fn search(&self, pattern: &[u8]) -> (usize, usize) {
        let mut left = 0;
        let mut right = self.len;

        for &c in pattern.iter().rev() {
            let occ_left = self.occ.occ(c, left);
            let occ_right = self.occ.occ(c, right);
            // F[c] is 1-indexed (first position of character c in sorted SA)
            // Subtract 1 to convert to 0-indexed before adding occ offset
            let f_c = (self.f_column[c as usize] as usize).saturating_sub(1);
            left = f_c + occ_left as usize;
            right = f_c + occ_right as usize;

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
        right.saturating_sub(left)
    }

    pub fn find_all(&self, pattern: &[u8]) -> Vec<u32> {
        let (left, right) = self.search(pattern);
        (left..right).filter_map(|i| self.sa.get(i)).collect()
    }

    pub fn save(&self, path: &Path) -> Result<(), BwaError> {
        let mut file = File::create(path)?;

        file.write_all(MAGIC)?;
        file.write_all(&[VERSION])?;
        file.write_all(&(self.len as u32).to_le_bytes())?;

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

        let mut f_column = [0u32; 5];
        for fc in &mut f_column {
            let mut bytes = [0u8; 4];
            file.read_exact(&mut bytes)?;
            *fc = u32::from_le_bytes(bytes);
        }

        let bwt = BWT::read_from(&mut file, len)?;
        let sa = SuffixArray::read_from(&mut file, len)?;
        let occ = CompactOccTable::from_bwt(bwt.as_slice());

        Ok(Self {
            bwt,
            sa,
            occ,
            f_column,
            len,
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
        assert_eq!(index.len, 9); // n+1 with sentinel

        // Check search range for single character T
        let (left, right) = index.search(&[3]); // T
        assert!(
            left < right,
            "Search for T should find matches: ({}, {})",
            left,
            right
        );

        let positions = index.find_all(&[3]);
        assert!(
            !positions.is_empty(),
            "FM-index should find T, got: {:?}",
            positions
        );

        // Now search for G
        let (left_g, right_g) = index.search(&[2]); // G
        assert!(
            left_g < right_g,
            "Search for G should find matches: ({}, {})",
            left_g,
            right_g
        );

        // Search for GT (last two chars)
        let (left_tg, right_tg) = index.search(&[2, 3]); // GT
        assert!(
            left_tg < right_tg,
            "Search for GT should find matches: ({}, {})",
            left_tg,
            right_tg
        );

        // Search for AC (first two chars)
        let (left_ac, right_ac) = index.search(&[0, 1]); // AC
        assert!(
            left_ac < right_ac,
            "Search for AC should find matches: ({}, {})",
            left_ac,
            right_ac
        );

        // Search for CG (middle two chars)
        let (left_cg, right_cg) = index.search(&[1, 2]); // CG
        assert!(
            left_cg < right_cg,
            "Search for CG should find matches: ({}, {})",
            left_cg,
            right_cg
        );
    }

    #[test]
    fn test_simple_search() {
        // Test with a simple sequence ACAC
        let ref_seq = Reference::parse_fasta(">test\nACAC").unwrap();
        let index = FMIndex::build(&ref_seq);

        // Manually trace through the search for AC
        // First search for C (pattern[1] = 1)
        let (c_left, c_right) = index.search(&[1]);
        assert!(
            c_left < c_right,
            "C search failed: ({}, {})",
            c_left,
            c_right
        );

        // Now search for AC - should use C result as starting point
        let pattern = vec![0, 1]; // AC
        let (left, right) = index.search(&pattern);
        assert!(left < right, "AC search failed: ({}, {})", left, right);

        let positions = index.find_all(&pattern);
        assert!(
            !positions.is_empty(),
            "FM-index should find AC in ACAC, got: {:?}",
            positions
        );
    }

    #[test]
    fn test_occ_directly() {
        let ref_seq = Reference::parse_fasta(">test\nACAC").unwrap();
        let index = FMIndex::build(&ref_seq);

        // BWT for ACAC is [$, C, A, A, C] = [4, 1, 0, 0, 1] with sentinel
        let bwt_slice = index.bwt.as_slice();
        assert_eq!(
            bwt_slice.len(),
            5,
            "BWT should have 5 entries (n+1 with sentinel)"
        );
        assert_eq!(bwt_slice[0], 4, "BWT[0] should be 4 (sentinel)");
        assert_eq!(bwt_slice[1], 1, "BWT[1] should be 1 (C)");
        assert_eq!(bwt_slice[2], 0, "BWT[2] should be 0 (A)");
        assert_eq!(bwt_slice[3], 0, "BWT[3] should be 0 (A)");
        assert_eq!(bwt_slice[4], 1, "BWT[4] should be 1 (C)");

        // Check occ at intermediate positions using the public API
        // A (0) appears at positions 2 and 3
        assert_eq!(
            index.occ.occ(0, 1),
            0,
            "occ(A, 1) should be 0 (BWT[0] is $)"
        );
        assert_eq!(index.occ.occ(0, 2), 0, "occ(A, 2) should be 0");
        assert_eq!(index.occ.occ(0, 3), 1, "occ(A, 3) should be 1");
        assert_eq!(index.occ.occ(0, 4), 2, "occ(A, 4) should be 2");
        assert_eq!(index.occ.occ(0, 5), 2, "occ(A, 5) should be 2");

        // C (1) appears at positions 1 and 4
        assert_eq!(
            index.occ.occ(1, 1),
            0,
            "occ(C, 1) should be 0 (BWT[0] is $)"
        );
        assert_eq!(index.occ.occ(1, 2), 1, "occ(C, 2) should be 1");
        assert_eq!(index.occ.occ(1, 4), 1, "occ(C, 4) should be 1");
        assert_eq!(index.occ.occ(1, 5), 2, "occ(C, 5) should be 2");
        assert_eq!(index.occ.occ(1, 0), 0, "occ(C, 0) should be 0");
    }

    #[test]
    fn test_reference_length() {
        let ref_seq = Reference::parse_fasta(">test\nACGT").unwrap();
        let index = FMIndex::build(&ref_seq);
        assert_eq!(index.len, 5); // n+1 with sentinel
    }

    #[test]
    fn test_save_load_roundtrip() {
        let reference = Reference::parse_fasta(">test\nACGTACGT").unwrap();
        let original = FMIndex::build(&reference);

        let pattern = [0, 1]; // AC
        let original_positions = original.find_all(&pattern);
        assert!(
            !original_positions.is_empty(),
            "Original should find AC pattern"
        );

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
        let reference =
            Reference::parse_fasta(&(">test\n".to_string() + &String::from_utf8(seq).unwrap()))
                .unwrap();
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
        std::fs::write(
            &temp_path,
            b"BWAIDX\x01\x00\x00\x00\x00\x00\x00\x00\x00\x00",
        )
        .ok();

        let result = FMIndex::load(&temp_path);
        assert!(result.is_err());

        std::fs::remove_file(temp_path).ok();
    }
}
