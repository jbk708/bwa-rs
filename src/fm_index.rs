//! FM-Index implementation for fast substring search.

use std::path::Path;
use std::io::{Read, Write};
use std::fs::File;
use std::io;

use crate::error::BwaError;
use crate::reference::Reference;

const MAGIC: &[u8; 6] = b"BWAIDX";
const VERSION: u8 = 1;

#[derive(Clone, Debug)]
pub struct SuffixArray {
    sa: Vec<u32>,
    len: usize,
}

impl SuffixArray {
    pub fn build(sequence: &[u8]) -> Self {
        let len = sequence.len();
        let mut indices: Vec<u32> = (0..len as u32).collect();
        indices.sort_by(|&a, &b| {
            sequence[a as usize..].cmp(&sequence[b as usize..])
        });
        Self { sa: indices, len }
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

    fn write_to(&self, writer: &mut impl Write) -> io::Result<()> {
        for &idx in &self.sa {
            writer.write_all(&idx.to_le_bytes())?;
        }
        Ok(())
    }

    fn read_from(reader: &mut impl Read, count: usize) -> io::Result<Self> {
        let mut sa = Vec::with_capacity(count);
        for _ in 0..count {
            let mut bytes = [0u8; 4];
            reader.read_exact(&mut bytes)?;
            sa.push(u32::from_le_bytes(bytes));
        }
        Ok(Self { sa, len: count })
    }
}

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
}

impl OccTable {
    pub fn new(bwt: &BWT, sample_rate: usize) -> Self {
        let mut counts = Vec::new();
        let mut total = [0u32; 5];
        let mut last_sample = [0u32; 5];

        for (i, &c) in bwt.as_slice().iter().enumerate() {
            total[c as usize] += 1;
            if i % sample_rate == 0 {
                counts.push(last_sample);
                last_sample = total;
            }
        }

        Self { counts, sample_rate }
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

        let start = sample_idx * self.sample_rate + if sample_idx > 0 { 1 } else { 0 };
        for _i in start..idx {
            total[c as usize] += 1;
        }

        total[c as usize]
    }

    fn write_to(&self, writer: &mut impl Write) -> io::Result<()> {
        for count in &self.counts {
            for &c in count {
                writer.write_all(&c.to_le_bytes())?;
            }
        }
        Ok(())
    }

    fn read_from(reader: &mut impl Read, sample_rate: usize, bwt_len: usize) -> io::Result<Self> {
        let count_len = bwt_len.div_ceil(sample_rate);
        let mut counts = Vec::with_capacity(count_len);
        for _ in 0..count_len {
            let mut count = [0u32; 5];
            for c in &mut count {
                let mut bytes = [0u8; 4];
                reader.read_exact(&mut bytes)?;
                *c = u32::from_le_bytes(bytes);
            }
            counts.push(count);
        }
        Ok(Self { counts, sample_rate })
    }
}

#[derive(Clone, Debug)]
pub struct FMIndex {
    bwt: BWT,
    sa: SuffixArray,
    occ: OccTable,
    f_column: [u32; 5],
    len: usize,
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
        }
    }

    pub fn search(&self, pattern: &[u8]) -> (usize, usize) {
        let mut left = 0;
        let mut right = self.len;

        for &c in pattern.iter().rev() {
            left = self.f_column[c as usize] as usize + self.occ.occ(c, left) as usize;
            right = self.f_column[c as usize] as usize + self.occ.occ(c, right) as usize;

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
    }

    #[test]
    fn test_bwt() {
        let seq = b"GACGTAC$";
        let sa = SuffixArray::build(seq);
        let bwt = BWT::from_sa(seq, &sa);
        assert_eq!(bwt.len(), 8);
    }

    #[test]
    fn test_fm_index_search() {
        let ref_seq = Reference::parse_fasta(">test\nACGTACGT").unwrap();
        let index = FMIndex::build(&ref_seq);
        assert_eq!(index.len, 8);
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

        let pattern = [0, 1, 2, 3];
        let original_positions = original.find_all(&pattern);

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
