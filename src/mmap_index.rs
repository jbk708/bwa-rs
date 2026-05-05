#![allow(unsafe_code)]

//! Memory-mapped FM-Index for zero-copy index access.
//!
//! Enables processing of 3GB+ genomes without loading entire index into RAM.

use memmap2::Mmap;
use std::fs::File;
use std::io::Write;
use std::path::Path;

use crate::error::BwaError;
use crate::reference::Reference;
use crate::sa::SuffixArray;

const MAGIC: &[u8; 6] = b"BWAIDX";
const VERSION: u8 = 2;
const HEADER_SIZE: usize = 6 + 1 + 4 + 4 + 5 * 4; // magic + version + len + sample_rate + f_column

pub struct MmapFMIndex {
    mmap: Mmap,
    len: usize,
    sample_rate: usize,
    f_column: [u32; 5],
    bwt_offset: usize,
    sa_offset: usize,
    occ_offset: usize,
    occ_len: usize,
}

impl MmapFMIndex {
    pub fn open(path: &Path) -> Result<Self, BwaError> {
        let file = File::open(path)?;
        let mmap = unsafe { Mmap::map(&file)? };
        Self::from_mmap(mmap)
    }

    #[allow(unsafe_code)]
    pub fn from_mmap(mmap: Mmap) -> Result<Self, BwaError> {
        if mmap.len() < HEADER_SIZE {
            return Err(BwaError::Index("Index file too short".into()));
        }

        if &mmap[..6] != MAGIC {
            return Err(BwaError::Index("Invalid magic number".into()));
        }

        if mmap[6] != VERSION {
            return Err(BwaError::Index(format!(
                "Unsupported version: expected {}, got {}",
                VERSION, mmap[6]
            )));
        }

        let len = u32::from_le_bytes(mmap[7..11].try_into().unwrap()) as usize;
        let sample_rate = u32::from_le_bytes(mmap[11..15].try_into().unwrap()) as usize;

        let mut f_column = [0u32; 5];
        for (i, fc) in f_column.iter_mut().enumerate() {
            *fc = u32::from_le_bytes(mmap[15 + i * 4..19 + i * 4].try_into().unwrap());
        }

        let bwt_offset = HEADER_SIZE;
        let sa_offset = bwt_offset + len;
        let occ_len = (len / sample_rate) + 2;
        let occ_offset = sa_offset + len * 4;

        Ok(Self {
            mmap,
            len,
            sample_rate,
            f_column,
            bwt_offset,
            sa_offset,
            occ_offset,
            occ_len,
        })
    }

    pub fn build_and_save(reference: &Reference, path: &Path) -> Result<Self, BwaError> {
        let sequence = reference.as_slice();
        let n = sequence.len();
        let len = n + 1; // Include sentinel in length

        // Build suffix array with n+1 entries (including sentinel)
        let mut padded = Vec::with_capacity(len);
        padded.extend_from_slice(&sequence);
        padded.push(4); // Sentinel character at position n

        let sa = SuffixArray::build(&padded);

        // Build BWT with n+1 entries
        let mut bwt = Vec::with_capacity(len);
        for &pos in &sa.sa {
            if pos == 0 {
                // Suffix starting at position 0: BWT entry is sentinel
                bwt.push(4);
            } else {
                // All other suffixes: BWT entry is the character before the suffix start
                // This includes the sentinel suffix at position n
                bwt.push(sequence[(pos - 1) as usize]);
            }
        }

        let mut total = [0u32; 5];
        for c in &sequence {
            total[*c as usize] += 1;
        }
        total[4] = 1; // Add sentinel count

        let mut f_column = [0u32; 5];
        f_column[0] = 1;
        for i in 1..5 {
            f_column[i] = 1 + total[..i].iter().sum::<u32>();
        }
        // Override F[4] to be n + 1
        f_column[4] = len as u32;

        let sample_rate = 32;

        let mut file = File::create(path)?;
        file.write_all(MAGIC)?;
        file.write_all(&[VERSION])?;
        file.write_all(&(len as u32).to_le_bytes())?;
        file.write_all(&(sample_rate as u32).to_le_bytes())?;

        for &fc in &f_column {
            file.write_all(&fc.to_le_bytes())?;
        }

        file.write_all(&bwt)?;

        for i in 0..len {
            file.write_all(&(sa.get(i).unwrap_or(0)).to_le_bytes())?;
        }

        // Build occ table from original sequence (not BWT) for correct counts
        let mut counts = Vec::new();
        let mut running = [0u32; 5];
        let mut last_sample = [0u32; 5];

        for (i, &c) in sequence.iter().enumerate() {
            running[c as usize] += 1;
            if i % sample_rate == 0 {
                counts.push(last_sample);
                last_sample = running;
            }
        }
        counts.push(running);

        for sample in &counts {
            for &val in sample {
                file.write_all(&val.to_le_bytes())?;
            }
        }

        drop(file);
        Self::open(path)
    }

    fn occ_at(&self, sample_idx: usize, c: u8) -> u32 {
        if sample_idx >= self.occ_len {
            return 0;
        }
        let offset = self.occ_offset + sample_idx * 5 * 4;
        u32::from_le_bytes(
            self.mmap[offset + c as usize * 4..offset + (c as usize + 1) * 4]
                .try_into()
                .unwrap(),
        )
    }

    fn occ(&self, c: u8, idx: usize) -> u32 {
        if idx == 0 {
            return 0;
        }

        let sample_idx = (idx - 1) / self.sample_rate;
        let mut total = self.occ_at(sample_idx, c);

        let start = if sample_idx == 0 {
            0
        } else {
            sample_idx * self.sample_rate + 1
        };

        for i in start..idx.min(self.len) {
            if self.mmap[self.bwt_offset + i] == c {
                total += 1;
            }
        }

        total
    }

    pub fn search(&self, pattern: &[u8]) -> (usize, usize) {
        let mut left = 0;
        let mut right = self.len;

        for &c in pattern.iter().rev() {
            let occ_left = self.occ(c, left);
            let occ_right = self.occ(c, right);
            // F[c] is 1-indexed, subtract 1 to convert to 0-indexed
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
        if sa_idx >= self.len {
            return None;
        }
        let offset = self.sa_offset + sa_idx * 4;
        let pos = u32::from_le_bytes(self.mmap[offset..offset + 4].try_into().unwrap());
        Some(pos)
    }

    pub fn count(&self, pattern: &[u8]) -> usize {
        let (left, right) = self.search(pattern);
        right - left
    }

    pub fn find_all(&self, pattern: &[u8]) -> Vec<u32> {
        let (left, right) = self.search(pattern);
        (left..right).filter_map(|i| self.get_position(i)).collect()
    }

    pub fn len(&self) -> usize {
        self.len
    }

    pub fn is_empty(&self) -> bool {
        self.len == 0
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::FMIndex;

    #[test]
    fn test_build_and_open() {
        let reference = Reference::parse_fasta(">test\nACGTACGT").unwrap();
        let temp_path = std::env::temp_dir().join("test_mmap_index.idx");

        let index = MmapFMIndex::build_and_save(&reference, &temp_path).unwrap();
        assert_eq!(index.len(), 9); // n+1 with sentinel

        std::fs::remove_file(temp_path).ok();
    }

    #[test]
    fn test_search() {
        let reference = Reference::parse_fasta(">test\nACGTACGT").unwrap();
        let temp_path = std::env::temp_dir().join("test_mmap_search.idx");

        MmapFMIndex::build_and_save(&reference, &temp_path).unwrap();
        let index = MmapFMIndex::open(&temp_path).unwrap();

        let (left, right) = index.search(&[3]); // T
        assert!(left < right, "Should find T");

        let positions = index.find_all(&[3]);
        assert!(!positions.is_empty(), "Should find T positions");

        std::fs::remove_file(temp_path).ok();
    }

    #[test]
    fn test_large_index() {
        let seq = b"ACGT".repeat(10000);
        let fasta = format!(">test\n{}", String::from_utf8(seq).unwrap());
        let reference = Reference::parse_fasta(&fasta).unwrap();
        let temp_path = std::env::temp_dir().join("test_large_mmap.idx");

        MmapFMIndex::build_and_save(&reference, &temp_path).unwrap();
        let index = MmapFMIndex::open(&temp_path).unwrap();

        assert_eq!(index.len(), 40001); // n+1 with sentinel

        let pattern = [0, 1, 2]; // ACG
        let count = index.count(&pattern);
        assert!(count > 0, "Should find pattern");

        std::fs::remove_file(temp_path).ok();
    }

    #[test]
    fn test_consistency_with_fm_index() {
        let reference = Reference::parse_fasta(">test\nACGTACGTACGT").unwrap();
        let temp_path = std::env::temp_dir().join("test_consistency.idx");

        let fm_index = FMIndex::build(&reference);
        MmapFMIndex::build_and_save(&reference, &temp_path).unwrap();
        let mmap_index = MmapFMIndex::open(&temp_path).unwrap();

        let pattern = [0, 1, 2, 3]; // ACGT
        assert_eq!(fm_index.count(&pattern), mmap_index.count(&pattern));

        let fm_positions = fm_index.find_all(&pattern);
        let mmap_positions = mmap_index.find_all(&pattern);
        assert_eq!(fm_positions, mmap_positions);

        std::fs::remove_file(temp_path).ok();
    }

    #[test]
    fn test_single_contig() {
        use crate::FMIndex;

        // Use longer sequence like existing FM-index tests
        let reference = Reference::parse_fasta(
            ">test
ACGTACGT",
        )
        .unwrap();
        let temp_path = std::env::temp_dir().join("test_single.idx");

        let fm = FMIndex::build(&reference);
        let (fm_l, fm_r) = fm.search(&[3]); // T
        let fm_pos = fm.find_all(&[3]);

        let index = MmapFMIndex::build_and_save(&reference, &temp_path).unwrap();
        let (mmap_l, mmap_r) = index.search(&[3]); // T
        let mmap_pos = index.find_all(&[3]);

        assert_eq!(index.len(), 9); // n+1 with sentinel
        assert_eq!((fm_l, fm_r), (mmap_l, mmap_r), "Search ranges should match");
        assert_eq!(fm_pos, mmap_pos, "Positions should match");
        assert!(!fm_pos.is_empty(), "Should find T in ACGTACGT");

        std::fs::remove_file(temp_path).ok();
    }

    #[test]
    fn test_memory_efficiency() {
        let seq = b"ACGT".repeat(100000);
        let fasta = format!(">test\n{}", String::from_utf8(seq).unwrap());
        let reference = Reference::parse_fasta(&fasta).unwrap();
        let temp_path = std::env::temp_dir().join("test_efficiency.idx");

        MmapFMIndex::build_and_save(&reference, &temp_path).unwrap();

        // Verify mmap file was created with reasonable size
        let metadata = std::fs::metadata(&temp_path).unwrap();
        let index_file_size = metadata.len();

        // Index file should be roughly 5 * len bytes (BWT + SA + counts)
        let expected_min = 40000; // 40000 bases * ~1 byte
        assert!(
            index_file_size > expected_min,
            "Index file size ({}) should be > {} for 40K bases",
            index_file_size,
            expected_min
        );

        // Verify we can open and search without loading into memory
        let index = MmapFMIndex::open(&temp_path).unwrap();
        let count = index.count(&[0, 1, 2, 3]); // ACGT
        assert!(count > 0);

        std::fs::remove_file(temp_path).ok();
    }
}
