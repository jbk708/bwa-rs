//! FM-Index implementation for fast substring search.

use crate::reference::Reference;

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

    pub fn get(&self, idx: usize) -> u8 {
        self.bwt[idx]
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
}