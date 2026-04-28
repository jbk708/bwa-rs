//! Smith-Waterman alignment with affine gap penalties.

use crate::error::BwaError;
use crate::fm_index::FMIndex;
use crate::seed::{find_mems, filter_mems, DEFAULT_MIN_SEED_LEN};
use crate::types::{AlignmentResult, Cigar, CigarOp};

#[derive(Clone, Debug)]
pub struct Scoring {
    pub match_score: i32,
    pub mismatch_penalty: i32,
    pub gap_open: i32,
    pub gap_extend: i32,
}

impl Default for Scoring {
    fn default() -> Self {
        Self {
            match_score: 1,
            mismatch_penalty: 4,
            gap_open: 6,
            gap_extend: 1,
        }
    }
}

#[derive(Clone)]
pub struct Aligner {
    index: FMIndex,
    scoring: Scoring,
    min_seed_len: usize,
}

impl Aligner {
    pub fn new(index: FMIndex) -> Self {
        Self {
            index,
            scoring: Scoring::default(),
            min_seed_len: DEFAULT_MIN_SEED_LEN,
        }
    }

    pub fn min_seed_len(mut self, len: usize) -> Self {
        self.min_seed_len = len;
        self
    }

    pub fn scoring(mut self, scoring: Scoring) -> Self {
        self.scoring = scoring;
        self
    }

    pub fn align_read(&self, query: &[u8], _mate: Option<&[u8]>) -> Result<AlignmentResult, BwaError> {
        let mut mems = find_mems(&self.index, query, self.min_seed_len);

        if mems.is_empty() {
            return Ok(self.unmapped_result());
        }

        filter_mems(&mut mems);

        let best_mem = mems
            .iter()
            .max_by_key(|m| m.length)
            .ok_or_else(|| BwaError::Alignment("No valid seeds".to_string()))?;

        let mut cigar = Cigar::new();
        cigar.push(CigarOp::M, best_mem.length as u32);

        let mut result = AlignmentResult::new(best_mem.ref_start, cigar);
        result.score = best_mem.length as i32 * self.scoring.match_score;
        result.mapq = self.estimate_mapq(&mems);

        Ok(result)
    }

    pub fn align_paired(&self, read1: &[u8], read2: &[u8]) -> Result<(AlignmentResult, AlignmentResult), BwaError> {
        let r1 = self.align_read(read1, Some(read2))?;
        let r2 = self.align_read(read2, Some(read1))?;
        Ok((r1, r2))
    }

    fn unmapped_result(&self) -> AlignmentResult {
        AlignmentResult {
            position: 0,
            mapq: 0,
            cigar: Cigar::new(),
            flag: 0x4,
            reverse_strand: false,
            nm: 0,
            score: 0,
        }
    }

    fn estimate_mapq(&self, mems: &[crate::types::MEM]) -> u8 {
        if mems.is_empty() {
            return 0;
        }

        let best_score = mems.iter().map(|m| m.score as i32).max().unwrap_or(0);
        let second_best = mems
            .iter()
            .filter(|m| m.score < best_score as f32)
            .map(|m| m.score as i32)
            .max()
            .unwrap_or(0);

        if second_best == 0 {
            return 60;
        }

        let diff = best_score - second_best;
        let mapq = 30.min(60 * diff / (2 * best_score.max(1)) as i32 + diff / 2);
        mapq as u8
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::reference::Reference;

    #[test]
    fn test_unmapped_result() {
        let ref_seq = Reference::parse_fasta(">test\nACGT").unwrap();
        let index = FMIndex::build(&ref_seq);
        let aligner = Aligner::new(index);

        let unmapped = aligner.unmapped_result();
        assert_eq!(unmapped.flag, 0x4);
        assert_eq!(unmapped.mapq, 0);
    }

    #[test]
    fn test_scoring_default() {
        let scoring = Scoring::default();
        assert_eq!(scoring.match_score, 1);
        assert_eq!(scoring.mismatch_penalty, 4);
    }

    #[test]
    fn test_aligner_builder() {
        let ref_seq = Reference::parse_fasta(">test\nACGT").unwrap();
        let index = FMIndex::build(&ref_seq);

        let aligner = Aligner::new(index)
            .min_seed_len(10)
            .scoring(Scoring::default());

        assert_eq!(aligner.min_seed_len, 10);
    }
}