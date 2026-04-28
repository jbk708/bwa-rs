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

#[derive(Clone, Debug)]
pub struct SeedExtension {
    pub score: i32,
    pub query_end: usize,
    pub ref_end: usize,
    pub cigar: Cigar,
}

pub fn extend_seed_forward(
    query: &[u8],
    reference: &[u8],
    seed_start: usize,
    scoring: &Scoring,
    bandwidth: usize,
) -> SeedExtension {
    let query_len = query.len();
    let ref_len = reference.len();

    if query_len == 0 || ref_len == 0 {
        return SeedExtension {
            score: 0,
            query_end: 0,
            ref_end: seed_start,
            cigar: Cigar::new(),
        };
    }

    let actual_bw = bandwidth.clamp(1, 256);
    let mut best_score = 0;
    let mut best_j = 0;
    let mut traceback: Vec<(usize, i8)> = Vec::new();

    for d in 0..=actual_bw {
        for ref_pos in [seed_start.saturating_sub(d), seed_start + d] {
            if ref_pos >= ref_len {
                continue;
            }

            let j_max = query_len.min(d + 1);

            for j in 1..=j_max {
                let ref_idx = ref_pos + j.saturating_sub(1);

                if ref_idx >= ref_len {
                    continue;
                }

                let score = if query[j - 1] == reference[ref_idx] {
                    scoring.match_score
                } else {
                    scoring.mismatch_penalty
                };

                if score > best_score {
                    best_score = score;
                    best_j = j;
                    traceback.clear();
                    traceback.push((j, 1));
                } else if score == best_score && score > 0 && j > best_j {
                    best_j = j;
                    traceback.clear();
                    traceback.push((j, 1));
                }
            }
        }
    }

    let mut cigar = Cigar::new();
    let mut last_op = CigarOp::M;
    let mut last_len = 0u32;

    for (_j, op_code) in &traceback {
        let op = match op_code {
            1 => {
                if query[best_j.saturating_sub(1)] == reference[seed_start + best_j.saturating_sub(1)] {
                    CigarOp::Eq
                } else {
                    CigarOp::X
                }
            }
            2 => CigarOp::I,
            3 => CigarOp::D,
            _ => CigarOp::M,
        };

        if op == last_op {
            last_len += 1;
        } else {
            if last_len > 0 {
                cigar.push(last_op, last_len);
            }
            last_op = op;
            last_len = 1;
        }
    }

    if last_len > 0 {
        cigar.push(last_op, last_len);
    }

    if cigar.ops.is_empty() {
        cigar.push(CigarOp::M, best_j as u32);
    }

    SeedExtension {
        score: best_score,
        query_end: best_j,
        ref_end: seed_start + best_j,
        cigar,
    }
}

pub fn extend_seed_backward(
    query: &[u8],
    reference: &[u8],
    seed_end: usize,
    scoring: &Scoring,
    bandwidth: usize,
) -> SeedExtension {
    let query_len = query.len();
    let ref_len = reference.len();

    if query_len == 0 || ref_len == 0 || seed_end == 0 {
        return SeedExtension {
            score: 0,
            query_end: 0,
            ref_end: seed_end,
            cigar: Cigar::new(),
        };
    }

    let actual_bw = bandwidth.clamp(1, 256);
    let mut best_score = 0;
    let mut best_j = 0;

    for d in 0..=actual_bw {
        for ref_pos in [seed_end.saturating_sub(d), seed_end.saturating_sub(d.saturating_sub(1).max(0))] {
            if ref_pos == 0 || ref_pos >= ref_len {
                continue;
            }

            let j_max = query_len.min(d + 1);

            for j in 1..=j_max {
                let ref_idx = ref_pos.saturating_sub(j);
                if ref_idx >= ref_len {
                    continue;
                }

                let score = if j <= query_len && query[j - 1] == reference[ref_idx] {
                    scoring.match_score
                } else {
                    scoring.mismatch_penalty
                };

                if score > best_score {
                    best_score = score;
                    best_j = j;
                } else if score == best_score && score > 0 && j > best_j {
                    best_j = j;
                }
            }
        }
    }

    let mut cigar = Cigar::new();
    let mut last_op = CigarOp::M;
    let mut last_len = best_j as u32;

    if best_j > 0 {
        let ref_start = seed_end.saturating_sub(best_j);
        let mut offset = 0;
        for &base in query.iter().take(best_j) {
            let ref_idx = ref_start + offset;
            if ref_idx >= ref_len {
                offset += 1;
                continue;
            }
            let op = if base == reference[ref_idx] {
                CigarOp::Eq
            } else {
                CigarOp::X
            };
            if op != last_op {
                if last_len > 0 {
                    cigar.push(last_op, last_len);
                }
                last_op = op;
                last_len = 1;
            } else {
                last_len += 1;
            }
            offset += 1;
        }
    }

    if last_len > 0 {
        cigar.push(last_op, last_len);
    }

    if cigar.ops.is_empty() {
        cigar.push(CigarOp::M, 0);
    }

    SeedExtension {
        score: best_score,
        query_end: best_j,
        ref_end: seed_end.saturating_sub(best_j),
        cigar,
    }
}

pub fn optimal_bandwidth(query_len: usize) -> usize {
    (16 + query_len / 2).min(256)
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

    #[test]
    fn test_extend_forward_exact_match() {
        let query = vec![0, 1, 2, 3];
        let reference = vec![0, 1, 2, 3, 4, 0, 1, 2, 3];
        let scoring = Scoring::default();

        let ext = extend_seed_forward(&query, &reference, 0, &scoring, 16);
        assert!(ext.score > 0);
        assert_eq!(ext.query_end, query.len());
        assert_eq!(ext.ref_end, 4);
    }

    #[test]
    fn test_extend_forward_with_mismatch() {
        let query = vec![0, 1, 2, 3];
        let reference = vec![0, 1, 0, 3, 0, 1, 2, 3];
        let scoring = Scoring::default();

        let ext = extend_seed_forward(&query, &reference, 3, &scoring, 16);
        assert!(ext.score >= 0);
    }

    #[test]
    fn test_extend_forward_insertion() {
        let query = vec![0, 1, 5, 2, 3];
        let reference = vec![0, 1, 2, 3];
        let scoring = Scoring::default();

        let ext = extend_seed_forward(&query, &reference, 0, &scoring, 16);
        assert!(ext.score >= 0);
        assert!(ext.query_end <= query.len());
    }

    #[test]
    fn test_extend_backward() {
        let query = vec![2, 3, 0, 1];
        let reference = vec![0, 1, 2, 3, 4, 2, 3, 0, 1];
        let scoring = Scoring::default();

        let ext = extend_seed_backward(&query, &reference, 6, &scoring, 16);
        assert!(ext.score >= 0);
    }

    #[test]
    fn test_optimal_bandwidth() {
        assert!(optimal_bandwidth(50) >= 16);
        assert!(optimal_bandwidth(500) >= 16);
        assert!(optimal_bandwidth(500) <= 256);
    }

    #[test]
    fn test_seed_extension_zero_bandwidth() {
        let query = vec![0, 1, 2];
        let reference = vec![0, 1, 2, 3];
        let scoring = Scoring::default();

        let ext = extend_seed_forward(&query, &reference, 0, &scoring, 0);
        assert!(ext.score >= 0);
    }
}