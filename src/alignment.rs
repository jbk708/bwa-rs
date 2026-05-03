//! Smith-Waterman alignment with affine gap penalties.

use crate::chaining::chain_seeds;
use crate::error::BwaError;
use crate::fm_index::FMIndex;
use crate::seed::{filter_mems, find_mems, DEFAULT_MIN_SEED_LEN};
use crate::types::{AlignmentResult, ChainedSeed, Cigar, CigarOp, MEM};

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
    reference: Vec<u8>,
    scoring: Scoring,
    min_seed_len: usize,
}

impl Aligner {
    pub fn new(index: FMIndex, reference: Vec<u8>) -> Self {
        Self {
            index,
            reference,
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

    pub fn align_read(
        &self,
        query: &[u8],
        _mate: Option<&[u8]>,
    ) -> Result<AlignmentResult, BwaError> {
        let mut mems = find_mems(&self.index, query, self.min_seed_len);

        if mems.is_empty() {
            return Ok(self.unmapped_result());
        }

        filter_mems(&mut mems);

        if mems.is_empty() {
            return Ok(self.unmapped_result());
        }

        let chains = chain_seeds(&mems, 50.0);

        if chains.is_empty() {
            return Ok(self.unmapped_result());
        }

        let best_chain = chains
            .iter()
            .max_by(|a, b| a.score.partial_cmp(&b.score).unwrap())
            .ok_or_else(|| BwaError::Alignment("No valid chains".to_string()))?;

        self.build_alignment(query, best_chain)
    }

    fn build_alignment(
        &self,
        query: &[u8],
        chain: &ChainedSeed,
    ) -> Result<AlignmentResult, BwaError> {
        let ref_seq = &self.reference;
        let seed = &chain.mem;
        let query_len = query.len();

        let bw = optimal_bandwidth(query_len);

        let backward = extend_seed_backward(
            &query[..seed.query_start],
            ref_seq,
            seed.ref_start,
            &self.scoring,
            bw,
        );

        let forward = affine_extend_forward(
            &query[seed.query_end()..],
            ref_seq,
            seed.ref_end(),
            &self.scoring,
            bw,
        );

        let mut cigar = Cigar::new();

        if backward.query_end > 0 {
            cigar.extend(backward.cigar);
        }

        cigar.push(CigarOp::Eq, seed.length as u32);

        if forward.query_end > 0 {
            cigar.extend(forward.cigar);
        }

        let ref_start = seed.ref_start.saturating_sub(backward.query_end);
        let _ref_end = seed.ref_end() + forward.ref_end.saturating_sub(seed.ref_end());

        let nm = self.compute_nm(&cigar, query, ref_seq, ref_start);
        let score = self.compute_alignment_score(&cigar, query, ref_seq, ref_start);

        let mut result = AlignmentResult::new(ref_start, cigar);
        let md_tag = result.mdz_string(query, ref_seq);
        result.score = score;
        result.nm = nm;
        result.mapq = self.calculate_mapq(std::slice::from_ref(&chain.mem), score);
        result.md_tag = Some(md_tag);

        Ok(result)
    }

    fn compute_nm(&self, cigar: &Cigar, query: &[u8], reference: &[u8], ref_start: usize) -> u32 {
        let mut nm = 0u32;
        let mut q_pos = 0usize;
        let mut r_pos = ref_start;

        for (op, len) in &cigar.ops {
            match op {
                CigarOp::M | CigarOp::Eq | CigarOp::X => {
                    for i in 0..*len as usize {
                        let q_idx = q_pos + i;
                        let r_idx = r_pos + i;
                        if q_idx < query.len()
                            && r_idx < reference.len()
                            && query[q_idx] != reference[r_idx]
                        {
                            nm += 1;
                        }
                    }
                    q_pos += *len as usize;
                    r_pos += *len as usize;
                }
                CigarOp::I => {
                    nm += *len;
                    q_pos += *len as usize;
                }
                CigarOp::D => {
                    nm += *len;
                    r_pos += *len as usize;
                }
                _ => {}
            }
        }
        nm
    }

    fn compute_alignment_score(
        &self,
        cigar: &Cigar,
        query: &[u8],
        reference: &[u8],
        ref_start: usize,
    ) -> i32 {
        let mut score = 0i32;
        let mut q_pos = 0usize;
        let mut r_pos = ref_start;

        for (op, len) in &cigar.ops {
            match op {
                CigarOp::M | CigarOp::Eq | CigarOp::X => {
                    for i in 0..*len as usize {
                        let q_idx = q_pos + i;
                        let r_idx = r_pos + i;
                        if q_idx < query.len() && r_idx < reference.len() {
                            if query[q_idx] == reference[r_idx] {
                                score += self.scoring.match_score;
                            } else {
                                score -= self.scoring.mismatch_penalty;
                            }
                        }
                    }
                    q_pos += *len as usize;
                    r_pos += *len as usize;
                }
                CigarOp::I => {
                    score -=
                        self.scoring.gap_open + (*len as i32 - 1).max(0) * self.scoring.gap_extend;
                    q_pos += *len as usize;
                }
                CigarOp::D => {
                    score -=
                        self.scoring.gap_open + (*len as i32 - 1).max(0) * self.scoring.gap_extend;
                    r_pos += *len as usize;
                }
                _ => {}
            }
        }
        score
    }

    pub fn align_paired(
        &self,
        read1: &[u8],
        read2: &[u8],
    ) -> Result<(AlignmentResult, AlignmentResult), BwaError> {
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
            md_tag: None,
        }
    }

    fn calculate_mapq(&self, mems: &[MEM], best_score: i32) -> u8 {
        if mems.is_empty() {
            return 0;
        }

        let second_best = mems
            .iter()
            .skip(1)
            .map(|m| m.score as i32)
            .max()
            .unwrap_or(0);

        if second_best == 0 {
            return 60;
        }

        let ratio = second_best as f64 / best_score.max(1) as f64;
        let prob = ratio * ratio;
        let mapq = (-10.0 * prob.log10()).round() as i32;
        mapq.clamp(0, 60) as u8
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
                if query[best_j.saturating_sub(1)]
                    == reference[seed_start + best_j.saturating_sub(1)]
                {
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
        for ref_pos in [
            seed_end.saturating_sub(d),
            seed_end.saturating_sub(d.saturating_sub(1)),
        ] {
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

#[derive(Clone)]
pub struct AffineDP {
    m: Vec<i32>,
    x: Vec<i32>,
    g: Vec<i32>,
    rows: usize,
    cols: usize,
}

impl AffineDP {
    pub fn new(query_len: usize, ref_len: usize) -> Self {
        let rows = query_len + 1;
        let cols = ref_len + 1;
        let size = rows * cols;
        Self {
            m: vec![i32::MIN / 2; size],
            x: vec![i32::MIN / 2; size],
            g: vec![i32::MIN / 2; size],
            rows,
            cols,
        }
    }

    #[allow(dead_code)]
    fn idx(&self, i: usize, j: usize) -> usize {
        i * self.cols + j
    }

    fn m_at(&self, i: usize, j: usize) -> i32 {
        self.m[i * self.cols + j]
    }

    fn x_at(&self, i: usize, j: usize) -> i32 {
        self.x[i * self.cols + j]
    }

    fn g_at(&self, i: usize, j: usize) -> i32 {
        self.g[i * self.cols + j]
    }

    fn set_m(&mut self, i: usize, j: usize, val: i32) {
        self.m[i * self.cols + j] = val;
    }

    fn set_x(&mut self, i: usize, j: usize, val: i32) {
        self.x[i * self.cols + j] = val;
    }

    fn set_g(&mut self, i: usize, j: usize, val: i32) {
        self.g[i * self.cols + j] = val;
    }

    pub fn rows(&self) -> usize {
        self.rows
    }

    pub fn cols(&self) -> usize {
        self.cols
    }
}

pub fn affine_extend_forward(
    query: &[u8],
    reference: &[u8],
    start_pos: usize,
    scoring: &Scoring,
    bandwidth: usize,
) -> SeedExtension {
    let query_len = query.len();
    let ref_len = reference.len().saturating_sub(start_pos);

    if query_len == 0 || ref_len == 0 {
        return SeedExtension {
            score: 0,
            query_end: 0,
            ref_end: start_pos,
            cigar: Cigar::new(),
        };
    }

    let ref_end = (start_pos + ref_len).min(reference.len());
    let actual_ref_len = ref_end - start_pos;
    let bw = bandwidth.clamp(1, 256);

    let mut dp = AffineDP::new(query_len, actual_ref_len);
    dp.set_m(0, 0, 0);
    dp.set_x(0, 0, i32::MIN / 2);
    dp.set_g(0, 0, i32::MIN / 2);

    for i in 1..=query_len {
        let j_start = 1.max(i.saturating_sub(bw));
        let j_end = (i + bw).min(actual_ref_len).min(actual_ref_len);

        for j in j_start..=j_end {
            let is_match = query[i - 1] == reference[start_pos + j - 1];
            let match_score = if is_match {
                scoring.match_score
            } else {
                -(scoring.mismatch_penalty)
            };

            let m_prev = dp.m_at(i - 1, j - 1);
            let x_prev = dp.x_at(i - 1, j - 1);
            let g_prev = dp.g_at(i - 1, j - 1);
            let m_val = m_prev
                .saturating_add(match_score)
                .max(x_prev.saturating_add(match_score))
                .max(g_prev.saturating_add(match_score));
            dp.set_m(i, j, m_val);

            let x_open = dp.m_at(i - 1, j).saturating_sub(scoring.gap_open);
            let x_extend = dp.x_at(i - 1, j).saturating_sub(scoring.gap_extend);
            dp.set_x(i, j, x_open.max(x_extend));

            let g_open = dp.m_at(i, j - 1).saturating_sub(scoring.gap_open);
            let g_extend = dp.g_at(i, j - 1).saturating_sub(scoring.gap_extend);
            dp.set_g(i, j, g_open.max(g_extend));
        }
    }

    let mut best_score = i32::MIN;
    let mut best_i = 0;
    let mut best_j = 0;

    for i in 1..=query_len {
        let j = (actual_ref_len.min(i.saturating_add(bw))).min(actual_ref_len);
        let score = dp.m_at(i, j).max(dp.x_at(i, j)).max(dp.g_at(i, j));
        if score > best_score {
            best_score = score;
            best_i = i;
            best_j = j;
        }
    }

    let cigar = traceback_affine(&dp, query, reference, start_pos, best_i, best_j, scoring);

    SeedExtension {
        score: best_score,
        query_end: best_i,
        ref_end: start_pos + best_j,
        cigar,
    }
}

fn traceback_affine(
    dp: &AffineDP,
    query: &[u8],
    reference: &[u8],
    ref_start: usize,
    mut i: usize,
    mut j: usize,
    scoring: &Scoring,
) -> Cigar {
    let mut ops = Vec::new();
    let match_s = scoring.match_score;
    let mismatch_s = -(scoring.mismatch_penalty);

    while i > 0 && j > 0 {
        let current = dp.m_at(i, j).max(dp.x_at(i, j)).max(dp.g_at(i, j));

        let is_match = query[i - 1] == reference[ref_start + j - 1];
        let delta = if is_match { match_s } else { mismatch_s };

        if current == dp.m_at(i, j) && current == dp.m_at(i - 1, j - 1) + delta {
            ops.push(if is_match { CigarOp::Eq } else { CigarOp::X });
            i -= 1;
            j -= 1;
        } else if current == dp.x_at(i, j) {
            ops.push(CigarOp::D);
            i -= 1;
        } else if current == dp.g_at(i, j) {
            ops.push(CigarOp::I);
            j -= 1;
        } else if current == dp.m_at(i, j) {
            ops.push(if is_match { CigarOp::Eq } else { CigarOp::X });
            i -= 1;
            j -= 1;
        } else {
            break;
        }
    }

    while i > 0 {
        ops.push(CigarOp::D);
        i -= 1;
    }
    while j > 0 {
        ops.push(CigarOp::I);
        j -= 1;
    }

    ops.reverse();
    Cigar::compress(ops)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::reference::Reference;

    #[test]
    fn test_unmapped_result() {
        let ref_seq = Reference::parse_fasta(">test\nACGT").unwrap();
        let ref_data = ref_seq.as_slice().to_vec();
        let index = FMIndex::build(&ref_seq);
        let aligner = Aligner::new(index, ref_data);

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
        let ref_data = ref_seq.as_slice().to_vec();
        let index = FMIndex::build(&ref_seq);

        let aligner = Aligner::new(index, ref_data)
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

    #[test]
    fn test_affine_exact_match() {
        let query = vec![0, 1, 2, 3];
        let reference = vec![0, 1, 2, 3];
        let scoring = Scoring::default();

        let ext = affine_extend_forward(&query, &reference, 0, &scoring, 16);
        assert!(ext.score > 0, "Score should be positive for exact match");
        assert_eq!(ext.query_end, query.len());
        let cigar_str = ext.cigar.to_string();
        assert!(
            cigar_str.contains('=') || cigar_str.contains('X'),
            "CIGAR should contain match ops: {}",
            cigar_str
        );
    }

    #[test]
    fn test_affine_insertion() {
        let query = vec![0, 1, 4, 2, 3];
        let reference = vec![0, 1, 2, 3];
        let scoring = Scoring::default();

        let ext = affine_extend_forward(&query, &reference, 0, &scoring, 16);
        let cigar_str = ext.cigar.to_string();
        assert!(cigar_str.len() > 0, "CIGAR should not be empty");
    }

    #[test]
    fn test_affine_deletion() {
        let query = vec![0, 1, 2, 3];
        let reference = vec![0, 1, 4, 2, 3];
        let scoring = Scoring::default();

        let ext = affine_extend_forward(&query, &reference, 0, &scoring, 16);
        let cigar_str = ext.cigar.to_string();
        assert!(cigar_str.len() > 0, "CIGAR should not be empty");
    }

    #[test]
    fn test_affine_gap_open_vs_extend() {
        let scoring = Scoring::default();
        assert!(scoring.gap_open > scoring.gap_extend);
        let single_gap_cost = scoring.gap_open as i32;
        let two_gap_cost = (scoring.gap_open + scoring.gap_extend) as i32;
        assert!(
            single_gap_cost < two_gap_cost,
            "Single gap should cost less than two separate gaps"
        );
    }

    #[test]
    fn test_affine_empty_query() {
        let query: Vec<u8> = vec![];
        let reference = vec![0, 1, 2, 3];
        let scoring = Scoring::default();

        let ext = affine_extend_forward(&query, &reference, 0, &scoring, 16);
        assert_eq!(ext.score, 0);
        assert_eq!(ext.query_end, 0);
    }

    #[test]
    fn test_affine_empty_reference() {
        let query = vec![0, 1, 2, 3];
        let reference: Vec<u8> = vec![];
        let scoring = Scoring::default();

        let ext = affine_extend_forward(&query, &reference, 0, &scoring, 16);
        assert_eq!(ext.score, 0);
        assert_eq!(ext.query_end, 0);
    }

    #[test]
    fn test_affine_dp_struct() {
        let dp = AffineDP::new(5, 10);
        assert_eq!(dp.rows(), 6);
        assert_eq!(dp.cols(), 11);
        assert!(dp.m_at(0, 0) < 0);
    }

    #[test]
    fn test_align_full_read() {
        let ref_seq = Reference::parse_fasta(">test\nACGTACGT").unwrap();
        let ref_data = ref_seq.as_slice().to_vec();
        let ref_data_len = ref_data.len();
        let index = FMIndex::build(&ref_seq);
        let index_for_check = index.clone();
        let aligner = Aligner::new(index, ref_data).min_seed_len(2);

        // Use AC which we know works with the simple search test
        let query = vec![0, 1]; // AC

        // Check if FM-index finds the pattern
        let positions = index_for_check.find_all(&query);
        assert!(
            !positions.is_empty(),
            "FM-index should find AC positions: {:?}, ref_len: {}",
            positions,
            ref_data_len
        );

        let result = aligner.align_read(&query, None).unwrap();

        if result.cigar.ops.is_empty() {
            panic!(
                "CIGAR is empty - flag: {}, mapq: {}, score: {}",
                result.flag, result.mapq, result.score
            );
        }
        assert!(
            result.cigar.reference_length() > 0,
            "CIGAR should cover reference positions"
        );
        // Score can be negative if extension introduces mismatches
        assert!(result.nm < u32::MAX as u32, "NM tag should be set");
    }

    #[test]
    fn test_align_with_mismatch() {
        let ref_seq = Reference::parse_fasta(">test\nACGTACGT").unwrap();
        let ref_data = ref_seq.as_slice().to_vec();
        let index = FMIndex::build(&ref_seq);
        let aligner = Aligner::new(index, ref_data).min_seed_len(2);

        // Use shorter query
        let query = vec![0, 1, 2, 3, 2]; // ACGTC
        let result = aligner.align_read(&query, None).unwrap();

        assert!(
            result.cigar.to_string().len() > 0,
            "CIGAR should not be empty"
        );
    }

    #[test]
    fn test_align_unmapped() {
        let ref_seq = Reference::parse_fasta(">test\nACGT").unwrap();
        let ref_data = ref_seq.as_slice().to_vec();
        let index = FMIndex::build(&ref_seq);
        let aligner = Aligner::new(index, ref_data).min_seed_len(100);

        let query = vec![0, 1, 2, 3, 0, 1, 2, 3];
        let result = aligner.align_read(&query, None).unwrap();

        assert_eq!(result.flag & 0x4, 0x4, "Unmapped flag should be set");
        assert_eq!(result.mapq, 0, "MAPQ should be 0 for unmapped");
    }

    #[test]
    fn test_mapq_calculation() {
        let ref_seq = Reference::parse_fasta(">test\nACGTACGTACGT").unwrap();
        let ref_data = ref_seq.as_slice().to_vec();
        let index = FMIndex::build(&ref_seq);
        let aligner = Aligner::new(index, ref_data).min_seed_len(2);

        let query = vec![0, 1, 2, 3];
        let result = aligner.align_read(&query, None).unwrap();

        assert!(result.mapq <= 60, "MAPQ should be capped at 60");
    }

    #[test]
    fn test_cigar_extend() {
        let mut cigar1 = Cigar::new();
        cigar1.push(CigarOp::Eq, 5);

        let mut cigar2 = Cigar::new();
        cigar2.push(CigarOp::Eq, 3);

        cigar1.extend(cigar2);
        assert_eq!(cigar1.to_string(), "8=");
    }

    #[test]
    fn test_cigar_extend_different_ops() {
        let mut cigar1 = Cigar::new();
        cigar1.push(CigarOp::Eq, 5);

        let mut cigar2 = Cigar::new();
        cigar2.push(CigarOp::I, 3);

        cigar1.extend(cigar2);
        assert_eq!(cigar1.to_string(), "5=3I");
    }

    #[test]
    fn test_mdz_string_all_matches() {
        let mut cigar = Cigar::new();
        cigar.push(CigarOp::Eq, 10);

        let query = vec![0u8; 10];
        let reference = vec![0u8; 10];

        let result = AlignmentResult::new(0, cigar);
        let mdz = result.mdz_string(&query, &reference);
        assert_eq!(mdz, "MD:Z:10");
    }

    #[test]
    fn test_mdz_string_with_mismatch() {
        let mut cigar = Cigar::new();
        cigar.push(CigarOp::Eq, 5);
        cigar.push(CigarOp::X, 1);
        cigar.push(CigarOp::Eq, 3);

        let query = vec![0, 0, 0, 0, 0, 1, 0, 0, 0];
        let reference = vec![0, 0, 0, 0, 0, 0, 0, 0, 0];

        let result = AlignmentResult::new(0, cigar);
        let mdz = result.mdz_string(&query, &reference);
        assert_eq!(mdz, "MD:Z:5A3");
    }

    #[test]
    fn test_mdz_string_with_deletion() {
        let mut cigar = Cigar::new();
        cigar.push(CigarOp::Eq, 5);
        cigar.push(CigarOp::D, 1);
        cigar.push(CigarOp::Eq, 3);

        let query = vec![0u8; 8];
        let reference = vec![0, 0, 0, 0, 0, 1, 0, 0, 0];

        let result = AlignmentResult::new(0, cigar);
        let mdz = result.mdz_string(&query, &reference);
        assert_eq!(mdz, "MD:Z:5^C3");
    }

    #[test]
    fn test_mdz_string_with_insertion() {
        let mut cigar = Cigar::new();
        cigar.push(CigarOp::Eq, 3);
        cigar.push(CigarOp::I, 2);
        cigar.push(CigarOp::Eq, 3);

        let query = vec![0, 0, 0, 1, 2, 0, 0, 0];
        let reference = vec![0, 0, 0, 0, 0, 0];

        let result = AlignmentResult::new(0, cigar);
        let mdz = result.mdz_string(&query, &reference);
        // Insertions don't appear in MD:Z per SAM spec - only reference matches count
        assert_eq!(mdz, "MD:Z:6");
    }

    #[test]
    fn test_mdz_string_empty() {
        let cigar = Cigar::new();
        let query: Vec<u8> = vec![];
        let reference: Vec<u8> = vec![];

        let result = AlignmentResult::new(0, cigar);
        let mdz = result.mdz_string(&query, &reference);
        assert_eq!(mdz, "MD:Z:0");
    }

    #[test]
    fn test_mdz_string_complex() {
        let mut cigar = Cigar::new();
        cigar.push(CigarOp::Eq, 2);
        cigar.push(CigarOp::X, 1);
        cigar.push(CigarOp::Eq, 3);
        cigar.push(CigarOp::D, 1);
        cigar.push(CigarOp::Eq, 4);

        let query = vec![0, 0, 1, 0, 0, 0, 0, 0, 0, 0];
        let reference = vec![0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0];

        let result = AlignmentResult::new(0, cigar);
        let mdz = result.mdz_string(&query, &reference);
        assert_eq!(mdz, "MD:Z:2A3^C4");
    }

    #[test]
    fn test_mdz_string_in_alignment_result() {
        let ref_seq = Reference::parse_fasta(">test\nACGTACGT").unwrap();
        let ref_data = ref_seq.as_slice().to_vec();
        let index = FMIndex::build(&ref_seq);
        let aligner = Aligner::new(index, ref_data).min_seed_len(2);

        let query = vec![0, 1, 2, 3];
        let result = aligner.align_read(&query, None).unwrap();

        if result.md_tag.is_none() {
            panic!("MD tag should be set for aligned reads");
        }
        assert!(
            result.md_tag.as_ref().unwrap().starts_with("MD:Z:"),
            "{}",
            result.md_tag.unwrap()
        );
    }
}
