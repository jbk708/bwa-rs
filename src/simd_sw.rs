//! Smith-Waterman alignment with SIMD feature detection.
//!
//! Provides CPU feature detection for future SIMD optimization.

use crate::alignment::{Scoring, SeedExtension};
use crate::types::{Cigar, CigarOp};

#[derive(Clone, Copy, Debug)]
pub struct SimdConfig {
    pub enabled: bool,
    pub use_avx512: bool,
    pub use_avx2: bool,
    pub lanes: usize,
}

impl SimdConfig {
    pub fn detect() -> Self {
        #[cfg(any(target_arch = "x86_64", target_arch = "x86"))]
        {
            if is_x86_feature_detected!("avx512f") && is_x86_feature_detected!("avx512bw") {
                return Self {
                    enabled: true,
                    use_avx512: true,
                    use_avx2: false,
                    lanes: 16,
                };
            }
            if is_x86_feature_detected!("avx2") {
                return Self {
                    enabled: true,
                    use_avx512: false,
                    use_avx2: true,
                    lanes: 8,
                };
            }
        }
        Self {
            enabled: false,
            use_avx512: false,
            use_avx2: false,
            lanes: 1,
        }
    }

    pub fn lane_count(&self) -> usize {
        self.lanes
    }
}

impl Default for SimdConfig {
    fn default() -> Self {
        Self::detect()
    }
}

pub fn get_simd_config() -> SimdConfig {
    SimdConfig::detect()
}

#[inline]
pub fn nw_score(query: &[u8], reference: &[u8], scoring: &Scoring) -> i32 {
    scalar_nw_score(query, reference, scoring)
}

fn scalar_nw_score(query: &[u8], reference: &[u8], scoring: &Scoring) -> i32 {
    let (q_len, r_len) = (query.len(), reference.len());
    if q_len == 0 || r_len == 0 {
        return 0;
    }

    let cols = r_len + 1;
    let mut dp_prev = vec![0i32; cols];
    let mut dp_curr = vec![0i32; cols];

    for j in 1..cols {
        dp_prev[j] = dp_prev[j - 1].saturating_sub(scoring.gap_extend);
    }

    for i in 1..=q_len {
        dp_curr[0] = dp_curr[0].saturating_sub(scoring.gap_extend);
        for j in 1..=r_len {
            let is_match = query[i - 1] == reference[j - 1];
            let match_val = if is_match {
                scoring.match_score
            } else {
                -scoring.mismatch_penalty
            };
            let diag = dp_prev[j - 1].saturating_add(match_val);
            let up = dp_curr[j - 1].saturating_sub(scoring.gap_extend);
            let left = dp_prev[j].saturating_sub(scoring.gap_open);
            dp_curr[j] = diag.max(up).max(left);
        }
        std::mem::swap(&mut dp_prev, &mut dp_curr);
    }
    dp_prev[r_len]
}

#[inline]
pub fn extend_forward_simd(
    query: &[u8],
    reference: &[u8],
    start_pos: usize,
    scoring: &Scoring,
    bandwidth: usize,
) -> SeedExtension {
    scalar_extend_forward(query, reference, start_pos, scoring, bandwidth)
}

fn scalar_extend_forward(
    query: &[u8],
    reference: &[u8],
    start_pos: usize,
    scoring: &Scoring,
    bandwidth: usize,
) -> SeedExtension {
    let (q_len, ref_len) = (query.len(), reference.len().saturating_sub(start_pos));
    if q_len == 0 || ref_len == 0 {
        return SeedExtension {
            score: 0,
            query_end: 0,
            ref_end: start_pos,
            cigar: Cigar::new(),
        };
    }

    let bw = bandwidth.clamp(1, 256);
    let mut best_score = 0i32;
    let mut best_i = 0usize;
    let mut best_j = 0usize;

    let max_j = ref_len.min(q_len + bw);
    let cols = max_j + 1;

    let mut dp_prev = vec![i32::MIN / 4; cols];
    let mut dp_curr = vec![i32::MIN / 4; cols];

    dp_curr[0] = 0;
    for j in 1..=max_j {
        dp_curr[j] = dp_curr[j - 1].saturating_sub(scoring.gap_extend);
    }
    std::mem::swap(&mut dp_prev, &mut dp_curr);

    for i in 1..=q_len {
        dp_curr[0] = dp_prev[0].saturating_sub(scoring.gap_extend);
        let j_start = 1.max(i.saturating_sub(bw));
        let j_end = (i + bw).min(max_j);

        for j in j_start..=j_end {
            let ref_idx = start_pos + j - 1;
            if ref_idx >= reference.len() {
                break;
            }
            let is_match = query[i - 1] == reference[ref_idx];
            let match_val = if is_match {
                scoring.match_score
            } else {
                -scoring.mismatch_penalty
            };
            let diag = if j > j_start {
                dp_prev[j - 1].saturating_add(match_val)
            } else {
                match_val
            };
            let up = dp_prev[j].saturating_sub(scoring.gap_open);
            let left = dp_curr[j - 1].saturating_sub(scoring.gap_extend);
            let val = diag.max(up).max(left);
            dp_curr[j] = val;
            if val > best_score {
                best_score = val;
                best_i = i;
                best_j = j;
            }
        }
        std::mem::swap(&mut dp_prev, &mut dp_curr);
    }

    let cigar = traceback(
        &dp_prev, query, reference, start_pos, best_i, best_j, scoring,
    );
    SeedExtension {
        score: best_score,
        query_end: best_i,
        ref_end: start_pos + best_j,
        cigar,
    }
}

fn traceback(
    dp_prev: &[i32],
    query: &[u8],
    reference: &[u8],
    ref_start: usize,
    mut i: usize,
    mut j: usize,
    scoring: &Scoring,
) -> Cigar {
    let mut ops = Vec::new();
    let match_s = scoring.match_score;
    let mismatch_s = -scoring.mismatch_penalty;
    let gap_open = scoring.gap_open;
    let cols = dp_prev.len();

    while i > 0 && j > 0 {
        let rel_j = j - 1;
        if rel_j >= cols {
            break;
        }
        let is_match = query[i - 1] == reference[ref_start + rel_j];
        let delta = if is_match { match_s } else { mismatch_s };
        let curr = dp_prev[rel_j];
        let diag = if rel_j > 0 {
            dp_prev[rel_j - 1]
        } else {
            i32::MIN / 4
        };
        if curr == diag + delta {
            ops.push(if is_match { CigarOp::Eq } else { CigarOp::X });
            i -= 1;
            j -= 1;
        } else if curr == dp_prev[rel_j] - gap_open {
            ops.push(CigarOp::D);
            i -= 1;
        } else {
            ops.push(CigarOp::I);
            j -= 1;
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
    use crate::alignment::Scoring;

    #[test]
    fn test_simd_config_detection() {
        let config = SimdConfig::detect();
        assert!(config.lanes >= 1);
        #[cfg(target_arch = "x86_64")]
        {
            if is_x86_feature_detected!("avx2") {
                assert!(config.use_avx2 || config.use_avx512);
            }
        }
    }

    #[test]
    fn test_scalar_nw_exact_match() {
        let query = vec![0, 1, 2, 3];
        let reference = vec![0, 1, 2, 3];
        let scoring = Scoring::default();
        let score = scalar_nw_score(&query, &reference, &scoring);
        assert_eq!(score, 4);
    }

    #[test]
    fn test_scalar_nw_with_mismatch() {
        let query = vec![0, 1, 2, 3];
        let reference = vec![0, 1, 0, 3];
        let scoring = Scoring::default();
        let _score = scalar_nw_score(&query, &reference, &scoring);
    }

    #[test]
    fn test_nw_score_empty() {
        let scoring = Scoring::default();
        assert_eq!(nw_score(&[], &[0, 1, 2], &scoring), 0);
        assert_eq!(nw_score(&[0, 1, 2], &[], &scoring), 0);
    }

    #[test]
    fn test_extend_forward_empty() {
        let scoring = Scoring::default();
        let result = extend_forward_simd(&[], &[0, 1, 2], 0, &scoring, 16);
        assert_eq!(result.score, 0);
        assert_eq!(result.query_end, 0);
    }

    #[test]
    fn test_extend_forward_exact_match() {
        let query = vec![0, 1, 2, 3];
        let reference = vec![0, 1, 2, 3, 4, 5];
        let scoring = Scoring::default();
        let result = extend_forward_simd(&query, &reference, 0, &scoring, 16);
        assert!(result.score > 0);
        assert!(result.query_end > 0);
    }

    #[test]
    fn test_extend_forward_partial_match() {
        let query = vec![0, 1, 2, 3, 4];
        let reference = vec![0, 1, 2, 5, 6];
        let scoring = Scoring::default();
        let result = extend_forward_simd(&query, &reference, 0, &scoring, 16);
        assert!(result.query_end >= 3);
    }

    #[test]
    fn test_compress_cigar() {
        let ops = vec![
            CigarOp::Eq,
            CigarOp::Eq,
            CigarOp::Eq,
            CigarOp::X,
            CigarOp::X,
        ];
        let cigar = Cigar::compress(ops);
        assert_eq!(cigar.ops.len(), 2);
        assert_eq!(cigar.ops[0], (CigarOp::Eq, 3));
        assert_eq!(cigar.ops[1], (CigarOp::X, 2));
    }

    #[test]
    fn test_compress_cigar_empty() {
        let cigar = Cigar::compress(vec![]);
        assert!(cigar.ops.is_empty());
    }

    #[test]
    fn test_get_simd_config() {
        let config = get_simd_config();
        assert!(config.lanes >= 1);
    }
}
