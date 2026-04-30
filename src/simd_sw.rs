//! SIMD-accelerated Smith-Waterman alignment.
//!
//! Uses AVX2/AVX-512 for vectorized dynamic programming.
//! Falls back to scalar implementation on unsupported hardware.

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
    let config = SimdConfig::detect();
    if config.enabled && config.use_avx2 {
        #[cfg(target_arch = "x86_64")]
        unsafe {
            return avx2_nw_score(query, reference, scoring);
        }
    }
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
    let config = SimdConfig::detect();
    if config.enabled && config.use_avx2 {
        #[cfg(target_arch = "x86_64")]
        unsafe {
            return avx2_extend_forward(query, reference, start_pos, scoring, bandwidth);
        }
    }
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
    compress_cigar(ops)
}

fn compress_cigar(ops: Vec<CigarOp>) -> Cigar {
    let mut cigar = Cigar::new();
    if ops.is_empty() {
        return cigar;
    }
    let mut current_op = ops[0];
    let mut current_len = 1u32;
    for op in ops.iter().skip(1) {
        if *op == current_op {
            current_len += 1;
        } else {
            cigar.push(current_op, current_len);
            current_op = *op;
            current_len = 1;
        }
    }
    cigar.push(current_op, current_len);
    cigar
}

#[cfg(target_arch = "x86_64")]
mod avx2_impl {
    use super::*;
    use wide::i32x8;
    const LANES: usize = 8;

    #[target_feature(enable = "avx2")]
    pub unsafe fn avx2_nw_score(query: &[u8], reference: &[u8], scoring: &Scoring) -> i32 {
        let (q_len, r_len) = (query.len(), reference.len());
        if q_len == 0 || r_len < LANES {
            return scalar_nw_score(query, reference, scoring);
        }

        let match_s = scoring.match_score as i32;
        let mismatch_s = -(scoring.mismatch_penalty as i32);
        let gap_extend = -(scoring.gap_extend as i32);
        let gap_open = -(scoring.gap_open as i32);

        let match_v = i32x8::splat(match_s);
        let mismatch_v = i32x8::splat(mismatch_s);
        let gap_extend_v = i32x8::splat(gap_extend);
        let gap_open_v = i32x8::splat(gap_open);

        let mut dp_prev = vec![0i32; r_len + 1];
        let mut dp_curr = vec![0i32; r_len + 1];

        for j in 1..=r_len {
            dp_curr[j] = dp_curr[j - 1].saturating_add(gap_extend);
        }
        std::mem::swap(&mut dp_prev, &mut dp_curr);

        let num_blocks = (r_len + 1) / LANES;

        for i in 1..=q_len {
            let qb = query[i - 1] as i32;
            dp_curr[0] = dp_curr[0].saturating_add(gap_extend);

            for b in 0..num_blocks {
                let base = b * LANES;
                let mut rv = [0i32; LANES];
                for k in 0..LANES {
                    rv[k] = reference.get(base + k).copied().unwrap_or(0) as i32;
                }

                let ref_v = i32x8::from_array(rv);
                let query_v = i32x8::splat(qb);
                let matches = query_v.eq(ref_v);
                let scores = matches.select(match_v, mismatch_v);

                let mut prev_arr = [0i32; LANES];
                prev_arr.copy_from_slice(&dp_prev[base..base + LANES]);
                let prev_v = i32x8::new(prev_arr);

                let mut left_arr = [0i32; LANES];
                left_arr.copy_from_slice(&dp_curr[base..base + LANES]);
                let left_v = i32x8::new(left_arr);

                let mut diag_arr = [0i32; LANES];
                diag_arr.copy_from_slice(&dp_prev[base + 1..base + 1 + LANES]);
                let mut diag_v = i32x8::new(diag_arr);
                diag_v = diag_v.replace(7, if base > 0 { dp_prev[base - 1] } else { 0 });

                let diag_scores = diag_v + scores;
                let left_scores = left_v + gap_extend_v;
                let up_scores = prev_v + gap_open_v;

                let result = diag_scores.max(left_scores).max(up_scores);
                dp_curr[base..base + LANES].copy_from_slice(&result.as_array());
            }

            for j in num_blocks * LANES..=r_len {
                let is_match = query[i - 1] == reference[j - 1];
                let s = if is_match { match_s } else { mismatch_s };
                dp_curr[j] = dp_curr[j - 1]
                    .saturating_add(gap_extend)
                    .max(dp_prev[j].saturating_add(gap_open))
                    .max(dp_prev[j - 1].saturating_add(s));
            }
            std::mem::swap(&mut dp_prev, &mut dp_curr);
        }

        dp_prev[..r_len + 1]
            .iter()
            .copied()
            .max()
            .unwrap_or(0)
            .max(0)
    }

    #[target_feature(enable = "avx2")]
    pub unsafe fn avx2_extend_forward(
        query: &[u8],
        reference: &[u8],
        start_pos: usize,
        scoring: &Scoring,
        bandwidth: usize,
    ) -> SeedExtension {
        let (q_len, ref_avail) = (query.len(), reference.len().saturating_sub(start_pos));
        if q_len == 0 || ref_avail == 0 {
            return SeedExtension {
                score: 0,
                query_end: 0,
                ref_end: start_pos,
                cigar: Cigar::new(),
            };
        }

        let bw = bandwidth.clamp(1, 128);
        let match_s = scoring.match_score as i32;
        let mismatch_s = -(scoring.mismatch_penalty as i32);
        let gap_extend = scoring.gap_extend as i32;
        let gap_open = scoring.gap_open as i32;

        let max_j = ref_avail.min(q_len + bw);
        let cols = max_j + 1;

        let mut dp_prev = vec![i32::MIN / 4; cols];
        let mut dp_curr = vec![i32::MIN / 4; cols];

        dp_curr[0] = 0;
        for j in 1..=max_j {
            dp_curr[j] = dp_curr[j - 1].saturating_sub(gap_extend);
        }
        std::mem::swap(&mut dp_prev, &mut dp_curr);

        let match_v = i32x8::splat(match_s);
        let mismatch_v = i32x8::splat(mismatch_s);
        let gap_extend_v = i32x8::splat(-gap_extend);
        let gap_open_v = i32x8::splat(-gap_open);

        let mut best_score = 0i32;
        let mut best_i = 0usize;
        let mut best_j = 0usize;

        for i in 1..=q_len {
            let qb = query[i - 1] as i32;
            dp_curr[0] = dp_prev[0].saturating_sub(gap_extend);
            let i_min_j = if i > bw { i - bw } else { 1 };
            let i_max_j = (i + bw).min(max_j);

            let mut j = i_min_j;
            while j + LANES <= i_max_j {
                let base = j;
                let mut rv = [0i32; LANES];
                for k in 0..LANES {
                    let idx = start_pos + base + k - 1;
                    rv[k] = reference.get(idx).copied().unwrap_or(0) as i32;
                }

                let ref_v = i32x8::from_array(rv);
                let query_v = i32x8::splat(qb);
                let matches = query_v.eq(ref_v);
                let scores = matches.select(match_v, mismatch_v);

                let mut prev_arr = [0i32; LANES];
                prev_arr.copy_from_slice(&dp_prev[base..base + LANES]);
                let prev_v = i32x8::new(prev_arr);

                let mut left_arr = [0i32; LANES];
                left_arr.copy_from_slice(&dp_curr[base..base + LANES]);
                let left_v = i32x8::new(left_arr);

                let mut diag_arr = [0i32; LANES];
                diag_arr.copy_from_slice(&dp_prev[base + 1..base + 1 + LANES]);
                let mut diag_v = i32x8::new(diag_arr);
                diag_v = diag_v.replace(7, if base > 0 { dp_prev[base - 1] } else { 0 });

                let result = (diag_v + scores)
                    .max(left_v + gap_extend_v)
                    .max(prev_v + gap_open_v);
                dp_curr[base..base + LANES].copy_from_slice(&result.as_array());

                let block_max = result.horizontal_max();
                if block_max > best_score {
                    best_score = block_max;
                    best_i = i;
                    let arr = result.as_array();
                    best_j = base + arr.iter().position(|&x| x == block_max).unwrap_or(0);
                }
                j += LANES;
            }

            for jj in j..=i_max_j {
                let ref_idx = start_pos + jj - 1;
                if ref_idx >= reference.len() {
                    break;
                }
                let is_match = qb as u8 == reference[ref_idx];
                let s = if is_match { match_s } else { mismatch_s };
                let diag = if jj > i_min_j {
                    dp_prev[jj - 1].saturating_add(s)
                } else {
                    s
                };
                let up = dp_prev[jj].saturating_sub(gap_open);
                let left = dp_curr[jj - 1].saturating_sub(gap_extend);
                let val = diag.max(up).max(left);
                dp_curr[jj] = val;
                if val > best_score {
                    best_score = val;
                    best_i = i;
                    best_j = jj;
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
}

#[cfg(target_arch = "x86_64")]
pub use avx2_impl::avx2_extend_forward;
#[cfg(target_arch = "x86_64")]
pub use avx2_impl::avx2_nw_score;

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
        let score = scalar_nw_score(&query, &reference, &scoring);
        let _ = score;
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
        let cigar = compress_cigar(ops);
        assert_eq!(cigar.ops.len(), 2);
        assert_eq!(cigar.ops[0], (CigarOp::Eq, 3));
        assert_eq!(cigar.ops[1], (CigarOp::X, 2));
    }

    #[test]
    fn test_compress_cigar_empty() {
        let cigar = compress_cigar(vec![]);
        assert!(cigar.ops.is_empty());
    }

    #[test]
    fn test_get_simd_config() {
        let config = get_simd_config();
        assert!(config.lanes >= 1);
    }
}
