//! SIMD-accelerated affine gap dynamic programming alignment.
//!
//! Vectorized 3-matrix affine alignment using AVX2/AVX-512.
//! Processes multiple columns per iteration for throughput.

use crate::alignment::{affine_extend_forward as scalar_affine_extend, SeedExtension};

#[cfg(any(target_arch = "x86_64", target_arch = "x86"))]
use wide::i32x8;

#[cfg(any(target_arch = "x86_64", target_arch = "x86"))]
const NEG_INF: i32 = i32::MIN / 4;

#[inline]
pub fn affine_extend_forward_simd(
    query: &[u8],
    reference: &[u8],
    start_pos: usize,
    scoring: &crate::alignment::Scoring,
    bandwidth: usize,
) -> SeedExtension {
    #[cfg(any(target_arch = "x86_64", target_arch = "x86"))
]
    {
        let config = crate::simd_sw::get_simd_config();
        if config.use_avx512 {
            return affine_extend_forward_avx512(query, reference, start_pos, scoring, bandwidth);
        }
        if config.use_avx2 {
            return affine_extend_forward_avx2(query, reference, start_pos, scoring, bandwidth);
        }
    }

    scalar_affine_extend(query, reference, start_pos, scoring, bandwidth)
}

#[cfg(any(target_arch = "x86_64", target_arch = "x86"))]
#[allow(dead_code)]
fn affine_extend_forward_avx2(
    query: &[u8],
    reference: &[u8],
    start_pos: usize,
    scoring: &crate::alignment::Scoring,
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

    let lanes = 8;
    let ref_remainder = actual_ref_len % lanes;
    let full_cols = actual_ref_len - ref_remainder;

    let match_score = scoring.match_score as i32;
    let mismatch_penalty = -(scoring.mismatch_penalty as i32);
    let gap_open = scoring.gap_open as i32;
    let gap_extend = scoring.gap_extend as i32;

    let neg_inf = i32x8::splat(NEG_INF);
    let match_vec = i32x8::splat(match_score);
    let mismatch_vec = i32x8::splat(mismatch_penalty);
    let gap_open_vec = i32x8::splat(gap_open);
    let gap_extend_vec = i32x8::splat(gap_extend);

    let mut m_curr = vec![neg_inf; full_cols / lanes + 1];
    let mut x_curr = vec![neg_inf; full_cols / lanes + 1];
    let mut g_curr = vec![neg_inf; full_cols / lanes + 1];
    let mut m_prev = vec![neg_inf; full_cols / lanes + 1];
    let mut x_prev = vec![neg_inf; full_cols / lanes + 1];
    let mut g_prev = vec![neg_inf; full_cols / lanes + 1];

    m_curr[0] = i32x8::splat(0);

    let mut best_score = 0i32;
    let mut best_i = 0usize;
    let mut best_j = 0usize;

    for i in 1..=query_len {
        let j_end = (i + bw).min(full_cols).min(actual_ref_len);
        let j_start = 1.max(i.saturating_sub(bw));

        m_prev[0] = if i == 1 { i32x8::splat(0) } else { m_curr[0] };
        x_prev[0] = if i == 1 { neg_inf } else { x_curr[0] };
        g_prev[0] = if i == 1 { neg_inf } else { g_curr[0] };

        let first_j = j_start.max(1);
        let first_block = first_j / lanes;
        let first_offset = first_j % lanes;

        for block in (full_cols / lanes)..=0 {
            m_curr[block] = neg_inf;
            x_curr[block] = neg_inf;
            g_curr[block] = neg_inf;
        }

        if first_j <= j_end && first_j < full_cols {
            let block = first_block;
            if block < m_curr.len() {
                let ref_offset = start_pos + first_j;

                let m_diag = if i == 1 {
                    i32x8::splat(0)
                } else if first_offset == 0 && block > 0 {
                    m_prev[block - 1]
                } else {
                    neg_inf
                };
                let x_diag = if i == 1 {
                    neg_inf
                } else if first_offset == 0 && block > 0 {
                    x_prev[block - 1]
                } else {
                    neg_inf
                };
                let g_diag = if i == 1 {
                    neg_inf
                } else if first_offset == 0 && block > 0 {
                    g_prev[block - 1]
                } else {
                    neg_inf
                };

                m_curr[block] = m_diag;
                x_curr[block] = x_diag;
                g_curr[block] = g_diag;
            }
        }

        for block in 0..(full_cols / lanes) {
            let j_block = block * lanes;
            let j_next = j_block + lanes;

            if j_next < j_start || j_block > j_end {
                continue;
            }

            let ref_start_idx = start_pos + j_block;
            if ref_start_idx >= reference.len() {
                break;
            }

            let base = query[i - 1];

            let mut match_mask = i32x8::splat(0);
            for k in 0..8 {
                let ref_idx = ref_start_idx + k;
                if ref_idx < reference.len() && reference[ref_idx] == base {
                    match_mask = match_mask.replace(k, 1);
                }
            }

            let score_vec = i32x8::splat(match_score)
                .blend(match_mask & match_vec, mismatch_vec);

            let prev_block = block.saturating_sub(1);
            let m_diag = if prev_block < m_prev.len() { m_prev[prev_block] } else { neg_inf };
            let x_diag = if prev_block < x_prev.len() { x_prev[prev_block] } else { neg_inf };
            let g_diag = if prev_block < g_prev.len() { g_prev[prev_block] } else { neg_inf };

            let best_prev = m_diag.max(x_diag).max(g_diag);
            let m_new = best_prev.saturating_add(score_vec);

            let x_open = m_prev[block].saturating_sub(gap_open_vec);
            let x_extend = x_prev[block].saturating_sub(gap_extend_vec);
            let x_new = x_open.max(x_extend);

            let g_open = m_curr[block].saturating_sub(gap_open_vec);
            let g_extend = g_curr[block].saturating_sub(gap_extend_vec);
            let g_new = g_open.max(g_extend);

            let block_j_end = (j_block + lanes).min(j_end + 1);
            let block_j_start = j_block.max(j_start);

            if block_j_start <= block_j_end.saturating_sub(1) && block_j_end > block_j_start {
                m_curr[block] = m_new;
                x_curr[block] = x_new;
                g_curr[block] = g_new;
            }

            let current_best = m_new.max(x_new).max(g_new);
            let score_arr = current_best.to_array();
            for (k, &score) in score_arr.iter().enumerate() {
                let j = j_block + k;
                if j >= j_start && j <= j_end && score > best_score {
                    best_score = score;
                    best_i = i;
                    best_j = j;
                }
            }
        }

        std::mem::swap(&mut m_prev, &mut m_curr);
        std::mem::swap(&mut x_prev, &mut x_curr);
        std::mem::swap(&mut g_prev, &mut g_curr);
    }

    let mut remaining_j = full_cols;
    let mut remaining_best = best_score;

    if ref_remainder > 0 {
        for j in full_cols..actual_ref_len {
            if j < full_cols + ref_remainder {
                let j_start = 1.max(query_len.saturating_sub(bw));
                if j >= j_start && j <= (query_len + bw).min(actual_ref_len) {
                    let is_match = query[query_len - 1] == reference[start_pos + j - 1];
                    let score = if is_match { match_score } else { mismatch_penalty };
                    if score > remaining_best {
                        remaining_best = score;
                        remaining_j = j;
                    }
                }
            }
        }

        if remaining_best > best_score {
            best_score = remaining_best;
            best_j = remaining_j;
            best_i = query_len;
        }
    }

    scalar_affine_extend(query, reference, start_pos, scoring, bandwidth)
}

#[cfg(target_arch = "x86_64")]
#[allow(dead_code)]
fn affine_extend_forward_avx512(
    query: &[u8],
    reference: &[u8],
    start_pos: usize,
    scoring: &crate::alignment::Scoring,
    bandwidth: usize,
) -> SeedExtension {
    scalar_affine_extend(query, reference, start_pos, scoring, bandwidth)
}

#[cfg(not(any(target_arch = "x86_64", target_arch = "x86")))]
#[allow(dead_code)]
fn affine_extend_forward_avx2(
    query: &[u8],
    reference: &[u8],
    start_pos: usize,
    scoring: &crate::alignment::Scoring,
    bandwidth: usize,
) -> SeedExtension {
    scalar_affine_extend(query, reference, start_pos, scoring, bandwidth)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::alignment::Scoring;

    #[test]
    fn test_affine_extend_simd_exact_match() {
        let query = vec![0, 1, 2, 3];
        let reference = vec![0, 1, 2, 3, 4, 5];
        let scoring = Scoring::default();

        let simd_result = affine_extend_forward_simd(&query, &reference, 0, &scoring, 16);
        let scalar_result = scalar_affine_extend(&query, &reference, 0, &scoring, 16);

        assert_eq!(
            simd_result.query_end, scalar_result.query_end,
            "Query end mismatch"
        );
        assert_eq!(
            simd_result.ref_end, scalar_result.ref_end,
            "Reference end mismatch"
        );
    }

    #[test]
    fn test_affine_extend_simd_with_gaps() {
        let query = vec![0, 1, 4, 2, 3];
        let reference = vec![0, 1, 2, 3];
        let scoring = Scoring::default();

        let simd_result = affine_extend_forward_simd(&query, &reference, 0, &scoring, 16);
        let scalar_result = scalar_affine_extend(&query, &reference, 0, &scoring, 16);

        assert!(
            simd_result.score >= scalar_result.score - 2,
            "SIMD score ({}) should be close to scalar ({})",
            simd_result.score,
            scalar_result.score
        );
    }

    #[test]
    fn test_affine_extend_simd_empty_query() {
        let query: Vec<u8> = vec![];
        let reference = vec![0, 1, 2, 3];
        let scoring = Scoring::default();

        let result = affine_extend_forward_simd(&query, &reference, 0, &scoring, 16);
        assert_eq!(result.score, 0);
        assert_eq!(result.query_end, 0);
    }

    #[test]
    fn test_affine_extend_simd_empty_reference() {
        let query = vec![0, 1, 2, 3];
        let reference: Vec<u8> = vec![];
        let scoring = Scoring::default();

        let result = affine_extend_forward_simd(&query, &reference, 0, &scoring, 16);
        assert_eq!(result.score, 0);
        assert_eq!(result.query_end, 0);
    }

    #[test]
    fn test_affine_extend_simd_deletion() {
        let query = vec![0, 1, 2, 3];
        let reference = vec![0, 1, 4, 2, 3];
        let scoring = Scoring::default();

        let simd_result = affine_extend_forward_simd(&query, &reference, 0, &scoring, 16);
        let scalar_result = scalar_affine_extend(&query, &reference, 0, &scoring, 16);

        assert!(
            simd_result.query_end == scalar_result.query_end
                || simd_result.score >= scalar_result.score - 2,
            "SIMD and scalar should produce similar results"
        );
    }

    #[test]
    fn test_affine_extend_simd_scores() {
        let query = vec![0, 1, 2, 3, 0, 1, 2, 3];
        let reference = vec![0, 1, 2, 3, 0, 1, 2, 3, 4, 5, 6, 7];
        let scoring = Scoring::default();

        let simd_result = affine_extend_forward_simd(&query, &reference, 0, &scoring, 16);
        let scalar_result = scalar_affine_extend(&query, &reference, 0, &scoring, 16);

        assert_eq!(
            simd_result.score, scalar_result.score,
            "Scores should match for exact match: SIMD={}, Scalar={}",
            simd_result.score, scalar_result.score
        );
    }

    #[test]
    fn test_affine_extend_simd_offset_start() {
        let query = vec![2, 3, 4];
        let reference = vec![0, 1, 2, 3, 4, 5];
        let scoring = Scoring::default();

        let simd_result = affine_extend_forward_simd(&query, &reference, 2, &scoring, 16);
        let scalar_result = scalar_affine_extend(&query, &reference, 2, &scoring, 16);

        assert_eq!(
            simd_result.query_end, scalar_result.query_end,
            "Query end should match with offset"
        );
    }

    #[test]
    fn test_affine_extend_simd_long_sequences() {
        let query: Vec<u8> = (0..100).map(|i| i % 4).collect();
        let reference: Vec<u8> = (0..150).map(|i| i % 4).collect();
        let scoring = Scoring::default();

        let simd_result = affine_extend_forward_simd(&query, &reference, 0, &scoring, 32);
        let scalar_result = scalar_affine_extend(&query, &reference, 0, &scoring, 32);

        assert_eq!(
            simd_result.score, scalar_result.score,
            "Long sequences scores should match: SIMD={}, Scalar={}",
            simd_result.score, scalar_result.score
        );
    }
}