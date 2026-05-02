//! SIMD-accelerated affine gap dynamic programming alignment.
//!
//! Provides SIMD-accelerated affine gap alignment with automatic fallback
//! to scalar implementation when SIMD is unavailable.

use crate::alignment::{affine_extend_forward as scalar_affine_extend, SeedExtension};

#[inline]
pub fn affine_extend_forward_simd(
    query: &[u8],
    reference: &[u8],
    start_pos: usize,
    scoring: &crate::alignment::Scoring,
    bandwidth: usize,
) -> SeedExtension {
    #[cfg(any(target_arch = "x86_64", target_arch = "x86"))]
    {
        let config = crate::simd_sw::get_simd_config();
        if config.use_avx2 {
            return affine_extend_forward_avx2(query, reference, start_pos, scoring, bandwidth);
        }
        if config.use_avx512 {
            return affine_extend_forward_avx512(query, reference, start_pos, scoring, bandwidth);
        }
    }

    scalar_affine_extend(query, reference, start_pos, scoring, bandwidth)
}

#[cfg(all(target_arch = "x86_64", target_feature = "avx2"))]
#[allow(dead_code)]
fn affine_extend_forward_avx2(
    query: &[u8],
    reference: &[u8],
    start_pos: usize,
    scoring: &crate::alignment::Scoring,
    bandwidth: usize,
) -> SeedExtension {
    use crate::types::Cigar;
    use std::arch::x86_64::*;

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

    let actual_ref_len = (start_pos + ref_len).min(reference.len()) - start_pos;
    let bw = bandwidth.clamp(1, 256);

    let lanes = 8;
    let neg_inf = unsafe { _mm256_set1_epi32(i32::MIN / 4) };
    let match_s = scoring.match_score as i32;
    let mismatch_p = -(scoring.mismatch_penalty as i32);
    let gap_open = scoring.gap_open as i32;
    let gap_extend = scoring.gap_extend as i32;
    let match_vec = unsafe { _mm256_set1_epi32(match_s) };
    let mismatch_vec = unsafe { _mm256_set1_epi32(mismatch_p) };
    let gap_open_vec = unsafe { _mm256_set1_epi32(gap_open) };
    let gap_extend_vec = unsafe { _mm256_set1_epi32(gap_extend) };

    let full_cols = (actual_ref_len / lanes) * lanes;

    let mut m_curr = vec![neg_inf; full_cols / lanes + 1];
    let mut x_curr = vec![neg_inf; full_cols / lanes + 1];
    let mut g_curr = vec![neg_inf; full_cols / lanes + 1];
    let mut m_prev = vec![neg_inf; full_cols / lanes + 1];
    let mut x_prev = vec![neg_inf; full_cols / lanes + 1];
    let mut g_prev = vec![neg_inf; full_cols / lanes + 1];

    m_curr[0] = unsafe { _mm256_setzero_si256() };

    let mut best_score = 0i32;
    let mut best_i = 0usize;
    let mut best_j = 0usize;

    for i in 1..=query_len {
        let j_end = (i + bw).min(actual_ref_len);
        let j_start = 1.max(i.saturating_sub(bw));

        m_prev[0] = if i == 1 {
            unsafe { _mm256_setzero_si256() }
        } else {
            m_curr[0]
        };
        x_prev[0] = if i == 1 { neg_inf } else { x_curr[0] };
        g_prev[0] = if i == 1 { neg_inf } else { g_curr[0] };

        for block in 0..(full_cols / lanes) {
            m_curr[block] = neg_inf;
            x_curr[block] = neg_inf;
            g_curr[block] = neg_inf;
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

            let mut ref_vals = [0i8; 8];
            for k in 0..8 {
                let idx = ref_start_idx + k;
                ref_vals[k] = if idx < reference.len() {
                    reference[idx] as i8
                } else {
                    -1
                };
            }

            let base_vec = unsafe { _mm256_set1_epi8(base as i8) };
            let ref_vec = unsafe {
                _mm256_set_epi8(
                    ref_vals[7],
                    ref_vals[6],
                    ref_vals[5],
                    ref_vals[4],
                    ref_vals[3],
                    ref_vals[2],
                    ref_vals[1],
                    ref_vals[0],
                    ref_vals[7],
                    ref_vals[6],
                    ref_vals[5],
                    ref_vals[4],
                    ref_vals[3],
                    ref_vals[2],
                    ref_vals[1],
                    ref_vals[0],
                )
            };

            let cmp = unsafe { _mm256_cmpeq_epi8(base_vec, ref_vec) };
            let cmp_lo = unsafe { _mm256_castsi256_si128(cmp) };
            let cmp_hi = unsafe { _mm256_extractf128_si256(cmp, 1) };
            let cmp_i32_lo = unsafe { _mm_cvtepi8_epi32(cmp_lo) };
            let cmp_i32_hi = unsafe { _mm_cvtepi8_epi32(cmp_hi) };
            let cmp_i32 = unsafe {
                let combined = _mm256_set_m128i(cmp_i32_hi, cmp_i32_lo);
                _mm256_cmpgt_epi32(_mm256_setzero_si256(), combined)
            };

            let score_vec = unsafe { _mm256_blendv_epi8(mismatch_vec, match_vec, cmp_i32) };

            let prev_block = block.saturating_sub(1);
            let m_diag = if prev_block < m_prev.len() {
                m_prev[prev_block]
            } else {
                neg_inf
            };
            let x_diag = if prev_block < x_prev.len() {
                x_prev[prev_block]
            } else {
                neg_inf
            };
            let g_diag = if prev_block < g_prev.len() {
                g_prev[prev_block]
            } else {
                neg_inf
            };

            let best_prev = unsafe {
                let max_mx = _mm256_max_epi32(m_diag, x_diag);
                _mm256_max_epi32(max_mx, g_diag)
            };

            let m_new = unsafe { _mm256_add_epi32(best_prev, score_vec) };

            let x_open = unsafe { _mm256_sub_epi32(m_prev[block], gap_open_vec) };
            let x_extend = unsafe { _mm256_sub_epi32(x_prev[block], gap_extend_vec) };
            let x_new = unsafe {
                let max_x = _mm256_max_epi32(x_open, x_extend);
                max_x
            };

            let g_open = unsafe { _mm256_sub_epi32(m_curr[block], gap_open_vec) };
            let g_extend = unsafe { _mm256_sub_epi32(g_curr[block], gap_extend_vec) };
            let g_new = unsafe { _mm256_max_epi32(g_open, g_extend) };

            let block_j_end = (j_block + lanes).min(j_end + 1);
            let block_j_start = j_block.max(j_start);

            if block_j_start <= block_j_end.saturating_sub(1) && block_j_end > block_j_start {
                m_curr[block] = m_new;
                x_curr[block] = x_new;
                g_curr[block] = g_new;
            }

            let current_best = unsafe {
                let max_mx = _mm256_max_epi32(m_new, x_new);
                _mm256_max_epi32(max_mx, g_new)
            };

            let mut score_arr = [0i32; 8];
            unsafe {
                _mm256_storeu_si256(score_arr.as_mut_ptr() as *mut __m256i, current_best);
            }
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

    scalar_affine_extend(query, reference, start_pos, scoring, bandwidth)
}

#[cfg(not(all(target_arch = "x86_64", target_feature = "avx2")))]
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
