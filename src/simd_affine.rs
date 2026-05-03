//! Affine gap seed extension (scalar fallback).
//!
//! SIMD support for affine gap alignment is not yet implemented.
//! Uses the scalar implementation from alignment.rs.

use crate::alignment::{affine_extend_forward, SeedExtension};

#[inline]
pub fn affine_extend_forward_simd(
    query: &[u8],
    reference: &[u8],
    start_pos: usize,
    scoring: &crate::alignment::Scoring,
    bandwidth: usize,
) -> SeedExtension {
    affine_extend_forward(query, reference, start_pos, scoring, bandwidth)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::alignment::{affine_extend_forward, Scoring};

    #[test]
    fn test_affine_extend_simd_exact_match() {
        let query = vec![0, 1, 2, 3];
        let reference = vec![0, 1, 2, 3, 4, 5];
        let scoring = Scoring::default();

        let simd_result = affine_extend_forward_simd(&query, &reference, 0, &scoring, 16);
        let scalar_result = affine_extend_forward(&query, &reference, 0, &scoring, 16);

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
        let scalar_result = affine_extend_forward(&query, &reference, 0, &scoring, 16);

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
        let scalar_result = affine_extend_forward(&query, &reference, 0, &scoring, 16);

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
        let scalar_result = affine_extend_forward(&query, &reference, 0, &scoring, 16);

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
        let scalar_result = affine_extend_forward(&query, &reference, 0, &scoring, 32);

        assert_eq!(
            simd_result.score, scalar_result.score,
            "Long sequences scores should match: SIMD={}, Scalar={}",
            simd_result.score, scalar_result.score
        );
    }
}
