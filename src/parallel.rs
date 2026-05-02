//! Multi-threaded alignment using Rayon.

use rayon::prelude::*;

use crate::alignment::{Aligner, Scoring};
use crate::error::BwaError;
use crate::fm_index::FMIndex;
use crate::types::AlignmentResult;

#[derive(Clone)]
pub struct ParallelAligner {
    inner: Aligner,
}

impl ParallelAligner {
    pub fn new(index: FMIndex, reference: Vec<u8>) -> Self {
        Self {
            inner: Aligner::new(index, reference),
        }
    }

    pub fn scoring(mut self, scoring: Scoring) -> Self {
        self.inner = self.inner.scoring(scoring);
        self
    }

    pub fn min_seed_len(mut self, len: usize) -> Self {
        self.inner = self.inner.min_seed_len(len);
        self
    }

    pub fn align_batch(&self, queries: &[&[u8]]) -> Vec<Result<AlignmentResult, BwaError>> {
        queries
            .par_iter()
            .map(|q| self.inner.align_read(q, None))
            .collect()
    }

    pub fn align_batch_with_mates(
        &self,
        reads: &[(Vec<u8>, Option<Vec<u8>>)],
    ) -> Vec<Result<AlignmentResult, BwaError>> {
        reads
            .par_iter()
            .map(|(read, mate)| self.inner.align_read(read, mate.as_deref()))
            .collect()
    }

    pub fn align_single(&self, query: &[u8]) -> Result<AlignmentResult, BwaError> {
        self.inner.align_read(query, None)
    }

    pub fn align_paired(
        &self,
        read1: &[u8],
        read2: &[u8],
    ) -> Result<(AlignmentResult, AlignmentResult), BwaError> {
        self.inner.align_paired(read1, read2)
    }

    pub fn align_batch_paired(
        &self,
        pairs: &[(Vec<u8>, Vec<u8>)],
    ) -> Vec<Result<(AlignmentResult, AlignmentResult), BwaError>> {
        pairs
            .par_iter()
            .map(|(r1, r2)| self.inner.align_paired(r1, r2))
            .collect()
    }
}

#[derive(Default)]
pub struct ThreadPoolConfig {
    pub num_threads: Option<usize>,
}

impl ThreadPoolConfig {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn num_threads(mut self, n: usize) -> Self {
        self.num_threads = Some(n);
        self
    }

    pub fn apply(&self) {
        if let Some(n) = self.num_threads {
            rayon::ThreadPoolBuilder::new()
                .num_threads(n)
                .build_global()
                .ok();
        }
    }
}

pub fn default_thread_count() -> usize {
    rayon::current_num_threads()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::fm_index::FMIndex;
    use crate::reference::Reference;

    fn create_test_aligner() -> (ParallelAligner, Vec<u8>) {
        let ref_seq = Reference::parse_fasta(">test\nACGTACGTACGT").unwrap();
        let ref_data = ref_seq.as_slice().to_vec();
        let index = FMIndex::build(&ref_seq);
        (
            ParallelAligner::new(index, ref_data.clone()).min_seed_len(2),
            ref_data,
        )
    }

    #[test]
    fn test_parallel_aligner_creation() {
        let (aligner, _ref_data) = create_test_aligner();
        assert_eq!(aligner.align_single(&[0, 1, 2, 3]).is_ok(), true);
    }

    #[test]
    fn test_align_batch_empty() {
        let (aligner, _ref_data) = create_test_aligner();
        let queries: Vec<&[u8]> = vec![];
        let results = aligner.align_batch(&queries);
        assert!(results.is_empty());
    }

    #[test]
    fn test_align_batch_single() {
        let (aligner, _ref_data) = create_test_aligner();
        let queries = vec![[0u8, 1, 2, 3].as_slice()];
        let results = aligner.align_batch(&queries);
        assert_eq!(results.len(), 1);
        assert!(results[0].is_ok());
    }

    #[test]
    fn test_align_batch_multiple() {
        let (aligner, _ref_data) = create_test_aligner();
        let q1 = vec![0, 1, 2, 3];
        let q2 = vec![2, 3, 0, 1];
        let q3 = vec![1, 2, 3, 0];
        let queries = vec![q1.as_slice(), q2.as_slice(), q3.as_slice()];
        let results = aligner.align_batch(&queries);
        assert_eq!(results.len(), 3);
        for r in &results {
            assert!(r.is_ok(), "Alignment should succeed: {:?}", r);
        }
    }

    #[test]
    fn test_align_batch_with_mates() {
        let (aligner, _ref_data) = create_test_aligner();
        let reads = vec![
            (vec![0, 1, 2, 3], Some(vec![3, 2, 1, 0])),
            (vec![1, 2, 3, 0], None),
        ];
        let results = aligner.align_batch_with_mates(&reads);
        assert_eq!(results.len(), 2);
        assert!(results[0].is_ok());
        assert!(results[1].is_ok());
    }

    #[test]
    fn test_align_batch_paired() {
        let (aligner, _ref_data) = create_test_aligner();
        let pairs = vec![
            (vec![0, 1, 2, 3], vec![3, 2, 1, 0]),
            (vec![1, 2, 3, 0], vec![0, 1, 2, 3]),
        ];
        let results = aligner.align_batch_paired(&pairs);
        assert_eq!(results.len(), 2);
        for r in &results {
            assert!(r.is_ok(), "Paired alignment should succeed: {:?}", r);
        }
    }

    #[test]
    fn test_align_batch_parallelism() {
        let (aligner, _ref_data) = create_test_aligner();
        let queries: Vec<Vec<u8>> = (0..100)
            .map(|i| vec![i as u8 % 4, (i + 1) as u8 % 4, (i + 2) as u8 % 4])
            .collect();
        let query_refs: Vec<&[u8]> = queries.iter().map(|v| v.as_slice()).collect();
        let results = aligner.align_batch(&query_refs);
        assert_eq!(results.len(), 100);
        let success_count = results.iter().filter(|r| r.is_ok()).count();
        assert!(success_count > 0, "At least some alignments should succeed");
    }

    #[test]
    fn test_align_batch_consistency() {
        let (aligner, _ref_data) = create_test_aligner();
        let query = vec![0, 1, 2, 3];

        let single_result = aligner.align_single(&query);
        let batch_result = aligner.align_batch(&[&query]);

        assert!(single_result.is_ok());
        assert!(batch_result[0].is_ok());

        let single = single_result.unwrap();
        let batch = batch_result[0].as_ref().unwrap();

        assert_eq!(single.position, batch.position);
        assert_eq!(single.mapq, batch.mapq);
    }

    #[test]
    fn test_thread_pool_config_default() {
        let config = ThreadPoolConfig::new();
        assert!(config.num_threads.is_none());
    }

    #[test]
    fn test_thread_pool_config_with_threads() {
        let config = ThreadPoolConfig::new().num_threads(4);
        assert_eq!(config.num_threads, Some(4));
    }

    #[test]
    fn test_default_thread_count() {
        let count = default_thread_count();
        assert!(count >= 1);
    }

    #[test]
    fn test_scoring_builder() {
        let ref_seq = Reference::parse_fasta(">test\nACGT").unwrap();
        let ref_data = ref_seq.as_slice().to_vec();
        let index = FMIndex::build(&ref_seq);

        let scoring = Scoring {
            match_score: 2,
            mismatch_penalty: 5,
            gap_open: 8,
            gap_extend: 2,
        };

        let aligner = ParallelAligner::new(index, ref_data).scoring(scoring.clone());

        let result = aligner.align_single(&[0, 1, 2, 3]);
        assert!(result.is_ok());
    }

    #[test]
    fn test_min_seed_len_builder() {
        let ref_seq = Reference::parse_fasta(">test\nACGTACGT").unwrap();
        let ref_data = ref_seq.as_slice().to_vec();
        let index = FMIndex::build(&ref_seq);

        let aligner = ParallelAligner::new(index, ref_data).min_seed_len(3);
        let result = aligner.align_single(&[0, 1, 2, 3]);
        assert!(result.is_ok());
    }

    #[test]
    fn test_align_unmapped_in_batch() {
        let ref_seq = Reference::parse_fasta(">test\nACGT").unwrap();
        let ref_data = ref_seq.as_slice().to_vec();
        let index = FMIndex::build(&ref_seq);

        let aligner = ParallelAligner::new(index, ref_data).min_seed_len(100);

        let q = vec![0, 1, 2, 3, 0, 1, 2, 3];
        let queries = vec![q.as_slice()];
        let results = aligner.align_batch(&queries);
        assert_eq!(results.len(), 1);

        let result = results[0].as_ref().unwrap();
        assert_eq!(result.flag & 0x4, 0x4);
    }
}
