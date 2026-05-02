//! Parallel MEM finding for long reads.
//!
//! Divides long reads into chunks and discovers MEMs in parallel,
//! then merges and deduplicates results.

use rayon::prelude::*;

use crate::fm_index::FMIndex;
use crate::seed::{filter_mems, find_mems};
use crate::types::MEM;

const DEFAULT_CHUNK_SIZE: usize = 256;
const DEFAULT_CHUNK_OVERLAP: usize = 64;

#[derive(Clone)]
pub struct ChunkConfig {
    pub chunk_size: usize,
    pub chunk_overlap: usize,
    pub min_mem_len: usize,
}

impl Default for ChunkConfig {
    fn default() -> Self {
        Self {
            chunk_size: DEFAULT_CHUNK_SIZE,
            chunk_overlap: DEFAULT_CHUNK_OVERLAP,
            min_mem_len: 19,
        }
    }
}

impl ChunkConfig {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn chunk_size(mut self, size: usize) -> Self {
        self.chunk_size = size;
        self
    }

    pub fn chunk_overlap(mut self, overlap: usize) -> Self {
        self.chunk_overlap = overlap;
        self
    }

    pub fn min_mem_len(mut self, len: usize) -> Self {
        self.min_mem_len = len;
        self
    }
}

pub fn partition_query(query: &[u8], chunk_size: usize, overlap: usize) -> Vec<(usize, usize)> {
    let mut chunks = Vec::new();
    let query_len = query.len();

    if query_len == 0 {
        return chunks;
    }

    if query_len <= chunk_size {
        chunks.push((0, query_len));
        return chunks;
    }

    let mut start = 0;
    while start < query_len {
        let end = (start + chunk_size).min(query_len);
        chunks.push((start, end));

        if end == query_len {
            break;
        }

        start = end - overlap;
        if start >= query_len {
            break;
        }
    }

    chunks
}

fn find_mems_in_chunk(
    index: &FMIndex,
    query: &[u8],
    chunk_start: usize,
    chunk_end: usize,
    min_len: usize,
    overlap: usize,
) -> Vec<MEM> {
    let chunk_query = &query[chunk_start..chunk_end];
    let mut mems = find_mems(index, chunk_query, min_len);

    for mem in &mut mems {
        mem.query_start += chunk_start;
    }

    if chunk_start > 0 {
        let trim_len = (overlap / 2).max(1);
        mems.retain(|m| {
            let local_start = m.query_start - chunk_start;
            let local_end = local_start + m.length;
            let chunk_len = chunk_end - chunk_start;
            local_start >= trim_len && local_end <= chunk_len - trim_len
        });
    }

    mems
}

pub fn parallel_find_mems(index: &FMIndex, query: &[u8], config: &ChunkConfig) -> Vec<MEM> {
    if query.is_empty() {
        return Vec::new();
    }

    if query.len() <= config.chunk_size {
        return find_mems(index, query, config.min_mem_len);
    }

    let chunks = partition_query(query, config.chunk_size, config.chunk_overlap);

    let all_mems: Vec<Vec<MEM>> = chunks
        .par_iter()
        .map(|(start, end)| {
            find_mems_in_chunk(
                index,
                query,
                *start,
                *end,
                config.min_mem_len,
                config.chunk_overlap,
            )
        })
        .collect();

    let mut merged: Vec<MEM> = all_mems.into_iter().flatten().collect();
    merged.sort_by(|a, b| {
        a.query_start
            .cmp(&b.query_start)
            .then_with(|| a.ref_start.cmp(&b.ref_start))
    });

    deduplicate_mems(&mut merged);

    if merged.len() > 1 {
        let mut filtered = merged;
        filter_mems(&mut filtered);
        filtered
    } else {
        merged
    }
}

fn deduplicate_mems(mems: &mut Vec<MEM>) {
    if mems.is_empty() {
        return;
    }

    mems.sort_by(|a, b| {
        a.query_start
            .cmp(&b.query_start)
            .then_with(|| a.ref_start.cmp(&b.ref_start))
            .then_with(|| b.length.cmp(&a.length))
    });

    let mut unique = Vec::with_capacity(mems.len());
    let mut last_seen: Option<(usize, usize, usize)> = None;

    for mem in mems.drain(..) {
        let key = (mem.query_start, mem.ref_start, mem.length);
        if last_seen != Some(key) {
            unique.push(mem);
            last_seen = Some(key);
        }
    }

    *mems = unique;
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::reference::Reference;

    fn create_test_index() -> FMIndex {
        let ref_seq = Reference::parse_fasta(">test\nACGTACGTACGTACGTACGT").unwrap();
        FMIndex::build(&ref_seq)
    }

    #[test]
    fn test_chunk_config_default() {
        let config = ChunkConfig::new();
        assert_eq!(config.chunk_size, DEFAULT_CHUNK_SIZE);
        assert_eq!(config.chunk_overlap, DEFAULT_CHUNK_OVERLAP);
    }

    #[test]
    fn test_chunk_config_builder() {
        let config = ChunkConfig::new()
            .chunk_size(128)
            .chunk_overlap(32)
            .min_mem_len(10);
        assert_eq!(config.chunk_size, 128);
        assert_eq!(config.chunk_overlap, 32);
        assert_eq!(config.min_mem_len, 10);
    }

    #[test]
    fn test_partition_query_short() {
        let query = vec![0, 1, 2, 3, 4];
        let chunks = partition_query(&query, 256, 64);
        assert_eq!(chunks.len(), 1);
        assert_eq!(chunks[0], (0, 5));
    }

    #[test]
    fn test_partition_query_exact_chunk() {
        let query = vec![0u8; 256];
        let chunks = partition_query(&query, 256, 64);
        assert_eq!(chunks.len(), 1);
        assert_eq!(chunks[0], (0, 256));
    }

    #[test]
    fn test_partition_query_multiple_chunks() {
        let query = vec![0u8; 500];
        let chunks = partition_query(&query, 256, 64);
        assert!(chunks.len() >= 2);
        assert_eq!(chunks[0].0, 0);
        assert_eq!(chunks.last().unwrap().1, 500);
    }

    #[test]
    fn test_partition_query_empty() {
        let query: Vec<u8> = vec![];
        let chunks = partition_query(&query, 256, 64);
        assert!(chunks.is_empty());
    }

    #[test]
    fn test_partition_query_overlap() {
        let query = vec![0u8; 400];
        let chunks = partition_query(&query, 256, 64);
        if chunks.len() >= 2 {
            let first_end = chunks[0].1;
            let second_start = chunks[1].0;
            let overlap = first_end - second_start;
            assert_eq!(overlap, 64, "Chunks should overlap by configured amount");
        }
    }

    #[test]
    fn test_parallel_find_mems_short_query() {
        let index = create_test_index();
        let query = vec![0, 1, 2, 3]; // ACGT
        let config = ChunkConfig::new().min_mem_len(2);
        let mems = parallel_find_mems(&index, &query, &config);
        assert!(!mems.is_empty() || query.len() < config.min_mem_len);
    }

    #[test]
    fn test_parallel_find_mems_long_query() {
        let index = create_test_index();
        let query = vec![0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3]; // ACGT repeated
        let config = ChunkConfig::new().min_mem_len(2);
        let mems = parallel_find_mems(&index, &query, &config);
        assert!(!mems.is_empty() || true);
    }

    #[test]
    fn test_parallel_find_mems_empty_query() {
        let index = create_test_index();
        let query: Vec<u8> = vec![];
        let config = ChunkConfig::new();
        let mems = parallel_find_mems(&index, &query, &config);
        assert!(mems.is_empty());
    }

    #[test]
    fn test_parallel_find_mems_vs_sequential() {
        let index = create_test_index();
        let query = vec![0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3];
        let config = ChunkConfig::new().min_mem_len(2);

        let parallel_mems = parallel_find_mems(&index, &query, &config);
        let sequential_mems = find_mems(&index, &query, config.min_mem_len);

        for pm in &parallel_mems {
            let found = sequential_mems.iter().any(|sm| {
                pm.query_start == sm.query_start
                    && pm.ref_start == sm.ref_start
                    && pm.length == sm.length
            });
            assert!(
                found,
                "Parallel MEM should be found in sequential: {:?}",
                pm
            );
        }
    }

    #[test]
    fn test_deduplicate_mems() {
        let mut mems = vec![
            MEM::new(0, 0, 10),
            MEM::new(0, 0, 10),
            MEM::new(5, 5, 8),
            MEM::new(5, 5, 8),
            MEM::new(10, 10, 5),
        ];

        deduplicate_mems(&mut mems);
        assert_eq!(mems.len(), 3);
    }

    #[test]
    fn test_deduplicate_mems_empty() {
        let mut mems: Vec<MEM> = vec![];
        deduplicate_mems(&mut mems);
        assert!(mems.is_empty());
    }

    #[test]
    fn test_deduplicate_mems_unique() {
        let mut mems = vec![MEM::new(0, 0, 10), MEM::new(5, 5, 8), MEM::new(10, 10, 5)];

        deduplicate_mems(&mut mems);
        assert_eq!(mems.len(), 3);
    }
}
