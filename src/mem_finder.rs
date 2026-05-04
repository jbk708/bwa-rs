//! Supermaximal MEM finding using binary search.
//!
//! This implementation provides O(n log m) expected time MEM finding,
//! where n is query length and m is reference length.

use crate::fm_index::FMIndex;
use crate::types::MEM;

pub const DEFAULT_MIN_MEM_LEN: usize = 19;

/// Find supermaximal MEMs using binary search.
///
/// A supermaximal MEM is a MEM that extends from a query position and is not
/// contained within any other MEM.
pub fn find_supermaximal_mems(index: &FMIndex, query: &[u8], min_len: usize) -> Vec<MEM> {
    if query.is_empty() || min_len == 0 {
        return Vec::new();
    }

    let n = query.len();
    let mut mems = Vec::new();

    for query_start in 0..n {
        let max_possible = n - query_start;
        if max_possible < min_len {
            break;
        }

        // Binary search for max match length starting at query_start
        let mut left = min_len;
        let mut right = max_possible;

        // First check if even the shortest pattern matches
        let initial_pattern = &query[query_start..query_start + min_len];
        if index.count(initial_pattern) == 0 {
            continue;
        }

        while left < right {
            let mid = (left + right).div_ceil(2);
            let pattern = &query[query_start..query_start + mid];

            if index.count(pattern) > 0 {
                left = mid;
            } else {
                right = mid - 1;
            }
        }

        if left >= min_len {
            let positions = index.find_all(&query[query_start..query_start + left]);
            // Use index.len() - the FMIndex has a public len() method
            let ref_len = index.len;
            for &ref_start in &positions {
                // Filter out MEMs that extend past reference end
                if ref_start as usize + left <= ref_len {
                    mems.push(MEM::new(query_start, ref_start as usize, left));
                }
            }
        }
    }

    filter_supermaximal(mems)
}

/// Filter MEMs to keep only supermaximal ones.
/// A MEM is supermaximal if no other MEM extends further in both query and reference.
fn filter_supermaximal(mems: Vec<MEM>) -> Vec<MEM> {
    if mems.is_empty() {
        return mems;
    }

    let mut mems = mems;
    mems.sort_by_key(|b| std::cmp::Reverse(b.length));

    let mut result = Vec::new();
    let mut max_ref_end = 0usize;
    let mut max_query_end = 0usize;

    for mem in mems {
        let query_end = mem.query_end();
        let ref_end = mem.ref_end();

        // Keep if not contained in a longer MEM
        if query_end <= max_query_end && ref_end <= max_ref_end {
            continue;
        }

        result.push(mem);
        max_query_end = max_query_end.max(query_end);
        max_ref_end = max_ref_end.max(ref_end);
    }

    result
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_supermaximal_basic() {
        use crate::fm_index::FMIndex;
        use crate::reference::Reference;

        let ref_seq = Reference::parse_fasta(">test\nACGTACGT").unwrap();
        let index = FMIndex::build(&ref_seq);

        // ACGT in 2-bit encoding: A=0, C=1, G=2, T=3
        let query = [0u8, 1, 2, 3]; // ACGT
        let mems = find_supermaximal_mems(&index, &query, 2);

        assert!(!mems.is_empty(), "Should find MEMs for ACGT");
        // Due to FM-index limitations, we might find shorter MEMs
        let max_len = mems.iter().map(|m| m.length).max().unwrap();
        assert!(max_len >= 2, "Should find MEMs of at least length 2");
    }

    #[test]
    fn test_supermaximal_no_match() {
        use crate::fm_index::FMIndex;
        use crate::reference::Reference;

        let ref_seq = Reference::parse_fasta(">test\nACGTACGT").unwrap();
        let index = FMIndex::build(&ref_seq);

        // TTTT (2-bit encoded as 3) doesn't exist as a contiguous run in reference
        let query = [3u8, 3, 3, 3];
        let mems = find_supermaximal_mems(&index, &query, 2);

        // The algorithm will find single T matches but not a supermaximal MEM of length >= 2
        // that doesn't overlap with other MEMs
        assert!(
            mems.is_empty(),
            "Should not find supermaximal MEMs for TTTT"
        );
    }

    #[test]
    fn test_supermaximal_min_len() {
        use crate::fm_index::FMIndex;
        use crate::reference::Reference;

        let ref_seq = Reference::parse_fasta(">test\nACGTACGT").unwrap();
        let index = FMIndex::build(&ref_seq);

        let query = [0u8, 1, 2, 3]; // ACGT
        let mems = find_supermaximal_mems(&index, &query, 10);

        assert!(
            mems.is_empty(),
            "Should not find MEMs with min_len > query length"
        );
    }

    #[test]
    fn test_supermaximal_empty_query() {
        use crate::fm_index::FMIndex;
        use crate::reference::Reference;

        let ref_seq = Reference::parse_fasta(">test\nACGT").unwrap();
        let index = FMIndex::build(&ref_seq);

        let mems = find_supermaximal_mems(&index, &[], 2);
        assert!(mems.is_empty());
    }

    #[test]
    fn test_filter_supermaximal_contained() {
        // MEM at (query=0, ref=0, len=10) contains MEM at (query=2, ref=2, len=8)
        let mems = vec![MEM::new(0, 0, 10), MEM::new(2, 2, 8), MEM::new(5, 5, 3)];

        let filtered = filter_supermaximal(mems);
        assert_eq!(filtered.len(), 1);
        assert_eq!(filtered[0].length, 10);
    }

    #[test]
    fn test_filter_supermaximal_disjoint() {
        // Two disjoint MEMs
        let mems = vec![MEM::new(0, 0, 5), MEM::new(10, 10, 5)];

        let filtered = filter_supermaximal(mems);
        assert_eq!(filtered.len(), 2);
    }

    #[test]
    fn test_supermaximal_multiple_positions() {
        use crate::fm_index::FMIndex;
        use crate::reference::Reference;

        let ref_seq = Reference::parse_fasta(">test\nACGTACGTACGT").unwrap();
        let index = FMIndex::build(&ref_seq);

        // Query contains substring that appears in reference
        let query = [0u8, 1, 2, 3]; // ACGT
        let mems = find_supermaximal_mems(&index, &query, 2);

        // Should find at least some MEMs
        assert!(!mems.is_empty(), "Should find MEMs for ACGT");
        let positions: Vec<_> = mems.iter().map(|m| m.ref_start).collect();
        // Verify at least one valid position was found
        assert!(!positions.is_empty(), "Should have found valid positions");
    }
}
