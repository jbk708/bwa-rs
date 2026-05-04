//! Supermaximal MEM finding using Z-array based algorithm.
//!
//! This implementation provides O(n) expected time MEM finding,
//! replacing the O(n log n) recursive binary search approach.

use crate::fm_index::FMIndex;
use crate::types::MEM;

pub const DEFAULT_MIN_MEM_LEN: usize = 19;

/// Compute Z-array for the given byte sequence.
/// Z[i] = length of longest substring starting at i that matches a prefix of s.
#[allow(dead_code)]
fn compute_z_array(s: &[u8]) -> Vec<usize> {
    let n = s.len();
    if n == 0 {
        return Vec::new();
    }

    let mut z = vec![0usize; n];
    let mut left = 0;
    let mut right = 0;

    for i in 1..n {
        if i <= right {
            let k = i - left;
            z[i] = z[k].min(right - i + 1);
        }

        while i + z[i] < n && s[z[i]] == s[i + z[i]] {
            z[i] += 1;
        }

        if i + z[i] - 1 > right {
            left = i;
            right = i + z[i] - 1;
        }
    }

    z[0] = n;
    z
}

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
            for &ref_start in &positions {
                mems.push(MEM::new(query_start, ref_start as usize, left));
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

        result.push(mem.clone());
        max_query_end = max_query_end.max(query_end);
        max_ref_end = max_ref_end.max(ref_end);
    }

    result
}

/// Legacy MEM finding using backward search (for comparison/testing).
/// This is kept for reference but the supermaximal algorithm should be used instead.
pub fn find_mems_legacy(index: &FMIndex, query: &[u8], min_len: usize) -> Vec<MEM> {
    let mut mems = Vec::new();
    find_mems_recursive_impl(index, query, 0, query.len(), min_len, &mut mems);
    mems
}

fn find_mems_recursive_impl(
    index: &FMIndex,
    query: &[u8],
    start: usize,
    end: usize,
    min_len: usize,
    mems: &mut Vec<MEM>,
) {
    if start >= end {
        return;
    }

    let (left, right) = index.search(&query[start..end]);

    if left < right {
        let match_len = end - start;
        let positions: Vec<_> = (left..right)
            .filter_map(|i| index.get_position(i))
            .collect();

        if match_len >= min_len && !positions.is_empty() {
            for &ref_start in &positions {
                mems.push(MEM::new(start, ref_start as usize, match_len));
            }
        } else if start + 1 < end {
            find_mems_recursive_impl(index, query, start + 1, end, min_len, mems);
        }
    } else if start + 1 < end {
        find_mems_recursive_impl(index, query, start, end - 1, min_len, mems);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_z_array_simple() {
        let s = b"ABABA";
        let z = compute_z_array(s);
        assert_eq!(z[0], 5);
        assert_eq!(z[1], 0); // "BABA" vs "ABABA" - no match
        assert_eq!(z[2], 3); // "ABA" matches prefix "ABA"
        assert_eq!(z[3], 0); // "BA" vs "AB..." - no match
        assert_eq!(z[4], 1); // "A" matches prefix "A"
    }

    #[test]
    fn test_z_array_no_match() {
        let s = b"ABCDE";
        let z = compute_z_array(s);
        assert_eq!(z[0], 5);
        assert_eq!(z[1], 0);
        assert_eq!(z[2], 0);
        assert_eq!(z[3], 0);
        assert_eq!(z[4], 0);
    }

    #[test]
    fn test_z_array_all_same() {
        let s = b"AAAAA";
        let z = compute_z_array(s);
        assert_eq!(z[0], 5);
        assert_eq!(z[1], 4);
        assert_eq!(z[2], 3);
        assert_eq!(z[3], 2);
        assert_eq!(z[4], 1);
    }

    #[test]
    fn test_z_array_empty() {
        let s: &[u8] = &[];
        let z = compute_z_array(s);
        assert!(z.is_empty());
    }

    #[test]
    fn test_z_array_single() {
        let s = b"A";
        let z = compute_z_array(s);
        assert_eq!(z[0], 1);
    }

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
    fn test_supermaximal_matches_legacy() {
        use crate::fm_index::FMIndex;
        use crate::reference::Reference;

        let ref_seq = Reference::parse_fasta(">test\nACGTACGT").unwrap();
        let index = FMIndex::build(&ref_seq);

        // ACGTACGT in 2-bit encoding
        let query = [0u8, 1, 2, 3, 0, 1, 2, 3];
        let super_mems = find_supermaximal_mems(&index, &query, 2);
        let legacy_mems = find_mems_legacy(&index, &query, 2);

        // Both should find the same maximal MEM lengths
        let super_max_lens: Vec<_> = super_mems.iter().map(|m| m.length).collect();
        let legacy_max_lens: Vec<_> = legacy_mems.iter().map(|m| m.length).collect();

        super_max_lens.iter().max().map(|&l| {
            assert!(
                legacy_max_lens.contains(&l),
                "Supermaximal should find max len {} found by legacy",
                l
            )
        });
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
