//! MEM finding using FM-index.
//!
//! This module provides MEM (Maximal Exact Match) finding functionality.
//! The primary algorithm uses supermaximal MEM finding for O(n) expected time.

use crate::fm_index::FMIndex;
use crate::mem_finder::find_supermaximal_mems;
use crate::types::MEM;

pub const DEFAULT_MIN_SEED_LEN: usize = 19;

/// Find MEMs using supermaximal algorithm (O(n) expected time).
///
/// `max_occ` caps the number of reference occurrences a SMEM may have; SMEMs
/// with more hits than `max_occ` are discarded. Use `DEFAULT_MAX_OCC` for the
/// bwa-mem2-equivalent behaviour.
pub fn find_mems(index: &FMIndex, query: &[u8], min_len: usize, max_occ: usize) -> Vec<MEM> {
    find_supermaximal_mems(index, query, min_len, max_occ).0
}

/// Like [`find_mems`] but also returns bwa's `frac_rep` (the fraction of the read
/// covered by over-frequent repetitive SMEMs), used to down-weight MAPQ.
pub fn find_mems_with_frac_rep(
    index: &FMIndex,
    query: &[u8],
    min_len: usize,
    max_occ: usize,
) -> (Vec<MEM>, f32) {
    find_supermaximal_mems(index, query, min_len, max_occ)
}

pub fn filter_mems(mems: &mut Vec<MEM>) {
    mems.sort_by(|a, b| {
        b.length
            .cmp(&a.length)
            .then_with(|| a.ref_start.cmp(&b.ref_start))
    });

    let mut filtered = Vec::new();
    let mut used = Vec::new();

    for mem in mems.drain(..) {
        let overlaps = used
            .iter()
            .any(|u: &MEM| mem.ref_start < u.ref_end() && mem.ref_end() > u.ref_start);

        if !overlaps {
            filtered.push(mem);
            used.push(mem);
        }
    }

    *mems = filtered;
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_filter_overlapping() {
        let mut mems = vec![MEM::new(0, 0, 10), MEM::new(5, 5, 10), MEM::new(20, 20, 8)];

        filter_mems(&mut mems);
        assert_eq!(mems.len(), 2);
    }

    #[test]
    fn test_mem_positions() {
        let mem = MEM::new(5, 100, 20);
        assert_eq!(mem.query_end(), 25);
        assert_eq!(mem.ref_end(), 120);
    }

    #[test]
    fn test_default_min_seed_len() {
        assert_eq!(DEFAULT_MIN_SEED_LEN, 19);
    }
}
