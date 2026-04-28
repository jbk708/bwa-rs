//! MEM finding using FM-index.

use crate::fm_index::FMIndex;
use crate::types::MEM;

pub const DEFAULT_MIN_SEED_LEN: usize = 19;

pub fn find_mems(index: &FMIndex, query: &[u8], min_len: usize) -> Vec<MEM> {
    let mut mems = Vec::new();
    find_mems_recursive(index, query, 0, query.len(), min_len, &mut mems);
    mems
}

fn find_mems_recursive(
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
        let positions = (left..right)
            .filter_map(|i| index.get_position(i))
            .collect::<Vec<_>>();

        if match_len >= min_len && !positions.is_empty() {
            for &ref_start in &positions {
                mems.push(MEM::new(start, ref_start as usize, match_len));
            }
        } else if start + 1 < end {
            find_mems_recursive(index, query, start + 1, end, min_len, mems);
        }
    } else if start + 1 < end {
        find_mems_recursive(index, query, start, end - 1, min_len, mems);
    }
}

pub fn filter_mems(mems: &mut Vec<MEM>) {
    mems.sort_by(|a, b| {
        b.length.cmp(&a.length)
            .then_with(|| a.ref_start.cmp(&b.ref_start))
    });

    let mut filtered = Vec::new();
    let mut used = Vec::new();

    for mem in mems.drain(..) {
        let overlaps = used.iter().any(|u: &MEM| {
            mem.ref_start < u.ref_end() && mem.ref_end() > u.ref_start
        });

        if !overlaps {
            filtered.push(mem.clone());
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
        let mut mems = vec![
            MEM::new(0, 0, 10),
            MEM::new(5, 5, 10),
            MEM::new(20, 20, 8),
        ];

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