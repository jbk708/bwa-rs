//! Seed chaining for combining MEMs into full alignments.

use crate::types::{ChainedSeed, MEM};

/// bwa `mem_chain_flt` `drop_ratio` — a chain is dropped when its weight is
/// below this fraction of an overlapping, higher-weight chain.
const CHAIN_DROP_RATIO: f32 = 0.50;

/// bwa `mask_level` — the query-span overlap fraction above which two chains are
/// considered to significantly overlap.
const CHAIN_MASK_LEVEL: f32 = 0.50;

/// bwa `max_chain_gap` — the largest query/reference gap that still chains.
const MAX_CHAIN_GAP: i64 = 10_000;

/// bwa chaining band width (`opt->w`) — the diagonal tolerance when merging a
/// seed into a chain.
const CHAIN_BAND_WIDTH: i64 = 100;

/// A group of collinear seeds, with bwa's `mem_chain_weight` and query span
/// (`chn_beg`/`chn_end`) precomputed for `mem_chain_flt`.
#[derive(Clone, Debug)]
pub struct Chain {
    pub seeds: Vec<MEM>,
    pub weight: i32,
    pub qbeg: usize,
    pub qend: usize,
}

/// bwa `mem_chain_weight`: the smaller of the query-span and reference-span
/// coverage of the chain's seeds, where coverage is the length of the union of
/// the seed intervals on that axis.
pub fn mem_chain_weight(seeds: &[MEM]) -> i32 {
    // Length of the union of half-open intervals `[b, e)` over the seeds, in the
    // order given. `w += e - max(b, prev_end)` whenever the interval extends past
    // the running end covers all three bwa cases (disjoint, overlapping, nested).
    fn axis_cover(intervals: impl Iterator<Item = (i64, i64)>) -> i64 {
        let mut w = 0i64;
        let mut prev_end = 0i64;
        for (b, e) in intervals {
            if e > prev_end {
                w += e - b.max(prev_end);
            }
            prev_end = prev_end.max(e);
        }
        w
    }

    let mut order: Vec<&MEM> = seeds.iter().collect();

    order.sort_by_key(|s| s.query_start);
    let w_q = axis_cover(
        order
            .iter()
            .map(|s| (s.query_start as i64, s.query_end() as i64)),
    );

    order.sort_by_key(|s| s.ref_start);
    let w_r = axis_cover(
        order
            .iter()
            .map(|s| (s.ref_start as i64, s.ref_end() as i64)),
    );

    w_q.min(w_r).min((1 << 30) - 1) as i32
}

fn make_chain(seeds: Vec<MEM>) -> Chain {
    let qbeg = seeds.iter().map(|s| s.query_start).min().unwrap_or(0);
    let qend = seeds.iter().map(|s| s.query_end()).max().unwrap_or(0);
    let weight = mem_chain_weight(&seeds);
    Chain {
        seeds,
        weight,
        qbeg,
        qend,
    }
}

/// Group collinear seeds into chains, faithful to bwa's `mem_chain` /
/// `test_and_merge`: a seed joins a chain when it lies within the diagonal band
/// (`|x - y| <= w`) and the query/reference gaps stay under `max_chain_gap`.
/// Seeds are processed by ascending reference position; cross-strand and
/// cross-locus seeds fall out automatically because their reference offset
/// breaks the band test. Containment-only seeds are absorbed without growing the
/// chain span.
pub fn build_candidate_chains(mut sorted: Vec<MEM>) -> Vec<Chain> {
    if sorted.is_empty() {
        return Vec::new();
    }

    sorted.sort_by_key(|s| (s.ref_start, s.query_start));

    let mut open: Vec<Vec<MEM>> = Vec::new();
    for seed in sorted {
        let mut merged = false;
        for chain in open.iter_mut() {
            let first = chain[0];
            let last = chain[chain.len() - 1];
            let qend = last.query_end() as i64;
            let rend = last.ref_end() as i64;

            // Contained in the chain's span: absorb without extending.
            if seed.query_start >= first.query_start
                && seed.query_end() <= qend as usize
                && seed.ref_start >= first.ref_start
                && seed.ref_end() <= rend as usize
            {
                merged = true;
                break;
            }

            let x = seed.query_start as i64 - last.query_start as i64;
            let y = seed.ref_start as i64 - last.ref_start as i64;
            if y >= 0
                && (x - y).abs() <= CHAIN_BAND_WIDTH
                && x - last.length as i64 <= MAX_CHAIN_GAP
                && y - last.length as i64 <= MAX_CHAIN_GAP
                && (seed.query_start as i64 >= qend || seed.ref_start as i64 >= rend)
            {
                chain.push(seed);
                merged = true;
                break;
            }
        }
        if !merged {
            open.push(vec![seed]);
        }
    }

    open.into_iter().map(make_chain).collect()
}

/// bwa `mem_chain_flt`: drop chains that significantly overlap a higher-weight
/// chain in query space and are much weaker (`w < best_w * drop_ratio` and
/// `best_w - w >= 2 * min_seed_len`). Survivors are returned in descending
/// weight order; the highest-weight chain is always kept.
pub fn mem_chain_flt(mut chains: Vec<Chain>, min_seed_len: usize) -> Vec<Chain> {
    if chains.len() <= 1 {
        return chains;
    }

    chains.sort_by(|a, b| b.weight.cmp(&a.weight));

    let mut kept: Vec<Chain> = Vec::with_capacity(chains.len());
    for cand in chains.into_iter() {
        let mut contained = false;
        for k in &kept {
            let b_max = k.qbeg.max(cand.qbeg);
            let e_min = k.qend.min(cand.qend);
            if e_min <= b_max {
                continue;
            }
            let overlap = (e_min - b_max) as f32;
            let min_l = (k.qend - k.qbeg).min(cand.qend - cand.qbeg);
            if overlap >= min_l as f32 * CHAIN_MASK_LEVEL
                && (min_l as i64) < MAX_CHAIN_GAP
                && (cand.weight as f32) < k.weight as f32 * CHAIN_DROP_RATIO
                && k.weight - cand.weight >= (min_seed_len as i32) << 1
            {
                contained = true;
                break;
            }
        }
        if !contained {
            kept.push(cand);
        }
    }

    kept
}

pub fn chain_seeds(seeds: &[MEM], gap_penalty: f32) -> Vec<ChainedSeed> {
    if seeds.is_empty() {
        return Vec::new();
    }

    let mut sorted = seeds.to_vec();
    sorted.sort_by_key(|s| s.ref_start);

    let mut chains = Vec::new();
    let mut current_chain = Vec::new();

    for seed in sorted {
        let chained = ChainedSeed {
            mem: seed,
            score: seed.score,
            forward_score: seed.score,
            backward_score: 0.0,
        };

        if current_chain.is_empty() {
            current_chain.push(chained);
        } else {
            let prev = &current_chain[current_chain.len() - 1];
            let gap = seed.ref_start as f32 - prev.mem.ref_end() as f32;

            if gap <= gap_penalty && seed.query_start >= prev.mem.query_end() {
                let extended = ChainedSeed {
                    mem: seed,
                    score: prev.score + seed.score - gap,
                    forward_score: prev.forward_score + seed.score,
                    backward_score: prev.backward_score,
                };
                current_chain.push(extended);
            } else {
                if !current_chain.is_empty() {
                    chains.push(finalize_chain(current_chain));
                }
                current_chain = vec![chained];
            }
        }
    }

    if !current_chain.is_empty() {
        chains.push(finalize_chain(current_chain));
    }

    chains
}

fn finalize_chain(chain: Vec<ChainedSeed>) -> ChainedSeed {
    chain
        .into_iter()
        .max_by(|a, b| a.score.partial_cmp(&b.score).unwrap())
        .unwrap()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::types::MEM;

    #[test]
    fn test_chain_seeds() {
        let seeds = vec![
            MEM::new(0, 100, 20),
            MEM::new(25, 130, 15),
            MEM::new(50, 200, 18),
        ];

        let chains = chain_seeds(&seeds, 50.0);
        assert_eq!(chains.len(), 2);
    }

    #[test]
    fn test_chain_overlapping_seeds() {
        let seeds = vec![MEM::new(0, 100, 30), MEM::new(20, 120, 25)];

        let chains = chain_seeds(&seeds, 100.0);
        assert!(!chains.is_empty());
    }

    #[test]
    fn mem_chain_weight_single_seed_is_length() {
        assert_eq!(mem_chain_weight(&[MEM::new(0, 1000, 40)]), 40);
    }

    #[test]
    fn mem_chain_weight_is_min_of_axes() {
        // Two collinear seeds with a query gap but a ref overlap: query coverage
        // 30+30=60, reference coverage union 30 + (1090..1100) overlap-trimmed.
        let seeds = vec![MEM::new(0, 1000, 30), MEM::new(60, 1090, 30)];
        // query union: [0,30) + [60,90) = 60; ref union: [1000,1030)+[1090,1120) = 60.
        assert_eq!(mem_chain_weight(&seeds), 60);
        // Overlapping query intervals trim the union to the min axis.
        let overlap = vec![MEM::new(0, 1000, 40), MEM::new(20, 1020, 40)];
        // query union [0,60) = 60; ref union [1000,1060) = 60.
        assert_eq!(mem_chain_weight(&overlap), 60);
    }

    #[test]
    fn build_candidate_chains_groups_collinear() {
        let seeds = vec![MEM::new(0, 1000, 30), MEM::new(40, 1040, 30)];
        let chains = build_candidate_chains(seeds);
        assert_eq!(chains.len(), 1);
        assert_eq!(chains[0].seeds.len(), 2);
    }

    #[test]
    fn build_candidate_chains_separates_distinct_loci() {
        // Same query span, two far-apart reference loci → two separate chains.
        let seeds = vec![MEM::new(0, 1000, 30), MEM::new(0, 9_000_000, 30)];
        let chains = build_candidate_chains(seeds);
        assert_eq!(chains.len(), 2);
    }

    #[test]
    fn mem_chain_flt_drops_weak_contained_chain() {
        // A long primary chain (weight ~120) and a short secondary (weight ~20)
        // overlapping it in query space at a distinct locus: bwa drops the secondary.
        let primary = make_chain(vec![MEM::new(0, 1000, 120)]);
        let secondary = make_chain(vec![MEM::new(10, 5_000_000, 20)]);
        let kept = mem_chain_flt(vec![primary, secondary], 19);
        assert_eq!(kept.len(), 1);
        assert_eq!(kept[0].weight, 120);
    }

    #[test]
    fn mem_chain_flt_keeps_comparable_chains() {
        // Two chains of comparable weight overlapping in query space: both survive
        // (the drop test requires the loser be < 0.5x and >= 2*min_seed_len weaker).
        let a = make_chain(vec![MEM::new(0, 1000, 90)]);
        let b = make_chain(vec![MEM::new(0, 5_000_000, 80)]);
        let kept = mem_chain_flt(vec![a, b], 19);
        assert_eq!(kept.len(), 2);
    }

    #[test]
    fn mem_chain_flt_keeps_non_overlapping_chains() {
        // Disjoint query spans: even a weak chain survives (no query overlap).
        let a = make_chain(vec![MEM::new(0, 1000, 100)]);
        let b = make_chain(vec![MEM::new(120, 5_000_000, 20)]);
        let kept = mem_chain_flt(vec![a, b], 19);
        assert_eq!(kept.len(), 2);
    }
}
