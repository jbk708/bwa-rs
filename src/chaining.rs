//! Seed chaining for combining MEMs into full alignments.

use crate::types::{ChainedSeed, MEM};

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
}
