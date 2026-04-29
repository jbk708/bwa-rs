//! RRR bitvector using succinct crate's Rank9 for O(1) rank queries.

use succinct::rank::Rank9;
use succinct::BitRankSupport;

#[derive(Default)]
pub struct RrrBitvec {
    bitvecs: [Option<Rank9<Vec<u64>>>; 5],
    len: usize,
}

impl RrrBitvec {
    pub fn from_bwt(bwt: &[u8]) -> Self {
        let n = bwt.len();
        if n == 0 {
            return Self::default();
        }

        let num_blocks = n.div_ceil(64);
        let mut bitvec_data = [
            vec![0u64; num_blocks],
            vec![0u64; num_blocks],
            vec![0u64; num_blocks],
            vec![0u64; num_blocks],
            vec![0u64; num_blocks],
        ];

        for (i, &c) in bwt.iter().enumerate() {
            bitvec_data[c as usize][i / 64] |= 1u64 << (i % 64);
        }

        Self {
            bitvecs: bitvec_data.map(|data| Some(Rank9::new(data))),
            len: n,
        }
    }

    pub fn rank(&self, c: u8, idx: usize) -> usize {
        if idx == 0 || self.len == 0 || c >= 5 {
            return 0;
        }
        self.bitvecs[c as usize]
            .as_ref()
            .map(|rank| rank.rank1(idx.min(self.len) as u64 - 1) as usize)
            .unwrap_or(0)
    }
}

#[cfg(test)]
#[allow(dead_code)]
fn naive_rank(seq: &[u8], c: u8, idx: usize) -> usize {
    seq[..idx.min(seq.len())].iter().filter(|&&x| x == c).count()
}

#[cfg(test)]
mod tests {
    use super::*;

    fn check_rank(rrr: &RrrBitvec, bwt: &[u8]) {
        for c in 0..5u8 {
            for k in 0..=bwt.len() {
                assert_eq!(rrr.rank(c, k), naive_rank(bwt, c, k));
            }
        }
    }

    #[test]
    fn test_rrr_simple() {
        let bwt = vec![0u8, 1, 2, 3, 0, 1, 2, 3];
        check_rank(&RrrBitvec::from_bwt(&bwt), &bwt);
    }

    #[test]
    fn test_rrr_with_n() {
        let bwt = vec![0u8, 4, 1, 2, 4, 3, 0, 1];
        check_rank(&RrrBitvec::from_bwt(&bwt), &bwt);
    }

    #[test]
    fn test_rrr_empty() {
        let rrr = RrrBitvec::from_bwt(&[]);
        assert_eq!(rrr.rank(0, 0), 0);
        assert_eq!(rrr.rank(0, 1), 0);
    }

    #[test]
    fn test_rrr_single() {
        let bwt = vec![2u8];
        check_rank(&RrrBitvec::from_bwt(&bwt), &bwt);
    }

    #[test]
    fn test_rrr_all_same() {
        let bwt = vec![0u8; 8];
        check_rank(&RrrBitvec::from_bwt(&bwt), &bwt);
    }

    #[test]
    fn test_rrr_large() {
        let bwt: Vec<u8> = (0..1000).map(|i| (i % 5) as u8).collect();
        let rrr = RrrBitvec::from_bwt(&bwt);

        for i in 0..100 {
            let c = ((i * 17 + 3) % 5) as u8;
            let k = (i * 31 + 7) % 1001;
            assert_eq!(rrr.rank(c, k), naive_rank(&bwt, c, k));
        }
    }

    #[test]
    fn test_rrr_compression() {
        use succinct::SpaceUsage;

        let n = 100_000;
        let bwt: Vec<u8> = (0..n).map(|i| (i % 5) as u8).collect();
        let rrr = RrrBitvec::from_bwt(&bwt);

        let rrr_bytes: usize = rrr.bitvecs.iter()
            .filter_map(|bv| bv.as_ref())
            .map(|bv| bv.heap_bytes())
            .sum();

        let ratio = rrr_bytes as f64 / n as f64;
        assert!(ratio < 10.0, "Compression ratio too high: {ratio:.2}x");
    }
}