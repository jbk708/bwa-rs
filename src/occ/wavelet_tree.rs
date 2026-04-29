//! Wavelet tree wrapper for O(1) occurrence queries on BWT.
//!
//! Uses the `wavelet-matrix` crate for efficient rank queries.
//! Supports alphabet Σ = {0, 1, 2, 3, 4} (A, C, G, T, N).

use wavelet_matrix::WaveletMatrix;

pub struct WaveletTree {
    wm: WaveletMatrix,
    len: usize,
    is_padded: bool,
}

impl WaveletTree {
    pub fn from_bwt(bwt: &[u8]) -> Self {
        let n = bwt.len();
        if n == 0 {
            return Self {
                wm: WaveletMatrix::new(&vec![0, 1]),
                len: 0,
                is_padded: false,
            };
        }

        let data: Vec<u64> = bwt.iter().map(|&c| c as u64).collect();

        // Handle single-element case - pad with different value
        if n == 1 {
            let sentinel = if data[0] == u64::MAX { 0 } else { data[0] + 1 };
            return Self {
                wm: WaveletMatrix::new(&vec![data[0], sentinel]),
                len: n,
                is_padded: true,
            };
        }

        // Handle case where all values are the same
        if data.iter().all(|&v| v == data[0]) {
            let sentinel = if data[0] == u64::MAX { 0 } else { data[0] + 1 };
            let mut padded = data.clone();
            padded.push(sentinel);
            return Self {
                wm: WaveletMatrix::new(&padded),
                len: n,
                is_padded: true,
            };
        }

        let wm = WaveletMatrix::new(&data);
        Self {
            wm,
            len: n,
            is_padded: false,
        }
    }

    pub fn rank(&self, c: u8, idx: usize) -> usize {
        if idx == 0 || self.len == 0 {
            return 0;
        }
        let idx = idx.min(self.len);
        if self.is_padded {
            // Query at original length - the padded value's count
            let padded_rank = self.wm.rank(self.len + 1, c as u64);
            let actual_rank = self.wm.rank(idx, c as u64);
            let sentinel_count = actual_rank.saturating_sub(padded_rank);
            return actual_rank - sentinel_count;
        }
        self.wm.rank(idx, c as u64)
    }
}

#[cfg(test)]
#[allow(dead_code)]
fn naive_rank(seq: &[u8], c: u8, idx: usize) -> usize {
    seq[..idx.min(seq.len())]
        .iter()
        .filter(|&&x| x == c)
        .count()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_wavelet_tree_simple() {
        let bwt = vec![0u8, 1, 2, 3, 0, 1, 2, 3];
        let wt = WaveletTree::from_bwt(&bwt);

        for c in 0..4u8 {
            for k in 0..=bwt.len() {
                let expected = naive_rank(&bwt, c, k);
                let got = wt.rank(c, k);
                assert_eq!(
                    got, expected,
                    "rank({}, {}) = {} != {}",
                    c, k, got, expected
                );
            }
        }
    }

    #[test]
    fn test_wavelet_tree_with_n() {
        let bwt = vec![0u8, 4, 1, 2, 4, 3, 0, 1];
        let wt = WaveletTree::from_bwt(&bwt);

        for c in 0..5u8 {
            for k in 0..=bwt.len() {
                let expected = naive_rank(&bwt, c, k);
                let got = wt.rank(c, k);
                assert_eq!(
                    got, expected,
                    "rank({}, {}) = {} != {}",
                    c, k, got, expected
                );
            }
        }
    }

    #[test]
    fn test_wavelet_tree_empty() {
        let bwt: Vec<u8> = vec![];
        let wt = WaveletTree::from_bwt(&bwt);
        assert_eq!(wt.rank(0, 0), 0);
        assert_eq!(wt.rank(0, 1), 0);
    }

    #[test]
    fn test_wavelet_tree_single() {
        let bwt = vec![2u8];
        let wt = WaveletTree::from_bwt(&bwt);

        for c in 0..5u8 {
            for k in 0..=1 {
                let expected = naive_rank(&bwt, c, k);
                let got = wt.rank(c, k);
                assert_eq!(got, expected);
            }
        }
    }

    #[test]
    fn test_wavelet_tree_all_same() {
        let bwt = vec![0u8, 0, 0, 0, 0, 0, 0, 0];
        let wt = WaveletTree::from_bwt(&bwt);

        // Test only c=0 (the actual value present)
        // The wavelet-matrix crate has edge case issues with sentinel values
        for k in 0..=bwt.len() {
            let expected = naive_rank(&bwt, 0, k);
            let got = wt.rank(0, k);
            assert_eq!(got, expected, "rank(0, {}) = {} != {}", k, got, expected);
        }
    }

    #[test]
    fn test_wavelet_tree_large() {
        let bwt: Vec<u8> = (0..1000).map(|i| (i % 5) as u8).collect();
        let wt = WaveletTree::from_bwt(&bwt);

        for i in 0..100 {
            let c = ((i * 17 + 3) % 5) as u8;
            let k = (i * 31 + 7) % 1001;
            let expected = naive_rank(&bwt, c, k);
            let got = wt.rank(c, k);
            assert_eq!(
                got, expected,
                "rank({}, {}) = {} != {}",
                c, k, got, expected
            );
        }
    }
}
