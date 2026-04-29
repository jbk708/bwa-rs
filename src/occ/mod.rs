//! Occurrence table implementations using succinct data structures.

mod wavelet_tree;

pub use wavelet_tree::WaveletTree;

pub trait OccurrenceTable {
    fn occ(&self, c: u8, idx: usize) -> usize;
}

impl OccurrenceTable for WaveletTree {
    fn occ(&self, c: u8, idx: usize) -> usize {
        WaveletTree::rank(self, c, idx)
    }
}