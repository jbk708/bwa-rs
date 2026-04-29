//! Occurrence table implementations using succinct data structures.

mod rrr;
mod wavelet_tree;

pub use rrr::RrrBitvec;
pub use wavelet_tree::WaveletTree;

pub trait OccurrenceTable {
    fn occ(&self, c: u8, idx: usize) -> usize;
}

impl OccurrenceTable for WaveletTree {
    fn occ(&self, c: u8, idx: usize) -> usize {
        WaveletTree::rank(self, c, idx)
    }
}

impl OccurrenceTable for RrrBitvec {
    fn occ(&self, c: u8, idx: usize) -> usize {
        RrrBitvec::rank(self, c, idx)
    }
}
