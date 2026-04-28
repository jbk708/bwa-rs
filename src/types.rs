//! Core types for sequence representation and alignment records.

use std::fmt;

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Sequence {
    pub name: String,
    pub bases: Vec<u8>,
}

impl Sequence {
    pub fn new(name: impl Into<String>, bases: Vec<u8>) -> Self {
        Self {
            name: name.into(),
            bases,
        }
    }

    pub fn len(&self) -> usize {
        self.bases.len()
    }

    pub fn is_empty(&self) -> bool {
        self.bases.is_empty()
    }

    pub fn get(&self, pos: usize) -> Option<u8> {
        self.bases.get(pos).copied()
    }

    pub fn as_slice(&self) -> &[u8] {
        &self.bases
    }
}

#[derive(Clone, Debug)]
pub struct Quality {
    pub scores: Vec<u8>,
}

impl Quality {
    pub fn new(scores: Vec<u8>) -> Self {
        Self { scores }
    }

    pub fn len(&self) -> usize {
        self.scores.len()
    }

    pub fn is_empty(&self) -> bool {
        self.scores.is_empty()
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Orientation {
    FR,
    FF,
    RR,
    RF,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Strand {
    Forward,
    Reverse,
}

impl fmt::Display for Strand {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Strand::Forward => write!(f, "+"),
            Strand::Reverse => write!(f, "-"),
        }
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum CigarOp {
    M,
    I,
    D,
    N,
    S,
    H,
    P,
    Eq,
    X,
}

#[derive(Clone, Debug)]
pub struct Cigar {
    pub ops: Vec<(CigarOp, u32)>,
}

impl Cigar {
    pub fn new() -> Self {
        Self { ops: Vec::new() }
    }

    pub fn push(&mut self, op: CigarOp, len: u32) {
        self.ops.push((op, len));
    }

    pub fn to_string(&self) -> String {
        self.ops
            .iter()
            .map(|(op, len)| format!("{}{}", len, char_from_op(*op)))
            .collect()
    }

    pub fn reference_length(&self) -> usize {
        self.ops
            .iter()
            .filter(|(op, _)| matches!(op, CigarOp::M | CigarOp::D | CigarOp::N | CigarOp::Eq | CigarOp::X))
            .map(|(_, len)| *len as usize)
            .sum()
    }
}

impl Default for Cigar {
    fn default() -> Self {
        Self::new()
    }
}

fn char_from_op(op: CigarOp) -> char {
    match op {
        CigarOp::M => 'M',
        CigarOp::I => 'I',
        CigarOp::D => 'D',
        CigarOp::N => 'N',
        CigarOp::S => 'S',
        CigarOp::H => 'H',
        CigarOp::P => 'P',
        CigarOp::Eq => '=',
        CigarOp::X => 'X',
    }
}

#[derive(Clone, Debug)]
pub struct MEM {
    pub query_start: usize,
    pub ref_start: usize,
    pub length: usize,
    pub score: f32,
}

impl MEM {
    pub fn new(query_start: usize, ref_start: usize, length: usize) -> Self {
        Self {
            query_start,
            ref_start,
            length,
            score: length as f32,
        }
    }

    pub fn query_end(&self) -> usize {
        self.query_start + self.length
    }

    pub fn ref_end(&self) -> usize {
        self.ref_start + self.length
    }
}

#[derive(Clone, Debug)]
pub struct ChainedSeed {
    pub mem: MEM,
    pub score: f32,
    pub forward_score: f32,
    pub backward_score: f32,
}

#[derive(Clone, Debug)]
pub struct AlignmentResult {
    pub position: usize,
    pub mapq: u8,
    pub cigar: Cigar,
    pub flag: u16,
    pub reverse_strand: bool,
    pub nm: u32,
    pub score: i32,
}

impl AlignmentResult {
    pub fn new(position: usize, cigar: Cigar) -> Self {
        Self {
            position,
            mapq: 60,
            cigar,
            flag: 0,
            reverse_strand: false,
            nm: 0,
            score: 0,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_sequence_encoding() {
        let seq = Sequence::new("test", vec![0, 1, 2, 3]);
        assert_eq!(seq.len(), 4);
        assert_eq!(seq.get(0), Some(0));
    }

    #[test]
    fn test_cigar() {
        let mut cigar = Cigar::new();
        cigar.push(CigarOp::M, 10);
        cigar.push(CigarOp::I, 2);
        cigar.push(CigarOp::M, 5);
        assert_eq!(cigar.to_string(), "10M2I5M");
        assert_eq!(cigar.reference_length(), 15);
    }

    #[test]
    fn test_mem() {
        let mem = MEM::new(5, 100, 20);
        assert_eq!(mem.query_end(), 25);
        assert_eq!(mem.ref_end(), 120);
    }
}