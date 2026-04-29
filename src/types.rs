//! Core types for sequence representation and alignment records.

use std::fmt::{self, Display};

fn decode_base(base: u8) -> char {
    match base {
        0 => 'A',
        1 => 'C',
        2 => 'G',
        3 => 'T',
        _ => 'N',
    }
}

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

    pub fn reference_length(&self) -> usize {
        self.ops
            .iter()
            .filter(|(op, _)| matches!(op, CigarOp::M | CigarOp::D | CigarOp::N | CigarOp::Eq | CigarOp::X))
            .map(|(_, len)| *len as usize)
            .sum()
    }

    pub fn extend(&mut self, other: Cigar) {
        if other.ops.is_empty() {
            return;
        }
        if self.ops.is_empty() {
            *self = other;
            return;
        }
        if let Some(last) = self.ops.last_mut() {
            if let Some(first) = other.ops.first() {
                if last.0 == first.0 {
                    last.1 += first.1;
                    self.ops.extend(other.ops.iter().skip(1));
                    return;
                }
            }
        }
        self.ops.extend(other.ops);
    }
}

impl Default for Cigar {
    fn default() -> Self {
        Self::new()
    }
}

impl Display for Cigar {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        for (op, len) in &self.ops {
            write!(f, "{}{}", len, char_from_op(*op))?;
        }
        Ok(())
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
    pub md_tag: Option<String>,
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
            md_tag: None,
        }
    }

    pub fn mdz_string(&self, read: &[u8], reference: &[u8]) -> String {
        let mut result = String::new();
        let mut q_pos = 0usize;
        let mut r_pos = self.position;
        let mut match_count = 0usize;

        for (op, len) in &self.cigar.ops {
            match op {
                CigarOp::Eq | CigarOp::X | CigarOp::M => {
                    for i in 0..*len as usize {
                        let q_idx = q_pos + i;
                        let r_idx = r_pos + i;

                        if q_idx >= read.len() || r_idx >= reference.len() {
                            continue;
                        }

                        if read[q_idx] == reference[r_idx] {
                            match_count += 1;
                        } else {
                            if match_count > 0 {
                                result.push_str(&match_count.to_string());
                                match_count = 0;
                            }
                            result.push(decode_base(reference[r_idx]));
                        }
                    }
                    q_pos += *len as usize;
                    r_pos += *len as usize;
                }
                CigarOp::D => {
                    if match_count > 0 {
                        result.push_str(&match_count.to_string());
                        match_count = 0;
                    }
                    result.push('^');
                    for i in 0..*len as usize {
                        let r_idx = r_pos + i;
                        if r_idx < reference.len() {
                            result.push(decode_base(reference[r_idx]));
                        }
                    }
                    r_pos += *len as usize;
                }
                CigarOp::I => {
                    q_pos += *len as usize;
                }
                _ => {}
            }
        }

        if match_count > 0 {
            result.push_str(&match_count.to_string());
        }

        if result.is_empty() {
            result.push('0');
        }

        result
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

    #[test]
    fn test_mdz_perfect_match() {
        let mut cigar = Cigar::new();
        cigar.push(CigarOp::Eq, 5);
        let result = AlignmentResult {
            position: 0,
            mapq: 60,
            cigar,
            flag: 0,
            reverse_strand: false,
            nm: 0,
            score: 5,
            md_tag: None,
        };

        let read = vec![0, 1, 2, 3, 0];
        let reference = vec![0, 1, 2, 3, 0];
        let mdz = result.mdz_string(&read, &reference);
        assert_eq!(mdz, "5");
    }

    #[test]
    fn test_mdz_single_mismatch() {
        let mut cigar = Cigar::new();
        cigar.push(CigarOp::Eq, 3);
        cigar.push(CigarOp::X, 2);
        let result = AlignmentResult {
            position: 0,
            mapq: 60,
            cigar,
            flag: 0,
            reverse_strand: false,
            nm: 2,
            score: 3,
            md_tag: None,
        };

        let read = vec![0, 1, 2, 1, 2]; // A, C, G, C, G
        let reference = vec![0, 1, 2, 3, 2]; // A, C, G, T, G
        let mdz = result.mdz_string(&read, &reference);
        assert_eq!(mdz, "3T1");
    }

    #[test]
    fn test_mdz_deletion() {
        let mut cigar = Cigar::new();
        cigar.push(CigarOp::Eq, 2);
        cigar.push(CigarOp::D, 2);
        cigar.push(CigarOp::Eq, 2);
        let result = AlignmentResult {
            position: 0,
            mapq: 60,
            cigar,
            flag: 0,
            reverse_strand: false,
            nm: 2,
            score: 4,
            md_tag: None,
        };

        // CIGAR: 2= 2D 2=
        // Read: A C _ _ T A (read[0,1] match, read[2,3] match after deletion)
        // Ref:  A C G T T A (ref[2,3]=GT deleted, ref[4,5]=TA match read[2,3])
        let read = vec![0, 1, 3, 0]; // ACTA
        let reference = vec![0, 1, 2, 3, 3, 0]; // ACGTTA encoding: A=0,C=1,G=2,T=3,T=3,A=0
        let mdz = result.mdz_string(&read, &reference);
        assert_eq!(mdz, "2^GT2");
    }

    #[test]
    fn test_mdz_mixed() {
        let mut cigar = Cigar::new();
        cigar.push(CigarOp::Eq, 2);
        cigar.push(CigarOp::X, 1);
        cigar.push(CigarOp::Eq, 2);
        cigar.push(CigarOp::D, 1);
        let result = AlignmentResult {
            position: 0,
            mapq: 60,
            cigar,
            flag: 0,
            reverse_strand: false,
            nm: 3,
            score: 6,
            md_tag: None,
        };

        // CIGAR: 2= 1X 2= 1D
        // Read: A C G T A (length 5)
        // Ref:  A C T T A C (ref[2]=T mismatch, ref[3,4]=TA match, ref[5]=C del)
        let read = vec![0, 1, 2, 3, 0]; // ACGTA
        let reference = vec![0, 1, 3, 3, 0, 1]; // ACTTAC encoding: A=0,C=1,T=3,T=3,A=0,C=1
        // After 2=: ref[0,1]=AC match, q_pos=2, r_pos=2
        // 1X: ref[2]=T vs G mismatch, q_pos=3, r_pos=3
        // 2=: ref[3,4]=TA match, q_pos=5, r_pos=5
        // 1D: ref[5]=C deleted, r_pos=6
        let mdz = result.mdz_string(&read, &reference);
        assert_eq!(mdz, "2T2^C");
    }

    #[test]
    fn test_mdz_insertion() {
        let mut cigar = Cigar::new();
        cigar.push(CigarOp::Eq, 5);
        let result = AlignmentResult {
            position: 0,
            mapq: 60,
            cigar,
            flag: 0,
            reverse_strand: false,
            nm: 0,
            score: 6,
            md_tag: None,
        };

        // 5= means 5 matches against reference
        let read = vec![0, 1, 2, 3, 0]; // ACGTA
        let reference = vec![0, 1, 2, 3, 0, 1]; // ACGTAC
        let mdz = result.mdz_string(&read, &reference);
        assert_eq!(mdz, "5");
    }

    #[test]
    fn test_mdz_no_match() {
        let mut cigar = Cigar::new();
        cigar.push(CigarOp::X, 3);
        let result = AlignmentResult {
            position: 0,
            mapq: 60,
            cigar,
            flag: 0,
            reverse_strand: false,
            nm: 3,
            score: -12,
            md_tag: None,
        };

        let read = vec![0, 1, 2]; // ACG
        let reference = vec![3, 2, 1]; // TGC
        let mdz = result.mdz_string(&read, &reference);
        assert_eq!(mdz, "TGC");
    }
}