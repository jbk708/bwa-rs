//! Reference genome handling with 2-bit encoding.

use crate::error::BwaError;
use crate::types::Sequence;

#[derive(Clone, Debug)]
pub struct Reference {
    pub name: String,
    pub contigs: Vec<Sequence>,
}

impl Reference {
    pub fn new(name: impl Into<String>, contigs: Vec<Sequence>) -> Self {
        Self {
            name: name.into(),
            contigs,
        }
    }

    pub fn total_len(&self) -> usize {
        self.contigs.iter().map(|c| c.len()).sum()
    }

    pub fn as_slice(&self) -> Vec<u8> {
        self.contigs
            .iter()
            .flat_map(|c| c.bases.iter())
            .copied()
            .collect()
    }

    pub fn from_fasta(path: impl AsRef<std::path::Path>) -> Result<Self, BwaError> {
        let content = std::fs::read_to_string(path)?;
        Self::parse_fasta(&content)
    }

    pub fn parse_fasta(content: &str) -> Result<Self, BwaError> {
        let mut name = String::new();
        let mut contigs = Vec::new();
        let mut current_bases = Vec::new();
        let mut current_name = String::new();

        for line in content.lines() {
            if line.is_empty() {
                continue;
            }
            if let Some(stripped) = line.strip_prefix('>') {
                if !current_name.is_empty() {
                    contigs.push(Sequence::new(
                        &current_name,
                        Self::encode_bases(&current_bases)?,
                    ));
                }
                current_name = stripped.split_whitespace().next().unwrap_or("").to_string();
                if name.is_empty() {
                    name = current_name.clone();
                }
                current_bases.clear();
            } else {
                current_bases.extend(line.chars().filter(|c| !c.is_whitespace()));
            }
        }

        if !current_name.is_empty() {
            contigs.push(Sequence::new(
                &current_name,
                Self::encode_bases(&current_bases)?,
            ));
        }

        Ok(Self { name, contigs })
    }

    fn encode_bases(bases: &[char]) -> Result<Vec<u8>, BwaError> {
        bases.iter().map(|b| Self::encode_base(*b)).collect()
    }

    pub fn encode_base(base: char) -> Result<u8, BwaError> {
        match base.to_ascii_uppercase() {
            'A' => Ok(0b00),
            'C' => Ok(0b01),
            'G' => Ok(0b10),
            'T' => Ok(0b11),
            'N' => Ok(0b100),
            _ => Err(BwaError::Parse(format!("Invalid base: {}", base))),
        }
    }

    pub fn decode_base(base: u8) -> char {
        match base {
            0 => 'A',
            1 => 'C',
            2 => 'G',
            3 => 'T',
            _ => 'N',
        }
    }

    pub fn decode_sequence(bases: &[u8]) -> String {
        bases.iter().map(|&b| Self::decode_base(b)).collect()
    }

    pub fn subsequence(&self, contig_idx: usize, start: usize, len: usize) -> Option<Vec<u8>> {
        let contig = self.contigs.get(contig_idx)?;
        let end = (start + len).min(contig.len());
        Some(contig.bases[start..end].to_vec())
    }

    /// The forward concatenation followed by its reverse complement (the "2N"
    /// sequence), matching how bwa / bwa-mem2 build their FM-index. A single
    /// forward search over this sequence finds matches on both strands: a hit at
    /// SA position `p < N` is a forward-strand match, while `p >= N` is on the
    /// reverse strand. Layout: `B[N + i] = complement(fwd[N - 1 - i])`.
    pub fn as_slice_2n(&self) -> Vec<u8> {
        let fwd = self.as_slice();
        let n = fwd.len();
        let mut out = Vec::with_capacity(2 * n);
        out.extend_from_slice(&fwd);
        out.extend(reverse_complement(&fwd));
        out
    }
}

/// Reverse-complement a 2-bit-encoded sequence (A=0,C=1,G=2,T=3,N=4).
/// Complement swaps A<->T (0<->3) and C<->G (1<->2); N maps to N.
pub fn reverse_complement(seq: &[u8]) -> Vec<u8> {
    seq.iter()
        .rev()
        .map(|&b| match b {
            0 => 3,
            1 => 2,
            2 => 1,
            3 => 0,
            other => other, // N (4) and anything else stays unchanged
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_encode_decode() {
        assert_eq!(Reference::encode_base('A').unwrap(), 0);
        assert_eq!(Reference::encode_base('C').unwrap(), 1);
        assert_eq!(Reference::encode_base('G').unwrap(), 2);
        assert_eq!(Reference::encode_base('T').unwrap(), 3);
        assert_eq!(Reference::decode_base(0), 'A');
        assert_eq!(Reference::decode_base(3), 'T');
    }

    #[test]
    fn test_fasta_parsing() {
        let fasta = ">seq1\nACGT\n>seq2\nTGCA";
        let ref_seq = Reference::parse_fasta(fasta).unwrap();
        assert_eq!(ref_seq.contigs.len(), 2);
        assert_eq!(ref_seq.contigs[0].name, "seq1");
        assert_eq!(ref_seq.contigs[1].name, "seq2");
    }
}
