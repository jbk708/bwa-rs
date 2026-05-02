//! Utility functions for sequence encoding/decoding.
//!
//! 2-bit encoding: A=0, C=1, G=2, T=3, N=4

use crate::error::BwaError;

pub fn encode_base(base: char) -> Result<u8, BwaError> {
    match base.to_ascii_uppercase() {
        'A' => Ok(0),
        'C' => Ok(1),
        'G' => Ok(2),
        'T' => Ok(3),
        'N' => Ok(4),
        c => Err(BwaError::Parse(format!("Invalid base: {}", c))),
    }
}

pub fn encode_sequence(seq: &[u8]) -> Vec<u8> {
    seq.iter()
        .map(|&b| match b.to_ascii_uppercase() {
            b'A' => 0,
            b'C' => 1,
            b'G' => 2,
            b'T' => 3,
            b'N' => 4,
            _ => 4,
        })
        .collect()
}

pub fn encode_sequence_u16(seq: &[u8]) -> Vec<u16> {
    seq.iter()
        .map(|&b| match b.to_ascii_uppercase() {
            b'A' => 0,
            b'C' => 1,
            b'G' => 2,
            b'T' => 3,
            b'N' => 4,
            _ => 4,
        })
        .collect()
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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_encode_base() {
        assert_eq!(encode_base('A').unwrap(), 0);
        assert_eq!(encode_base('C').unwrap(), 1);
        assert_eq!(encode_base('G').unwrap(), 2);
        assert_eq!(encode_base('T').unwrap(), 3);
        assert_eq!(encode_base('N').unwrap(), 4);
        assert_eq!(encode_base('a').unwrap(), 0);
        assert_eq!(encode_base('c').unwrap(), 1);
        assert!(encode_base('X').is_err());
    }

    #[test]
    fn test_encode_sequence_basic() {
        let seq = b"ACGT";
        let encoded = encode_sequence(seq);
        assert_eq!(encoded, &[0, 1, 2, 3]);
    }

    #[test]
    fn test_encode_sequence_with_n() {
        let seq = b"ACGNACGT";
        let encoded = encode_sequence(seq);
        assert_eq!(encoded, &[0, 1, 2, 4, 0, 1, 2, 3]);
    }

    #[test]
    fn test_encode_sequence_lowercase() {
        let seq = b"acgt";
        let encoded = encode_sequence(seq);
        assert_eq!(encoded, &[0, 1, 2, 3]);
    }

    #[test]
    fn test_encode_sequence_empty() {
        let seq: &[u8] = b"";
        let encoded = encode_sequence(seq);
        assert!(encoded.is_empty());
    }

    #[test]
    fn test_encode_sequence_unknown_char() {
        let seq = b"ACXGT";
        let encoded = encode_sequence(seq);
        assert_eq!(encoded, &[0, 1, 4, 2, 3]);
    }

    #[test]
    fn test_encode_sequence_u16() {
        let seq = b"ACGT";
        let encoded = encode_sequence_u16(seq);
        assert_eq!(encoded, &[0u16, 1, 2, 3]);
    }

    #[test]
    fn test_decode_base() {
        assert_eq!(decode_base(0), 'A');
        assert_eq!(decode_base(1), 'C');
        assert_eq!(decode_base(2), 'G');
        assert_eq!(decode_base(3), 'T');
        assert_eq!(decode_base(4), 'N');
        assert_eq!(decode_base(255), 'N');
    }

    #[test]
    fn test_encode_decode_roundtrip() {
        for base in ['A', 'C', 'G', 'T', 'N', 'a', 'c'] {
            let encoded = encode_base(base).unwrap();
            let decoded = decode_base(encoded);
            assert_eq!(
                decoded,
                base.to_ascii_uppercase(),
                "Roundtrip failed for {}",
                base
            );
        }
    }
}