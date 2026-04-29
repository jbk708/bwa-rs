use std::io::{BufRead, BufReader, Read};
use std::path::Path;

use crate::types::Sequence;
use crate::BwaError;

pub struct FASTQReader<R: Read> {
    reader: BufReader<R>,
    buffer: String,
}

impl<R: Read> FASTQReader<R> {
    pub fn new(input: R) -> Self {
        Self {
            reader: BufReader::new(input),
            buffer: String::new(),
        }
    }
}

impl FASTQReader<BufReader<std::fs::File>> {
    pub fn from_path(path: &Path) -> Result<Self, BwaError> {
        let file = std::fs::File::open(path)?;
        Ok(Self::new(BufReader::new(file)))
    }
}

impl<R: Read> Iterator for FASTQReader<R> {
    type Item = Result<FastqRecord, BwaError>;

    fn next(&mut self) -> Option<Self::Item> {
        self.buffer.clear();
        
        let name = match self.reader.read_line(&mut self.buffer) {
            Ok(0) => return None,
            Ok(_) => self.buffer.trim().to_string(),
            Err(e) => return Some(Err(BwaError::Io(e))),
        };

        if !name.starts_with('@') {
            return Some(Err(BwaError::Parse(format!(
                "Expected '@' at record start, got: {}",
                name
            ))));
        }
        let qname = name.trim_start_matches('@').to_string();

        self.buffer.clear();
        if let Err(e) = self.reader.read_line(&mut self.buffer) {
            return Some(Err(BwaError::Io(e)));
        }
        let seq = self.buffer.trim().to_string();

        self.buffer.clear();
        if let Err(e) = self.reader.read_line(&mut self.buffer) {
            return Some(Err(BwaError::Io(e)));
        }
        if !self.buffer.trim().starts_with('+') {
            return Some(Err(BwaError::Parse(format!(
                "Expected '+' separator, got: {}",
                self.buffer.trim()
            ))));
        }

        self.buffer.clear();
        if let Err(e) = self.reader.read_line(&mut self.buffer) {
            return Some(Err(BwaError::Io(e)));
        }
        let qual = self.buffer.trim().to_string();

        if seq.len() != qual.len() {
            return Some(Err(BwaError::Parse(format!(
                "Sequence and quality length mismatch: {} vs {}",
                seq.len(),
                qual.len()
            ))));
        }

        Some(Ok(FastqRecord {
            qname,
            seq,
            qual,
        }))
    }
}

pub struct FastqRecord {
    pub qname: String,
    pub seq: String,
    pub qual: String,
}

impl FastqRecord {
    pub fn to_sequence(&self) -> Sequence {
        let bases: Vec<u8> = self.seq
            .bytes()
            .map(|b| match b {
                b'A' | b'a' => 0,
                b'C' | b'c' => 1,
                b'G' | b'g' => 2,
                b'T' | b't' => 3,
                _ => 4,
            })
            .collect();

        Sequence::new(&self.qname, bases)
    }

    pub fn to_quality(&self) -> Vec<u8> {
        self.qual.bytes().map(|b| b.wrapping_sub(33)).collect()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Cursor;

    #[test]
    fn test_parse_simple() {
        let data = b"@read1\nACGT\n+\nIIII\n";
        let reader = FASTQReader::new(Cursor::new(&data[..]));
        
        let records: Vec<_> = reader.collect();
        assert_eq!(records.len(), 1);
        let record = records[0].as_ref().unwrap();
        assert_eq!(record.qname, "read1");
        assert_eq!(record.seq, "ACGT");
    }

    #[test]
    fn test_to_sequence() {
        let record = FastqRecord {
            qname: "test".to_string(),
            seq: "ACGT".to_string(),
            qual: "IIII".to_string(),
        };
        let seq = record.to_sequence();
        assert_eq!(&seq.bases, &[0, 1, 2, 3]);
    }
}