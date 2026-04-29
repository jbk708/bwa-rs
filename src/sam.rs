//! SAM format output generation.

use crate::reference::Reference;
use crate::types::{AlignmentResult, Sequence};
use std::fmt;
use std::io::{self, Write};

#[derive(Clone)]
pub struct SAMHeader {
    pub version: String,
    pub sort_order: String,
    pub reference: Reference,
}

impl fmt::Display for SAMHeader {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        writeln!(f, "@HD\tVN:{}\tSO:{}", self.version, self.sort_order)?;
        for contig in &self.reference.contigs {
            writeln!(f, "@SQ\tSN:{}\tLN:{}", contig.name, contig.len())?;
        }
        writeln!(f, "@PG\tID:bwa-rs\tPN:bwa-rs\tVN:0.1.0\tCL:bwa-rs")
    }
}

impl SAMHeader {
    pub fn new(reference: Reference) -> Self {
        Self {
            version: "1.6".to_string(),
            sort_order: "coordinate".to_string(),
            reference,
        }
    }
}

#[derive(Clone, Debug)]
pub struct SAMRecord {
    pub qname: String,
    pub flag: u16,
    pub rname: String,
    pub pos: u32,
    pub mapq: u8,
    pub cigar: String,
    pub rnext: String,
    pub pnext: u32,
    pub tlen: i32,
    pub seq: String,
    pub qual: String,
    pub md_tag: Option<String>,
}

impl fmt::Display for SAMRecord {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            self.qname, self.flag, self.rname, self.pos, self.mapq, self.cigar,
            self.rnext, self.pnext, self.tlen, self.seq, self.qual,
        )?;
        if let Some(ref md) = self.md_tag {
            write!(f, "\t{}", md)?;
        }
        Ok(())
    }
}

impl SAMRecord {
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        qname: String,
        flag: u16,
        rname: String,
        pos: u32,
        mapq: u8,
        cigar: String,
        rnext: String,
        pnext: u32,
        tlen: i32,
        seq: String,
        qual: String,
    ) -> Self {
        Self {
            qname,
            flag,
            rname,
            pos,
            mapq,
            cigar,
            rnext,
            pnext,
            tlen,
            seq,
            qual,
            md_tag: None,
        }
    }

    pub fn from_alignment(
        qname: &str,
        alignment: &AlignmentResult,
        seq: &Sequence,
        qual: &[u8],
        contig_name: &str,
    ) -> Self {
        let flag = alignment.flag
            | if alignment.reverse_strand { 0x10 } else { 0 };

        let (rname, pos, cigar) = if flag & 0x4 != 0 {
            ("*".to_string(), 0, "*".to_string())
        } else {
            (
                contig_name.to_string(),
                (alignment.position + 1) as u32,
                alignment.cigar.to_string(),
            )
        };

        let mut record = Self::new(
            qname.to_string(), flag, rname, pos, alignment.mapq, cigar,
            "*".to_string(), 0, 0,
            Reference::decode_sequence(&seq.bases),
            if qual.is_empty() { "*".to_string() } else { String::from_utf8_lossy(qual).to_string() },
        );
        record.md_tag = if flag & 0x4 != 0 { None } else { alignment.md_tag.clone() };
        record
    }

    pub fn unmapped(qname: &str, seq: &[u8], qual: &[u8]) -> Self {
        Self::new(
            qname.to_string(), 0x4, "*".to_string(), 0, 0, "*".to_string(),
            "*".to_string(), 0, 0,
            Reference::decode_sequence(seq),
            if qual.is_empty() { "*".to_string() } else { String::from_utf8_lossy(qual).to_string() },
        )
    }
}

pub struct SAMWriter {
    writer: Box<dyn Write>,
    header: SAMHeader,
}

impl SAMWriter {
    pub fn new(write: Box<dyn Write>, reference: Reference) -> io::Result<Self> {
        let mut sam_writer = Self { writer: write, header: SAMHeader::new(reference) };
        sam_writer.write_header()?;
        Ok(sam_writer)
    }

    pub fn from_path(path: &std::path::Path, reference: Reference) -> io::Result<Self> {
        let file = std::fs::File::create(path)?;
        Self::new(Box::new(file), reference)
    }

    pub fn to_stdout(reference: Reference) -> io::Result<Self> {
        Self::new(Box::new(io::stdout()), reference)
    }

    fn write_header(&mut self) -> io::Result<()> {
        write!(self.writer, "{}", self.header)
    }

    pub fn write_record(&mut self, record: &SAMRecord) -> io::Result<()> {
        writeln!(self.writer, "{}", record)
    }

    pub fn flush(&mut self) -> io::Result<()> {
        self.writer.flush()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_header_generation() {
        let ref_seq = Reference::parse_fasta(">chr1\nACGT").unwrap();
        let header = SAMHeader::new(ref_seq);
        let output = header.to_string();
        assert!(output.contains("@HD"));
        assert!(output.contains("@SQ"));
        assert!(output.contains("chr1"));
    }

    #[test]
    fn test_record_to_string() {
        let record = SAMRecord::new(
            "read1".to_string(), 0, "chr1".to_string(), 100, 60,
            "10M".to_string(), "*".to_string(), 0, 0,
            "ACGTACGT".to_string(), "IIIIIIII".to_string(),
        );
        let line = record.to_string();
        assert!(line.starts_with("read1"));
        assert!(line.contains("chr1"));
    }

    #[test]
    fn test_sam_record_fields() {
        let record = SAMRecord::new(
            "read1".to_string(), 99, "chr1".to_string(), 1000, 30,
            "50M".to_string(), "=".to_string(), 1050, 100,
            "ACGT".to_string(), "!@#$".to_string(),
        );
        assert_eq!(record.flag, 99);
        assert_eq!(record.pos, 1000);
        assert_eq!(record.mapq, 30);
        assert_eq!(record.cigar, "50M");
    }
}