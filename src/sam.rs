//! SAM format output generation.

use crate::paired::MateFields;
use crate::reference::{reverse_complement, Reference};
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
            self.qname,
            self.flag,
            self.rname,
            self.pos,
            self.mapq,
            self.cigar,
            self.rnext,
            self.pnext,
            self.tlen,
            self.seq,
            self.qual,
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
        let flag = alignment.flag | if alignment.reverse_strand { 0x10 } else { 0 };

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
            qname.to_string(),
            flag,
            rname,
            pos,
            alignment.mapq,
            cigar,
            "*".to_string(),
            0,
            0,
            Reference::decode_sequence(&seq.bases),
            if qual.is_empty() {
                "*".to_string()
            } else {
                String::from_utf8_lossy(qual).to_string()
            },
        );
        record.md_tag = if flag & 0x4 != 0 {
            None
        } else {
            alignment.md_tag.clone()
        };
        record
    }

    pub fn unmapped(qname: &str, seq: &[u8], qual: &[u8]) -> Self {
        Self::new(
            qname.to_string(),
            0x4,
            "*".to_string(),
            0,
            0,
            "*".to_string(),
            "*".to_string(),
            0,
            0,
            Reference::decode_sequence(seq),
            if qual.is_empty() {
                "*".to_string()
            } else {
                String::from_utf8_lossy(qual).to_string()
            },
        )
    }
}

pub struct SAMWriter {
    writer: Box<dyn Write>,
    header: SAMHeader,
}

impl SAMWriter {
    pub fn new(write: Box<dyn Write>, reference: Reference) -> io::Result<Self> {
        let mut sam_writer = Self {
            writer: write,
            header: SAMHeader::new(reference),
        };
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

/// Returns the SAM SEQ and QUAL for a read oriented to the reference strand.
/// On the reverse strand SEQ is reverse-complemented and QUAL reversed. An empty
/// `qual` yields `*` (the SAM "no quality" sentinel).
pub fn oriented_seq_qual(bases: &[u8], qual: &str, reverse: bool) -> (String, String) {
    let seq = if bases.is_empty() {
        "*".to_string()
    } else if reverse {
        Reference::decode_sequence(&reverse_complement(bases))
    } else {
        Reference::decode_sequence(bases)
    };
    let qual_out = if qual.is_empty() {
        "*".to_string()
    } else if reverse {
        qual.chars().rev().collect()
    } else {
        qual.to_string()
    };
    (seq, qual_out)
}

/// Writes a single paired-end alignment record to any `Write` sink.
///
/// Uses `MateFields` (pre-assembled by [`crate::paired::mate_fields`]) for all flag/rnext/pnext/tlen
/// values so the caller does not need to recompute them. RNAME/POS follow the SAM
/// unmapped-mate convention: an unmapped read whose mate is mapped inherits the mate's
/// coordinate via `placed_pos` so coordinate-sorted tools can keep pairs adjacent.
pub fn write_paired_record<W: Write>(
    out: &mut W,
    reference: &Reference,
    qname: &str,
    result: &AlignmentResult,
    bases: &[u8],
    qual: &str,
    mf: &MateFields,
) -> io::Result<()> {
    let mapped = !result.is_unmapped();
    let (rname, pos, mapq): (&str, i64, u8) = if mapped {
        reference
            .locate(result.position)
            .map_or(("*", 0, 0), |(name, off)| {
                (name, off as i64 + 1, result.mapq)
            })
    } else if let Some(p) = mf.placed_pos {
        reference
            .locate(p)
            .map_or(("*", 0, 0), |(name, off)| (name, off as i64 + 1, 0))
    } else {
        ("*", 0, 0)
    };

    // mf.pnext is a 1-based global coordinate; recover the 0-based global
    // position so the mate's contig resolves the same way as the read's.
    let (rnext, pnext, tlen): (&str, i64, i64) = if mf.rnext == "*" || mf.pnext == 0 {
        ("*", 0, 0)
    } else {
        reference
            .locate((mf.pnext - 1) as usize)
            .map_or(("*", 0, 0), |(mate_name, mate_off)| {
                if mate_name == rname {
                    ("=", mate_off as i64 + 1, mf.tlen)
                } else {
                    (mate_name, mate_off as i64 + 1, 0)
                }
            })
    };

    let cigar = if mapped {
        result.cigar.to_sam_string()
    } else {
        "*".into()
    };

    let (seq, qual_out) = oriented_seq_qual(bases, qual, result.reverse_strand);

    write!(
        out,
        "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
        qname, mf.flag, rname, pos, mapq, cigar, rnext, pnext, tlen, seq, qual_out
    )?;
    if mapped {
        write!(out, "\tNM:i:{}", result.nm)?;
        if let Some(ref md) = result.md_tag {
            write!(out, "\t{}", md)?;
        }
        if let Some(ref mc) = mf.mc {
            write!(out, "\tMC:Z:{}", mc)?;
        }
        write!(out, "\tAS:i:{}", result.score)?;
        write!(out, "\tXS:i:{}", result.xs)?;
    }
    writeln!(out)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::paired::MateFields;
    use crate::types::{Cigar, CigarOp};

    #[test]
    fn oriented_seq_qual_forward() {
        // bases [0,1,2,2] = ACGG, qual "ABCD" → forward: seq "ACGG", qual "ABCD"
        let (seq, qual) = oriented_seq_qual(&[0u8, 1, 2, 2], "ABCD", false);
        assert_eq!(seq, "ACGG");
        assert_eq!(qual, "ABCD");
    }

    #[test]
    fn oriented_seq_qual_reverse() {
        // bases [0,1,2,2] = ACGG, reverse_complement = CCGT → seq "CCGT", qual reversed "DCBA"
        let (seq, qual) = oriented_seq_qual(&[0u8, 1, 2, 2], "ABCD", true);
        assert_eq!(seq, "CCGT");
        assert_eq!(qual, "DCBA");
    }

    #[test]
    fn oriented_seq_qual_empty_qual() {
        // Empty qual string → QUAL = "*"
        let (seq, qual) = oriented_seq_qual(&[0u8, 1, 2, 3], "", false);
        assert_eq!(seq, "ACGT");
        assert_eq!(qual, "*");
    }

    #[test]
    fn oriented_seq_qual_empty_bases() {
        // Empty bases → SEQ = "*"
        let (seq, _qual) = oriented_seq_qual(&[], "ABCD", false);
        assert_eq!(seq, "*");
    }

    #[test]
    fn write_paired_record_resolves_contig_name_and_pos() {
        // Two contigs: c1=ACGT (len 4), c2=GGGGG (len 5)
        // A read mapped at global position 5 falls in c2 at offset 1 → POS = 2 (1-based).
        // Mate at global pnext=1-based 6 → 0-based 5 → c2 offset 1 → same contig → RNEXT "=".
        let reference = Reference::parse_fasta(">c1\nACGT\n>c2\nGGGGG").unwrap();

        let mut cigar = Cigar::new();
        cigar.push(CigarOp::M, 3);
        let mut result = crate::types::AlignmentResult::new(5, cigar);
        result.mapq = 30;

        let mf = MateFields {
            flag: 0x1 | 0x20 | 0x40,
            rnext: "=",
            pnext: 6,
            tlen: 10,
            placed_pos: None,
            mc: Some("3M".to_string()),
        };

        let mut buf: Vec<u8> = Vec::new();
        write_paired_record(
            &mut buf,
            &reference,
            "read1",
            &result,
            &[0u8, 1, 2, 3],
            "IIII",
            &mf,
        )
        .unwrap();
        let line = String::from_utf8(buf).unwrap();

        let fields: Vec<&str> = line.trim().split('\t').collect();
        assert_eq!(fields[2], "c2", "RNAME should be second contig name");
        assert_eq!(
            fields[3], "2",
            "POS should be per-contig 1-based (offset 1 + 1 = 2)"
        );
        assert_eq!(
            fields[6], "=",
            "RNEXT should be '=' when mate on same contig"
        );
    }

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
            "read1".to_string(),
            0,
            "chr1".to_string(),
            100,
            60,
            "10M".to_string(),
            "*".to_string(),
            0,
            0,
            "ACGTACGT".to_string(),
            "IIIIIIII".to_string(),
        );
        let line = record.to_string();
        assert!(line.starts_with("read1"));
        assert!(line.contains("chr1"));
    }

    #[test]
    fn test_sam_record_fields() {
        let record = SAMRecord::new(
            "read1".to_string(),
            99,
            "chr1".to_string(),
            1000,
            30,
            "50M".to_string(),
            "=".to_string(),
            1050,
            100,
            "ACGT".to_string(),
            "!@#$".to_string(),
        );
        assert_eq!(record.flag, 99);
        assert_eq!(record.pos, 1000);
        assert_eq!(record.mapq, 30);
        assert_eq!(record.cigar, "50M");
    }

    #[test]
    fn write_paired_record_emits_nm_md_mc_as_in_order() {
        use crate::reference::Reference;
        use crate::types::CigarOp;

        let reference = Reference::parse_fasta(&format!(">ref\n{}", "A".repeat(210))).unwrap();

        let mut cigar = Cigar::new();
        cigar.push(CigarOp::M, 100);
        let mut r1 = crate::types::AlignmentResult::new(0, cigar.clone());
        r1.nm = 2;
        r1.score = 151;
        r1.md_tag = Some("MD:Z:50A49".to_string());

        let mut r2 = crate::types::AlignmentResult::new(200, cigar);
        r2.reverse_strand = true;

        let mf = MateFields {
            flag: 0x1 | 0x20 | 0x40,
            rnext: "=",
            pnext: 201,
            tlen: 300,
            placed_pos: None,
            mc: Some("100M".to_string()),
        };

        let mut buf: Vec<u8> = Vec::new();
        write_paired_record(
            &mut buf,
            &reference,
            "read1",
            &r1,
            &[0u8; 100],
            &"I".repeat(100),
            &mf,
        )
        .unwrap();
        let line = String::from_utf8(buf).unwrap();

        assert!(
            line.contains("\tNM:i:2\t"),
            "NM tag present with correct value"
        );
        assert!(line.contains("\tMD:Z:50A49\t"), "MD tag present");
        assert!(line.contains("\tMC:Z:100M\t"), "MC tag present");
        assert!(line.contains("\tAS:i:151"), "AS tag present");
        assert!(line.contains("\tXS:i:0"), "XS:i:0 emitted on mapped record");

        let nm_pos = line.find("\tNM:i:").unwrap();
        let md_pos = line.find("\tMD:Z:").unwrap();
        let mc_pos = line.find("\tMC:Z:").unwrap();
        let as_pos = line.find("\tAS:i:").unwrap();
        assert!(nm_pos < md_pos, "NM before MD");
        assert!(md_pos < mc_pos, "MD before MC");
        assert!(mc_pos < as_pos, "MC before AS");
    }

    #[test]
    fn write_paired_record_emits_xs_after_as_when_xs_nonzero() {
        use crate::reference::Reference;
        use crate::types::CigarOp;

        let reference = Reference::parse_fasta(&format!(">ref\n{}", "A".repeat(210))).unwrap();

        let mut cigar = Cigar::new();
        cigar.push(CigarOp::M, 100);
        let mut r1 = crate::types::AlignmentResult::new(0, cigar.clone());
        r1.nm = 0;
        r1.score = 100;
        r1.xs = 80;
        r1.md_tag = Some("MD:Z:100".to_string());

        let mf = MateFields {
            flag: 0x1 | 0x20 | 0x40,
            rnext: "=",
            pnext: 201,
            tlen: 300,
            placed_pos: None,
            mc: Some("100M".to_string()),
        };

        let mut buf: Vec<u8> = Vec::new();
        write_paired_record(
            &mut buf,
            &reference,
            "repeat_read",
            &r1,
            &[0u8; 100],
            &"I".repeat(100),
            &mf,
        )
        .unwrap();
        let line = String::from_utf8(buf).unwrap();

        assert!(line.contains("\tXS:i:80"), "XS tag emitted when xs > 0");
        let as_pos = line.find("\tAS:i:").unwrap();
        let xs_pos = line.find("\tXS:i:").unwrap();
        assert!(as_pos < xs_pos, "XS must appear after AS");
    }

    #[test]
    fn write_paired_record_unmapped_emits_no_tags() {
        use crate::reference::Reference;

        let reference = Reference::parse_fasta(&format!(">ref\n{}", "A".repeat(210))).unwrap();

        let mut r = crate::types::AlignmentResult::new(0, Cigar::new());
        r.flag = 0x4;
        r.nm = 0;
        r.md_tag = None;

        let mf = MateFields {
            flag: 0x1 | 0x4 | 0x80,
            rnext: "*",
            pnext: 0,
            tlen: 0,
            placed_pos: None,
            mc: None,
        };

        let mut buf: Vec<u8> = Vec::new();
        write_paired_record(
            &mut buf,
            &reference,
            "read_unmap",
            &r,
            &[0u8, 1, 2, 3],
            "IIII",
            &mf,
        )
        .unwrap();
        let line = String::from_utf8(buf).unwrap();

        let fields: Vec<&str> = line.trim().split('\t').collect();
        assert_eq!(
            fields.len(),
            11,
            "unmapped record should have exactly 11 fields, no tags"
        );
        assert!(!line.contains("NM:"), "no NM on unmapped");
        assert!(!line.contains("MD:"), "no MD on unmapped");
        assert!(!line.contains("MC:"), "no MC on unmapped");
        assert!(!line.contains("AS:"), "no AS on unmapped");
    }
}
