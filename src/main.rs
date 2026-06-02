//! BWA-MEM CLI

use clap::Parser;
use std::io::{BufWriter, Write};
use std::path::PathBuf;

use bwa_mem::paired::{is_proper_pair, mate_fields, pair_insert_size, InsertSizeDistribution};
use bwa_mem::sam::{oriented_seq_qual, write_paired_record as sam_write_paired};
use bwa_mem::types::AlignmentResult;
use bwa_mem::{
    fastq::FASTQReader, BwaError, FMIndex, ParallelAligner, Reference, ThreadPoolConfig,
};

#[derive(Parser)]
#[command(name = "bwa-rs")]
#[command(about = "Pure Rust BWA-MEM implementation")]
enum Cli {
    /// Build FM-index for reference
    Index(IndexArgs),
    /// Align reads to reference
    Mem(MemArgs),
}

#[derive(Parser)]
struct IndexArgs {
    /// Reference FASTA file
    #[arg(short = 'r')]
    pub reference: PathBuf,

    /// Output index prefix
    #[arg(short = 'p', default_value = "index")]
    pub prefix: PathBuf,
}

#[derive(Parser)]
struct MemArgs {
    /// Reference genome (indexed)
    #[arg(short = 'R')]
    pub reference: PathBuf,

    /// FASTQ file for read 1
    #[arg(short = '1')]
    pub read1: PathBuf,

    /// FASTQ file for read 2
    #[arg(short = '2')]
    pub read2: Option<PathBuf>,

    /// Output file (SAM format)
    #[arg(short = 'o', default_value = "-")]
    pub output: PathBuf,

    /// Minimum seed length (default: 10)
    #[arg(short = 'k', default_value = "10")]
    pub min_seed_len: u32,

    /// Number of threads (default: auto-detect)
    #[arg(short = 't', default_value = "0")]
    pub threads: u32,

    /// Minimum alignment score to output (default: 30)
    #[arg(short = 'T', default_value = "30")]
    pub min_score: i32,
}

fn main() -> Result<(), BwaError> {
    let cli = Cli::parse();

    match cli {
        Cli::Index(args) => {
            let reference = Reference::from_fasta(&args.reference)?;
            // 2N index (forward + reverse-complement), matching bwa / bwa-mem2.
            let _index = FMIndex::build_2n(&reference);
            println!("Indexed {} bases", reference.total_len());
            println!("Index prefix: {:?}", args.prefix);
            Ok(())
        }
        Cli::Mem(args) => run_mem(args),
    }
}

struct ReadAln {
    qname: String,
    bases: Vec<u8>,
    qual: String,
    result: AlignmentResult,
}

fn run_mem(args: MemArgs) -> Result<(), BwaError> {
    if args.threads > 0 {
        ThreadPoolConfig::new()
            .num_threads(args.threads as usize)
            .apply();
    }

    let reference = Reference::from_fasta(&args.reference)?;
    // Build the 2N index and extend against the 2N sequence so the aligner can
    // place reads on both strands (positions in the reverse-complement half are
    // mapped back to forward coordinates with the 0x10 flag set).
    let ref_data = reference.as_slice_2n();
    let index = FMIndex::build_2n(&reference);

    let aligner = ParallelAligner::new(index, ref_data.to_vec())
        .min_seed_len(args.min_seed_len as usize)
        .min_score(args.min_score);

    let mut output: Box<dyn Write> = if args.output.to_string_lossy() == "-" {
        Box::new(BufWriter::new(std::io::stdout()))
    } else {
        Box::new(BufWriter::new(std::fs::File::create(&args.output)?))
    };

    write_header(&mut output, &reference)?;

    let read1_reader = FASTQReader::from_path(&args.read1)?;
    let read2_reader = args
        .read2
        .as_ref()
        .map(|p| FASTQReader::from_path(p))
        .transpose()?;

    let mut count = 0u64;
    let iter = read1_reader.into_iter();

    if let Some(r2) = read2_reader {
        // Two-phase: align all pairs and estimate the insert-size distribution first,
        // then emit records so proper-pair flags reflect the full dataset.
        let mut r2_iter = r2.into_iter();
        let mut buf: Vec<(ReadAln, Option<ReadAln>)> = Vec::new();
        let mut dist = InsertSizeDistribution::new();

        for r in iter {
            let r1 = r?;
            let r2_record = r2_iter.next().transpose()?;
            let seq1 = r1.to_sequence();
            let result1 = aligner.align_single(&seq1.bases)?;
            let read1 = ReadAln {
                qname: r1.qname,
                bases: seq1.bases,
                qual: r1.qual,
                result: result1,
            };

            let pair = if let Some(r2_rec) = r2_record {
                let seq2 = r2_rec.to_sequence();
                let result2 = aligner.align_single(&seq2.bases)?;
                if let Some(isize) = pair_insert_size(&read1.result, &result2) {
                    dist.add(isize);
                }
                Some(ReadAln {
                    qname: r2_rec.qname,
                    bases: seq2.bases,
                    qual: r2_rec.qual,
                    result: result2,
                })
            } else {
                None
            };

            buf.push((read1, pair));

            count += 1;
            if count.is_multiple_of(10000) {
                eprintln!("Processed {} read pairs", count);
            }
        }

        for (read1, pair) in buf {
            if let Some(read2) = pair {
                let proper = is_proper_pair(&read1.result, &read2.result, &dist);
                let mf1 = mate_fields(&read1.result, &read2.result, true, proper);
                let mf2 = mate_fields(&read2.result, &read1.result, false, proper);
                let mut emit = |r: &ReadAln, mf: &_| {
                    sam_write_paired(
                        &mut output,
                        &reference,
                        &r.qname,
                        &r.result,
                        &r.bases,
                        &r.qual,
                        mf,
                    )
                };
                emit(&read1, &mf1)?;
                emit(&read2, &mf2)?;
            } else {
                write_sam_record(
                    &mut output,
                    &reference,
                    &read1.qname,
                    &read1.result,
                    &read1.bases,
                    &read1.qual,
                    false,
                )?;
            }
        }
    } else {
        for record in iter {
            let r1 = record?;
            let seq = r1.to_sequence();

            let result = aligner.align_single(&seq.bases)?;
            write_sam_record(
                &mut output,
                &reference,
                &r1.qname,
                &result,
                &seq.bases,
                &r1.qual,
                false,
            )?;
            count += 1;

            if count.is_multiple_of(10000) {
                eprintln!("Processed {} reads", count);
            }
        }
    }

    eprintln!("Aligned {} reads", count);
    Ok(())
}

fn write_header(output: &mut Box<dyn Write>, reference: &Reference) -> Result<(), BwaError> {
    writeln!(output, "@HD\tVN:1.0\tSO:unsorted")?;
    for contig in &reference.contigs {
        let len = contig.len();
        let name = if contig.name.is_empty() {
            "ref"
        } else {
            contig.name.as_str()
        };
        writeln!(output, "@SQ\tSN:{}\tLN:{}", name, len)?;
    }
    writeln!(output, "@PG\tID:bwa-rs\tPN:bwa-rs\tVN:0.1.0")?;
    Ok(())
}

fn write_sam_record(
    output: &mut Box<dyn Write>,
    reference: &Reference,
    qname: &str,
    result: &AlignmentResult,
    bases: &[u8],
    qual: &str,
    is_mate: bool,
) -> Result<(), BwaError> {
    let mapped = !result.cigar.ops.is_empty() && !result.is_unmapped();
    let mut flag = if mapped { 0 } else { 0x4 };
    if mapped && result.reverse_strand {
        flag |= 0x10; // SEQ is reverse-complemented relative to the reference
    }
    let flag = if is_mate { flag | 0x80 } else { flag };
    let (rname, pos) = if mapped {
        reference
            .locate(result.position)
            .map_or(("*", 0), |(name, off)| (name, off as i64 + 1))
    } else {
        ("*", 0)
    };
    let mapq = if mapped { result.mapq } else { 0 };
    let cigar_str = result.cigar.to_sam_string();
    let cigar = if mapped { cigar_str.as_str() } else { "*" };
    let rnext = "*";
    let pnext: i64 = 0;
    let tlen: i64 = 0;
    let (seq, qual_out) = oriented_seq_qual(bases, qual, result.reverse_strand);

    write!(
        output,
        "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
        qname, flag, rname, pos, mapq, cigar, rnext, pnext, tlen, seq, qual_out
    )?;
    if mapped {
        write!(output, "\tNM:i:{}", result.nm)?;
        if let Some(ref md) = result.md_tag {
            write!(output, "\t{}", md)?;
        }
        write!(output, "\tAS:i:{}", result.score)?;
    }
    writeln!(output)?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use bwa_mem::types::{Cigar, CigarOp};
    use bwa_mem::Reference;

    fn call_write_sam_record(
        reference: &Reference,
        qname: &str,
        result: &AlignmentResult,
        bases: &[u8],
        qual: &str,
    ) -> String {
        use std::sync::{Arc, Mutex};
        let shared: Arc<Mutex<Vec<u8>>> = Arc::new(Mutex::new(Vec::new()));
        let shared2 = Arc::clone(&shared);
        struct SharedWriter(Arc<Mutex<Vec<u8>>>);
        impl Write for SharedWriter {
            fn write(&mut self, buf: &[u8]) -> std::io::Result<usize> {
                self.0.lock().unwrap().write(buf)
            }
            fn flush(&mut self) -> std::io::Result<()> {
                Ok(())
            }
        }
        let mut output: Box<dyn Write> = Box::new(SharedWriter(shared2));
        write_sam_record(&mut output, reference, qname, result, bases, qual, false).unwrap();
        drop(output);
        let bytes = shared.lock().unwrap().clone();
        String::from_utf8(bytes).unwrap()
    }

    #[test]
    fn single_end_mapped_record_emits_nm_md_as_no_mc_no_xs() {
        let reference = Reference::parse_fasta(&format!(">ref\n{}", "A".repeat(210))).unwrap();
        let mut cigar = Cigar::new();
        cigar.push(CigarOp::M, 50);
        let mut result = AlignmentResult::new(0, cigar);
        result.nm = 3;
        result.score = 50;
        result.md_tag = Some("MD:Z:50".to_string());

        let line = call_write_sam_record(
            &reference,
            "single_read",
            &result,
            &[0u8; 50],
            &"I".repeat(50),
        );

        assert!(line.contains("\tNM:i:3\t"), "NM tag present");
        assert!(line.contains("\tMD:Z:50\t"), "MD tag present");
        assert!(line.contains("\tAS:i:50"), "AS tag present");
        assert!(!line.contains("MC:"), "no MC in single-end");
        assert!(!line.contains("XS:"), "no XS emitted");

        let nm_pos = line.find("\tNM:i:").unwrap();
        let md_pos = line.find("\tMD:Z:").unwrap();
        let as_pos = line.find("\tAS:i:").unwrap();
        assert!(nm_pos < md_pos, "NM before MD");
        assert!(md_pos < as_pos, "MD before AS");
    }

    #[test]
    fn single_end_unmapped_record_emits_no_tags() {
        let reference = Reference::parse_fasta(&format!(">ref\n{}", "A".repeat(210))).unwrap();
        let mut result = AlignmentResult::new(0, Cigar::new());
        result.flag = 0x4;

        let line = call_write_sam_record(
            &reference,
            "unmapped_read",
            &result,
            &[0u8, 1, 2, 3],
            "IIII",
        );

        let fields: Vec<&str> = line.trim().split('\t').collect();
        assert_eq!(fields.len(), 11, "unmapped single-end: exactly 11 fields");
        assert!(!line.contains("NM:"), "no NM on unmapped");
        assert!(!line.contains("MD:"), "no MD on unmapped");
        assert!(!line.contains("AS:"), "no AS on unmapped");
    }
}
