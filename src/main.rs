//! BWA-MEM CLI

use clap::Parser;
use std::path::PathBuf;
use std::io::Write;

use bwa_mem::{Aligner, BwaError, FMIndex, Reference, fastq::FASTQReader};
use bwa_mem::types::AlignmentResult;

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

    /// Minimum seed length (default: 19)
    #[arg(short = 'k', default_value = "19")]
    pub min_seed_len: u32,
}

fn main() -> Result<(), BwaError> {
    let cli = Cli::parse();

    match cli {
        Cli::Index(args) => {
            let reference = Reference::from_fasta(&args.reference)?;
            let _index = FMIndex::build(&reference);
            println!("Indexed {} bases", reference.total_len());
            println!("Index prefix: {:?}", args.prefix);
            Ok(())
        }
        Cli::Mem(args) => run_mem(args),
    }
}

fn run_mem(args: MemArgs) -> Result<(), BwaError> {
    let reference = Reference::from_fasta(&args.reference)?;
    let ref_data = reference.as_slice();
    let index = FMIndex::build(&reference);
    
    let aligner = Aligner::new(index, ref_data.to_vec())
        .min_seed_len(args.min_seed_len as usize);

    let mut output: Box<dyn Write> = if args.output.to_string_lossy() == "-" {
        Box::new(std::io::stdout())
    } else {
        Box::new(std::fs::File::create(&args.output)?)
    };

    write_header(&mut output, &reference)?;

    let read1_reader = FASTQReader::from_path(&args.read1)?;
    let read2_reader = args.read2.as_ref().map(|p| FASTQReader::from_path(p)).transpose()?;

    let mut count = 0u64;
    let mut iter = read1_reader.into_iter();
    
    if let Some(r2) = read2_reader {
        let mut r2_iter = r2.into_iter();
        for r in iter {
            let r1 = r?;
            let r2_record = r2_iter.next().transpose()?;
            
            let seq1 = r1.to_sequence();
            let qual1 = r1.to_quality();
            
            let result1 = aligner.align_read(&seq1.bases, Some(&qual1))?;
            write_sam_record(&mut output, &r1.qname, &result1, false)?;
            
            if let Some(r2) = r2_record {
                let seq2 = r2.to_sequence();
                let qual2 = r2.to_quality();
                let result2 = aligner.align_read(&seq2.bases, Some(&qual2))?;
                write_sam_record(&mut output, &r2.qname, &result2, true)?;
            }
            
            count += 1;
            if count % 10000 == 0 {
                eprintln!("Processed {} read pairs", count);
            }
        }
    } else {
        for record in iter {
            let r1 = record?;
            let seq = r1.to_sequence();
            let qual = r1.to_quality();
            
            let result = aligner.align_read(&seq.bases, Some(&qual))?;
            write_sam_record(&mut output, &r1.qname, &result, false)?;
            count += 1;
            
            if count % 10000 == 0 {
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
        let name = if contig.name.is_empty() { "ref" } else { contig.name.as_str() };
        writeln!(output, "@SQ\tSN:{}\tLN:{}", name, len)?;
    }
    writeln!(output, "@PG\tID:bwa-rs\tPN:bwa-rs\tVN:0.1.0")?;
    Ok(())
}

fn write_sam_record(
    output: &mut Box<dyn Write>,
    qname: &str,
    result: &AlignmentResult,
    is_mate: bool,
) -> Result<(), BwaError> {
    let flag = if result.position == 0 { 0x4 } else { 0 };
    let flag = if is_mate { flag | 0x80 } else { flag };
    let rname = if result.position == 0 { "*" } else { "ref" };
    let pos = if result.position == 0 { 0 } else { result.position as i64 + 1 };
    let mapq = result.mapq;
    let cigar_str = result.cigar.to_string();
    let cigar = if result.cigar.ops.is_empty() { "*" } else { cigar_str.as_str() };
    let rnext = "*";
    let pnext: i64 = 0;
    let tlen: i64 = 0;
    let seq = "*";
    let qual = "*";

    writeln!(output, "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
        qname, flag, rname, pos, mapq, cigar, rnext, pnext, tlen, seq, qual)?;
    Ok(())
}