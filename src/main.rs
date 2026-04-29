//! BWA-MEM CLI

use clap::Parser;
use std::path::PathBuf;

use bwa_mem::{Aligner, BwaError, FMIndex, Reference};

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

    /// Minimum seed length
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
        Cli::Mem(args) => {
            let reference = Reference::from_fasta(&args.reference)?;
            let ref_data = reference.as_slice().to_vec();
            let index = FMIndex::build(&reference);
            let _aligner = Aligner::new(index, ref_data).min_seed_len(args.min_seed_len as usize);

            println!("Aligning reads from {:?}", args.read1);
            println!("Reference: {:?}", args.reference);

            if let Some(r2) = args.read2 {
                println!("Paired-end mode: {:?}", r2);
            }

            Ok(())
        }
    }
}