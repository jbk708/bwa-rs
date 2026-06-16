//! BWA-MEM CLI

use clap::Parser;
use std::io::{BufWriter, Write};
use std::path::PathBuf;

use bwa_mem::paired::{
    mate_fields, mem_pair, pair_insert_size, paired_mapq, InsertSizeDistribution, PEN_UNPAIRED,
};
use bwa_mem::sam::{oriented_seq_qual, write_paired_record as sam_write_paired};
use bwa_mem::types::{AlignmentResult, Cigar};
use bwa_mem::{
    fastq::FASTQReader, BwaError, FMIndex, ParallelAligner, ReadRegions, Reference,
    ThreadPoolConfig,
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

    /// Output index prefix (default: reference path)
    #[arg(short = 'p')]
    pub prefix: Option<PathBuf>,
}

#[derive(Parser)]
struct MemArgs {
    /// Reference genome (indexed)
    #[arg(short = 'R')]
    pub reference: PathBuf,

    /// Index prefix (default: reference path)
    #[arg(short = 'p')]
    pub prefix: Option<PathBuf>,

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

    /// Number of threads (default: auto-detect)
    #[arg(short = 't', default_value = "0")]
    pub threads: u32,

    /// Minimum alignment score to output (default: 30)
    #[arg(short = 'T', default_value = "30")]
    pub min_score: i32,
}

/// Returns `prefix` with `.bwarsidx` appended to its filename component.
fn index_path_for(prefix: &std::path::Path) -> std::path::PathBuf {
    let mut name = prefix
        .file_name()
        .map(std::ffi::OsString::from)
        .unwrap_or_default();
    name.push(".bwarsidx");
    prefix.with_file_name(name)
}

fn main() -> Result<(), BwaError> {
    let cli = Cli::parse();

    match cli {
        Cli::Index(args) => {
            let reference = Reference::from_fasta(&args.reference)?;
            let prefix = args.prefix.unwrap_or_else(|| args.reference.clone());
            let index = FMIndex::build_2n(&reference);
            let path = index_path_for(&prefix);
            index.save(&path)?;
            println!("Indexed {} bases", reference.total_len());
            println!("Wrote index to {:?}", path);
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
    regions: Vec<AlignmentResult>,
    sub: i32,
    sub_n: u32,
    frac_rep: f32,
    csub: i32,
}

struct RawRead {
    qname: String,
    bases: Vec<u8>,
    qual: String,
}

const BATCH_SIZE: usize = 65536;

/// Align a batch of queries in parallel, surfacing the first error.
fn align_all(
    aligner: &ParallelAligner,
    queries: &[&[u8]],
) -> Result<Vec<AlignmentResult>, BwaError> {
    aligner.align_batch(queries).into_iter().collect()
}

fn align_all_regions(
    aligner: &ParallelAligner,
    queries: &[&[u8]],
) -> Result<Vec<ReadRegions>, BwaError> {
    aligner.align_batch_regions(queries).into_iter().collect()
}

fn unmapped_aln() -> AlignmentResult {
    let mut u = AlignmentResult::new(0, Cigar::new());
    u.flag = 0x4;
    u
}

/// Build a paired-end `ReadAln` from a raw read and its candidate regions, retaining
/// only the near-best placements (within `PEN_UNPAIRED` of the primary) so mate-aware
/// pairing considers only genuine competitors. `result` is seeded with the primary
/// (or an unmapped placeholder) and may be replaced later by pairing or mate rescue.
fn paired_read_aln(raw: RawRead, rr: ReadRegions) -> ReadAln {
    let regions: Vec<AlignmentResult> = if rr.regions.is_empty() {
        Vec::new()
    } else {
        let cutoff = rr.regions[0].score - PEN_UNPAIRED;
        rr.regions
            .into_iter()
            .filter(|r| r.score >= cutoff)
            .collect()
    };
    let result = regions.first().cloned().unwrap_or_else(unmapped_aln);
    ReadAln {
        qname: raw.qname,
        bases: raw.bases,
        qual: raw.qual,
        result,
        regions,
        sub: rr.sub,
        sub_n: rr.sub_n,
        frac_rep: rr.frac_rep,
        csub: rr.csub,
    }
}

fn run_mem(args: MemArgs) -> Result<(), BwaError> {
    if args.threads > 0 {
        ThreadPoolConfig::new()
            .num_threads(args.threads as usize)
            .apply();
    }

    let reference = Reference::from_fasta(&args.reference)?;
    let prefix = args
        .prefix
        .clone()
        .unwrap_or_else(|| args.reference.clone());
    let index_path = index_path_for(&prefix);
    let index = if index_path.exists() {
        let index = FMIndex::load(&index_path)?;
        if index.n_fwd != reference.total_len() {
            return Err(BwaError::Index(format!(
                "index {:?} does not match reference (n_fwd {} != reference length {})",
                index_path,
                index.n_fwd,
                reference.total_len()
            )));
        }
        eprintln!("Loaded index from {:?}", index_path);
        index
    } else {
        eprintln!(
            "Index not found at {:?}; building in memory (run `bwa-mem index` to persist)",
            index_path
        );
        FMIndex::build_2n(&reference)
    };
    let ref_data = reference.as_slice_2n();

    let aligner = ParallelAligner::new(index, ref_data.to_vec())
        .min_seed_len(args.min_seed_len as usize)
        .min_score(args.min_score);

    let mut output: Box<dyn Write> = if args.output.to_string_lossy() == "-" {
        Box::new(BufWriter::new(std::io::stdout()))
    } else {
        Box::new(BufWriter::new(std::fs::File::create(&args.output)?))
    };

    let cmdline = std::env::args().collect::<Vec<_>>().join(" ");
    write_header(&mut output, &reference, &cmdline)?;

    let read1_reader = FASTQReader::from_path(&args.read1)?;
    let read2_reader = args
        .read2
        .as_ref()
        .map(|p| FASTQReader::from_path(p))
        .transpose()?;

    let mut count = 0u64;
    let iter = read1_reader.into_iter();

    if let Some(r2) = read2_reader {
        // Two-phase: align every pair and estimate the insert-size distribution
        // first, then emit records so proper-pair flags reflect the full dataset.
        //
        // The whole file must be buffered here (it cannot be chunked): mem_pair and
        // paired_mapq below need the complete InsertSizeDistribution, which is only
        // final once every pair has been aligned.
        let mut r2_iter = r2.into_iter();

        // Step 1: read all pairs into owned storage, preserving input order.
        let mut raw_pairs: Vec<(RawRead, Option<RawRead>)> = Vec::new();
        for r in iter {
            let r1 = r?;
            let seq1 = r1.to_sequence();
            let read1 = RawRead {
                qname: r1.qname,
                bases: seq1.bases,
                qual: r1.qual,
            };
            let read2 = r2_iter.next().transpose()?.map(|r2_rec| {
                let seq2 = r2_rec.to_sequence();
                RawRead {
                    qname: r2_rec.qname,
                    bases: seq2.bases,
                    qual: r2_rec.qual,
                }
            });
            raw_pairs.push((read1, read2));
        }

        // Step 2: flat ordered query list — read1 then read2 per pair. align_batch_regions
        // preserves input order, so results map 1-to-1 back onto the pairs.
        let mut query_refs: Vec<&[u8]> = Vec::with_capacity(raw_pairs.len() * 2);
        for (r1, r2) in &raw_pairs {
            query_refs.push(r1.bases.as_slice());
            if let Some(r2) = r2 {
                query_refs.push(r2.bases.as_slice());
            }
        }

        // Step 3: align the whole batch in parallel, collecting full region arrays.
        let mut result_iter = align_all_regions(&aligner, &query_refs)?.into_iter();

        // Step 4: walk results in order to build buf and collect FR insert sizes from primaries.
        let mut buf: Vec<(ReadAln, Option<ReadAln>)> = Vec::with_capacity(raw_pairs.len());
        let mut insert_sizes: Vec<u32> = Vec::with_capacity(raw_pairs.len());
        let mut next_rr = || {
            result_iter
                .next()
                .ok_or_else(|| BwaError::Alignment("batch result count mismatch".into()))
        };
        for (r1, r2) in raw_pairs {
            let read1 = paired_read_aln(r1, next_rr()?);

            let read2 = if let Some(r2) = r2 {
                let read2 = paired_read_aln(r2, next_rr()?);
                if let Some(isize) = pair_insert_size(&read1.result, &read2.result) {
                    insert_sizes.push(isize);
                }
                Some(read2)
            } else {
                None
            };

            buf.push((read1, read2));

            count += 1;
            if count.is_multiple_of(10000) {
                eprintln!("Processed {} read pairs", count);
            }
        }

        let dist = InsertSizeDistribution::from_insert_sizes(&insert_sizes);

        for (mut read1, pair) in buf {
            if let Some(mut read2) = pair {
                // Mate rescue: when exactly one mate mapped, try to place the other by
                // banded SW within the insert-size window around the mapped mate.
                if read1.result.is_unmapped() && !read2.result.is_unmapped() {
                    if let Some(rescued) = aligner.rescue_mate(&read1.bases, &read2.result, &dist) {
                        read1.result = rescued.clone();
                        read1.regions = vec![rescued];
                    }
                } else if read2.result.is_unmapped() && !read1.result.is_unmapped() {
                    if let Some(rescued) = aligner.rescue_mate(&read2.bases, &read1.result, &dist) {
                        read2.result = rescued.clone();
                        read2.regions = vec![rescued];
                    }
                }

                let mut proper = false;
                if !read1.regions.is_empty() && !read2.regions.is_empty() {
                    if let Some(choice) = mem_pair(
                        &read1.regions,
                        &read2.regions,
                        &dist,
                        aligner.match_score(),
                        aligner.pair_margin(),
                    ) {
                        let score_un =
                            read1.regions[0].score + read2.regions[0].score - PEN_UNPAIRED;
                        if choice.score > score_un {
                            let mut place1 = read1.regions[choice.idx1].clone();
                            let mut place2 = read2.regions[choice.idx2].clone();
                            let q_se1 = aligner.region_se_mapq(
                                &place1,
                                read1.sub,
                                read1.sub_n,
                                read1.frac_rep,
                                read1.csub,
                            );
                            let q_se2 = aligner.region_se_mapq(
                                &place2,
                                read2.sub,
                                read2.sub_n,
                                read2.frac_rep,
                                read2.csub,
                            );
                            let (m1, m2) = paired_mapq(
                                &choice,
                                score_un,
                                q_se1,
                                q_se2,
                                read1.frac_rep,
                                read2.frac_rep,
                                aligner.match_score(),
                                read1.csub,
                                read2.csub,
                                place1.score,
                                place2.score,
                            );
                            place1.mapq = m1;
                            place2.mapq = m2;
                            read1.result = place1;
                            read2.result = place2;
                            proper = true;
                        }
                    }
                }

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
        // Chunked batch alignment keeps memory bounded while still parallelising
        // within each chunk; align_batch preserves input order so output is identical.
        let mut chunk: Vec<RawRead> = Vec::with_capacity(BATCH_SIZE);

        let flush_chunk = |chunk: &mut Vec<RawRead>,
                           output: &mut Box<dyn Write>,
                           count: &mut u64|
         -> Result<(), BwaError> {
            let query_refs: Vec<&[u8]> = chunk.iter().map(|r| r.bases.as_slice()).collect();
            let results = align_all(&aligner, &query_refs)?;
            for (r, result) in chunk.iter().zip(results.iter()) {
                write_sam_record(
                    output, &reference, &r.qname, result, &r.bases, &r.qual, false,
                )?;
                *count += 1;
                if count.is_multiple_of(10000) {
                    eprintln!("Processed {} reads", count);
                }
            }
            chunk.clear();
            Ok(())
        };

        for record in iter {
            let r1 = record?;
            let seq = r1.to_sequence();
            chunk.push(RawRead {
                qname: r1.qname,
                bases: seq.bases,
                qual: r1.qual,
            });
            if chunk.len() >= BATCH_SIZE {
                flush_chunk(&mut chunk, &mut output, &mut count)?;
            }
        }
        if !chunk.is_empty() {
            flush_chunk(&mut chunk, &mut output, &mut count)?;
        }
    }

    eprintln!("Aligned {} reads", count);
    Ok(())
}

fn header_text(reference: &Reference, cmdline: &str) -> String {
    let mut text = String::new();
    for contig in &reference.contigs {
        let name = if contig.name.is_empty() {
            "ref"
        } else {
            contig.name.as_str()
        };
        text.push_str(&format!("@SQ\tSN:{}\tLN:{}\n", name, contig.len()));
    }
    text.push_str(&format!(
        "@PG\tID:bwa-rs\tPN:bwa-rs\tVN:{}\tCL:{}\n",
        env!("CARGO_PKG_VERSION"),
        cmdline
    ));
    text
}

fn write_header(
    output: &mut Box<dyn Write>,
    reference: &Reference,
    cmdline: &str,
) -> Result<(), BwaError> {
    output.write_all(header_text(reference, cmdline).as_bytes())?;
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
        write!(output, "\tXS:i:{}", result.xs)?;
    }
    writeln!(output)?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use bwa_mem::types::{Cigar, CigarOp};
    use bwa_mem::Reference;

    #[test]
    fn index_path_for_appends_extension_to_path_with_existing_extension() {
        let p = std::path::Path::new("test-data/chr1_bm2.fasta");
        assert_eq!(
            index_path_for(p),
            std::path::PathBuf::from("test-data/chr1_bm2.fasta.bwarsidx")
        );
    }

    #[test]
    fn index_path_for_appends_extension_to_path_without_extension() {
        let p = std::path::Path::new("some/dir/myref");
        assert_eq!(
            index_path_for(p),
            std::path::PathBuf::from("some/dir/myref.bwarsidx")
        );
    }

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
    fn single_end_mapped_record_emits_nm_md_as_xs_no_mc() {
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
        assert!(line.contains("\tXS:i:0"), "XS:i:0 emitted on mapped record");

        let nm_pos = line.find("\tNM:i:").unwrap();
        let md_pos = line.find("\tMD:Z:").unwrap();
        let as_pos = line.find("\tAS:i:").unwrap();
        assert!(nm_pos < md_pos, "NM before MD");
        assert!(md_pos < as_pos, "MD before AS");
    }

    #[test]
    fn single_end_mapped_record_with_xs_emits_xs_after_as() {
        let reference = Reference::parse_fasta(&format!(">ref\n{}", "A".repeat(210))).unwrap();
        let mut cigar = Cigar::new();
        cigar.push(CigarOp::M, 50);
        let mut result = AlignmentResult::new(0, cigar);
        result.nm = 0;
        result.score = 50;
        result.xs = 30;
        result.md_tag = Some("MD:Z:50".to_string());

        let line = call_write_sam_record(
            &reference,
            "repeat_read",
            &result,
            &[0u8; 50],
            &"I".repeat(50),
        );

        assert!(line.contains("\tXS:i:30"), "XS tag emitted when xs > 0");
        let as_pos = line.find("\tAS:i:").unwrap();
        let xs_pos = line.find("\tXS:i:").unwrap();
        assert!(as_pos < xs_pos, "XS must appear after AS");
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

    #[test]
    fn header_omits_hd_and_emits_pg_with_cmdline() {
        let reference = Reference::parse_fasta(&format!(">1\n{}", "A".repeat(100))).unwrap();
        let header = header_text(&reference, "bwa-rs mem -R ref.fa -1 r1.fq -2 r2.fq");

        assert!(!header.contains("@HD"), "bwa-mem2 emits no @HD line");
        assert!(
            header.contains("@SQ\tSN:1\tLN:100\n"),
            "contig @SQ line present"
        );
        let pg = header
            .lines()
            .find(|l| l.starts_with("@PG"))
            .expect("@PG line present");
        assert!(
            pg.starts_with("@PG\tID:bwa-rs\tPN:bwa-rs\tVN:"),
            "bwa-style @PG"
        );
        assert!(
            pg.contains("\tCL:bwa-rs mem -R ref.fa -1 r1.fq -2 r2.fq"),
            "@PG records the command line"
        );
        assert!(
            pg.contains(&format!("\tVN:{}\t", env!("CARGO_PKG_VERSION"))),
            "@PG version comes from crate version"
        );
    }
}
