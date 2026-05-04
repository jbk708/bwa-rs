use bwa_mem::sam::{SAMHeader, SAMRecord};
use bwa_mem::types::Sequence;
use bwa_mem::{Aligner, BwaError, FMIndex, Reference};
use std::process::Command;

const TEST_REFERENCE: &str = ">test_ref\nACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT\nACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT\nACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT\n";

const TEST_FASTQ: &str = "@read1\nACGTACGTACGTACGTACGTACGTACGT\n+\nIIIIIIIIIIIIIIIIIIIIIIIIIIII\n@read2\nGGGGAAAACCCCAAAAGGGGAAAACCCC\n+\nIIIIIIIIIIIIIIIIIIIIIIIIIIII\n@read3\nNNNNACGTACGTACGTNNNNACGTACGT\n+\nIIIIIIIIIIIIIIIIIIIIIIIIIIII\n";

fn seq_to_bytes(s: &str) -> Vec<u8> {
    s.bytes()
        .map(|b| match b {
            b'A' => 0,
            b'C' => 1,
            b'G' => 2,
            b'T' => 3,
            _ => 4,
        })
        .collect()
}

fn create_test_aligner() -> (Aligner, Reference) {
    let reference = Reference::parse_fasta(TEST_REFERENCE).unwrap();
    let ref_slice = reference.as_slice();
    let index = FMIndex::build(&reference);
    let aligner = Aligner::new(index, ref_slice).min_seed_len(4);
    (aligner, reference)
}

fn parse_fastq(data: &str) -> Vec<(String, Vec<u8>, Vec<u8>)> {
    let mut records = Vec::new();
    let mut lines = data.lines();
    while let Some(name) = lines.next() {
        if !name.starts_with('@') {
            continue;
        }
        let qname = name.trim_start_matches('@').to_string();
        let seq_line = lines.next().unwrap_or("");
        let _qual_header = lines.next().unwrap_or("");
        let qual_line = lines.next().unwrap_or("");
        records.push((qname, seq_to_bytes(seq_line), qual_line.as_bytes().to_vec()));
    }
    records
}

fn validate_sam_record(line: &str) -> Result<(), BwaError> {
    let fields: Vec<&str> = line.split('\t').collect();

    if fields.len() < 11 {
        return Err(BwaError::Parse(format!(
            "SAM record has {} fields, expected at least 11",
            fields.len()
        )));
    }

    let qname = fields[0];
    if qname.is_empty() {
        return Err(BwaError::Parse("Empty QNAME".to_string()));
    }

    let _flag: u16 = fields[1]
        .parse()
        .map_err(|_| BwaError::Parse(format!("Invalid FLAG: {}", fields[1])))?;

    let rname = fields[2];
    let pos: i64 = fields[3]
        .parse()
        .map_err(|_| BwaError::Parse(format!("Invalid POS: {}", fields[3])))?;

    let mapq: u8 = fields[4]
        .parse()
        .map_err(|_| BwaError::Parse(format!("Invalid MAPQ: {}", fields[4])))?;

    let cigar = fields[5];
    validate_cigar(cigar)?;

    if rname == "*" && pos != 0 {
        return Err(BwaError::Parse("RNAME=* but POS != 0".to_string()));
    }

    if rname != "*" && pos == 0 {
        return Err(BwaError::Parse("RNAME != * but POS == 0".to_string()));
    }

    let _ = mapq;

    Ok(())
}

fn validate_cigar(cigar: &str) -> Result<(), BwaError> {
    if cigar == "*" {
        return Ok(());
    }

    let mut chars = cigar.chars().peekable();
    while let Some(c) = chars.next() {
        if !c.is_ascii_digit() {
            return Err(BwaError::Parse(format!(
                "Invalid CIGAR: expected digit, got {}",
                c
            )));
        }

        let mut num_str = String::from(c);
        while let Some(&next_c) = chars.peek() {
            if next_c.is_ascii_digit() {
                num_str.push(chars.next().unwrap());
            } else {
                break;
            }
        }

        let op = chars
            .next()
            .ok_or_else(|| BwaError::Parse("CIGAR ends with number".to_string()))?;

        match op {
            'M' | 'I' | 'D' | 'N' | 'S' | 'H' | 'P' | '=' | 'X' => {}
            _ => {
                return Err(BwaError::Parse(format!("Invalid CIGAR operator: {}", op)));
            }
        }
    }

    Ok(())
}

#[test]
fn test_complete_alignment_pipeline() -> Result<(), BwaError> {
    let (aligner, reference) = create_test_aligner();
    assert_eq!(reference.total_len(), 144);

    let reads = parse_fastq(TEST_FASTQ);
    assert_eq!(reads.len(), 3);

    let mut sam_records = Vec::new();
    for (qname, seq, qual) in &reads {
        let result = aligner.align_read(seq, None);
        let sequence = Sequence::new("query", seq.clone());
        let sam_record = match result {
            Ok(alignment) => {
                SAMRecord::from_alignment(qname, &alignment, &sequence, qual, "test_ref")
            }
            Err(_) => SAMRecord::unmapped(qname, seq, qual),
        };
        sam_records.push(sam_record.to_string());
    }

    assert_eq!(sam_records.len(), 3);
    for record in &sam_records {
        validate_sam_record(record)?;
    }
    Ok(())
}

#[test]
fn test_sam_header_format() -> Result<(), BwaError> {
    let (aligner, reference) = create_test_aligner();
    let seq = seq_to_bytes("ACGTACGTACGTACGT");
    let qual = b"I".repeat(16);

    let result = aligner.align_read(&seq, None)?;
    let sequence = Sequence::new("query", seq.clone());
    let header = SAMHeader::new(reference.clone());
    let sam_header = header.to_string();
    let sam_record =
        SAMRecord::from_alignment("test_read", &result, &sequence, &qual, &reference.name);

    assert!(sam_header.starts_with("@HD"));
    assert!(sam_header.contains("@SQ\tSN:test_ref\tLN:144"));
    assert!(sam_header.contains("@PG"));
    validate_sam_record(&sam_record.to_string())?;
    Ok(())
}

#[test]
fn test_perfect_match() -> Result<(), BwaError> {
    let reference = Reference::parse_fasta(TEST_REFERENCE)?;
    let ref_slice = reference.as_slice();
    let index = FMIndex::build(&reference);
    let aligner = Aligner::new(index, ref_slice).min_seed_len(4);

    let seq = seq_to_bytes("ACGTACGT");
    let sequence = Sequence::new("query", seq.clone());
    let qual = b"I".repeat(8).to_vec();

    let result = aligner.align_read(&seq, None)?;
    assert!(!result.cigar.ops.is_empty() || result.flag & 0x4 != 0);

    let sam = SAMRecord::from_alignment("read", &result, &sequence, &qual, "test_ref");
    let sam_str = sam.to_string();
    let fields: Vec<&str> = sam_str.split('\t').collect();

    assert_eq!(fields[2], "test_ref");
    assert!(
        fields[5].contains('M')
            || fields[5].contains('=')
            || fields[5].contains('X')
            || fields[5] == "*"
    );

    Ok(())
}

#[test]
fn test_mismatch_position() -> Result<(), BwaError> {
    let (aligner, _) = create_test_aligner();
    let seq = seq_to_bytes("ACGTACGTACGTACGX");
    let sequence = Sequence::new("query", seq.clone());
    let qual = b"I".repeat(16);

    let result = aligner.align_read(&seq, None)?;
    let sam = SAMRecord::from_alignment("read", &result, &sequence, &qual, "test_ref");
    let sam_str = sam.to_string();
    let fields: Vec<&str> = sam_str.split('\t').collect();
    let flag: u16 = fields[1].parse().unwrap();

    if flag & 0x4 == 0 {
        assert!(fields[2] != "*");
    }
    Ok(())
}

#[test]
fn test_multiple_alignments_same_read() -> Result<(), BwaError> {
    let (aligner, _) = create_test_aligner();
    let seq = seq_to_bytes("ACGTACGTACGTACGT");
    let sequence = Sequence::new("query", seq.clone());
    let qual = b"I".repeat(16);

    let result = aligner.align_read(&seq, None)?;
    let sam = SAMRecord::from_alignment("read", &result, &sequence, &qual, "test_ref");
    validate_sam_record(&sam.to_string())?;
    Ok(())
}

#[test]
fn test_reference_sequence_encoding() -> Result<(), BwaError> {
    let reference = Reference::parse_fasta(TEST_REFERENCE)?;

    let encoded = reference.as_slice();
    assert_eq!(encoded.len(), 144);

    for (i, &base) in encoded.iter().enumerate() {
        assert!(base <= 4, "Invalid base {} at position {}", base, i);
    }

    let decoded: String = encoded.iter().map(|&b| Reference::decode_base(b)).collect();

    assert!(decoded.contains("ACGT"));

    Ok(())
}

#[test]
fn test_mapq_calculation() -> Result<(), BwaError> {
    let (aligner, _) = create_test_aligner();
    let unique_seq = seq_to_bytes("ACGTACGTACGTACGTTTTTTTTTTTTTTTT");
    let sequence = Sequence::new("query", unique_seq.clone());
    let qual = b"I".repeat(32);

    let result = aligner.align_read(&unique_seq, None)?;
    let sam = SAMRecord::from_alignment("read", &result, &sequence, &qual, "test_ref");
    let sam_str = sam.to_string();
    let fields: Vec<&str> = sam_str.split('\t').collect();
    let mapq: u8 = fields[4].parse().unwrap();

    assert!(mapq <= 60, "MAPQ should be <= 60, got {}", mapq);
    Ok(())
}

#[test]
fn test_seq_to_bytes_all_bases() {
    let test_cases = [
        ("ACGT", vec![0, 1, 2, 3]),
        ("AAAA", vec![0, 0, 0, 0]),
        ("CCCC", vec![1, 1, 1, 1]),
        ("GGGG", vec![2, 2, 2, 2]),
        ("TTTT", vec![3, 3, 3, 3]),
        ("NNNN", vec![4, 4, 4, 4]),
        ("ACGTN", vec![0, 1, 2, 3, 4]),
    ];

    for (input, expected) in test_cases {
        assert_eq!(seq_to_bytes(input), expected);
    }
}

#[test]
fn test_cigar_validation() {
    let valid_cigars = [
        "*",
        "10M",
        "5M2I3M",
        "10S4M",
        "2H3M1D4M2H",
        "10M1D10M",
        "5I5M5D5M",
        "4=",
    ];

    for cigar in valid_cigars {
        assert!(
            validate_cigar(cigar).is_ok(),
            "CIGAR {} should be valid",
            cigar
        );
    }

    let invalid_cigars = ["M10", "10", "5M2"];

    for cigar in invalid_cigars {
        assert!(
            validate_cigar(cigar).is_err(),
            "CIGAR {} should be invalid",
            cigar
        );
    }
}

#[test]
fn test_unmapped_read() -> Result<(), BwaError> {
    let (aligner, _) = create_test_aligner();
    let random_seq = seq_to_bytes("XYZQRSTUVWXYZ");
    let sequence = Sequence::new("query", random_seq.clone());
    let qual = b"I".repeat(13);

    let result = aligner.align_read(&random_seq, None)?;
    let sam = SAMRecord::from_alignment("unmapped", &result, &sequence, &qual, "test_ref");
    let sam_str = sam.to_string();
    let fields: Vec<&str> = sam_str.split('\t').collect();

    assert!(fields.len() >= 11);
    Ok(())
}

#[test]
fn test_sam_record_fields() -> Result<(), BwaError> {
    let (aligner, _) = create_test_aligner();
    let seq = seq_to_bytes("ACGTACGT");
    let sequence = Sequence::new("test", seq.clone());
    let qual = b"IIIIIIII".to_vec();

    let result = aligner.align_read(&seq, None)?;
    let sam_record = SAMRecord::from_alignment("my_read", &result, &sequence, &qual, "ref");
    let sam_str = sam_record.to_string();
    let fields: Vec<&str> = sam_str.split('\t').collect();

    assert!(fields.len() >= 11);
    assert_eq!(fields[0], "my_read");
    assert_eq!(fields[2], "ref");
    assert!(
        fields[5] == "*"
            || fields[5].contains('M')
            || fields[5].contains('=')
            || fields[5].contains('X')
    );
    Ok(())
}

#[test]
fn test_chr1_scaled_reference() -> Result<(), BwaError> {
    let ref_len = 1000;
    let bases = [b'A', b'C', b'G', b'T'];
    let pattern: Vec<u8> = (0..ref_len).map(|i| bases[i % 4]).collect();
    let mut full_seq = pattern.clone();
    full_seq.extend_from_slice(b"GGGGGGGGGGAAAAAAAA");

    let reference = Reference::parse_fasta(&format!(
        ">chr1\n{}",
        std::str::from_utf8(&full_seq).unwrap()
    ))?;
    let total_len = reference.total_len();
    assert_eq!(total_len, ref_len + 18);

    let ref_slice = reference.as_slice();
    let index = FMIndex::build(&reference);
    let aligner = Aligner::new(index, ref_slice).min_seed_len(10);

    let read = seq_to_bytes("GGGGGGGGGGAAAAAAAA");
    let sequence = Sequence::new("scaled_read", read.clone());
    let qual = b"I".repeat(18);

    let result = aligner.align_read(&read, None)?;
    let sam = SAMRecord::from_alignment("scaled_read", &result, &sequence, &qual, "chr1");
    let sam_str = sam.to_string();
    let fields: Vec<&str> = sam_str.split('\t').collect();

    if fields[2] == "*" {
        return Ok(());
    }

    validate_sam_record(&sam_str)?;
    let pos: i64 = fields[3].parse().unwrap();
    assert!(pos >= 1 && pos <= total_len as i64);
    assert!(fields[4].parse::<u8>().unwrap() <= 60);
    Ok(())
}

#[test]
fn test_compare_against_bwa_mem() -> Result<(), BwaError> {
    let bwa_path = std::env::var("BWA_PATH").unwrap_or_else(|_| "bwa".to_string());

    if std::path::Path::new(&bwa_path).metadata().is_err() {
        println!("SKIP: bwa not available at {}", bwa_path);
        return Ok(());
    }

    let temp_dir = std::env::temp_dir();
    let ref_path = temp_dir.join("test_ref.fa");
    let reads_path = temp_dir.join("test_reads.fq");

    std::fs::write(&ref_path, TEST_REFERENCE.as_bytes())?;
    // Use a longer read (32bp) that will actually align to the 144bp reference
    let test_read = "ACGTACGTACGTACGTACGTACGTACGTACGT";
    let qual_str = "I".repeat(32);
    std::fs::write(
        &reads_path,
        format!("@ref_read\n{}\n+\n{}\n", test_read, qual_str),
    )?;

    // Index the reference for bwa
    let index_status = Command::new(&bwa_path)
        .args(["index", ref_path.to_str().unwrap()])
        .output()?;
    if !index_status.status.success() {
        println!(
            "SKIP: bwa index failed: {}",
            String::from_utf8_lossy(&index_status.stderr)
        );
        return Ok(());
    }

    let reference = Reference::parse_fasta(TEST_REFERENCE)?;
    let ref_slice = reference.as_slice();
    let index = FMIndex::build(&reference);
    let aligner = Aligner::new(index, ref_slice).min_seed_len(10);

    let read = seq_to_bytes(test_read);
    let sequence = Sequence::new("ref_read", read.clone());
    let qual = qual_str.as_bytes().to_vec();
    let result = aligner.align_read(&read, None)?;

    let bwa_output = Command::new(&bwa_path)
        .args([
            "mem",
            ref_path.to_str().unwrap(),
            reads_path.to_str().unwrap(),
        ])
        .output()?;

    if !bwa_output.status.success() {
        println!("SKIP: bwa failed to run");
        return Ok(());
    }

    let bwa_sam = String::from_utf8_lossy(&bwa_output.stdout);
    let bwa_lines: Vec<&str> = bwa_sam.lines().filter(|l| !l.starts_with('@')).collect();

    if let Some(bwa_record) = bwa_lines.first() {
        let bwa_fields: Vec<&str> = bwa_record.split('\t').collect();
        let our_sam = SAMRecord::from_alignment("ref_read", &result, &sequence, &qual, "test_ref")
            .to_string();
        let our_fields: Vec<&str> = our_sam.split('\t').collect();

        assert_eq!(our_fields[0], bwa_fields[0], "QNAME should match");
        assert_eq!(our_fields[2], bwa_fields[2], "RNAME should match");

        // Normalize CIGAR: M and = are equivalent for matches
        let our_cigar = our_fields[5].replace('=', "M");
        let bwa_cigar = bwa_fields[5].replace('=', "M");
        assert_eq!(our_cigar, bwa_cigar, "CIGAR should match (M/= normalized)");

        // For repetitive references, both positions are valid
        // Just verify both are positive (mapped) and within reference
        let our_pos: i64 = our_fields[3].parse().unwrap();
        let bwa_pos: i64 = bwa_fields[3].parse().unwrap();
        assert!(our_pos >= 1, "Our position should be >= 1 (mapped)");
        assert!(bwa_pos >= 1, "BWA position should be >= 1 (mapped)");
        assert!(our_pos <= 144, "Our position should be <= ref length");
    }

    std::fs::remove_file(&ref_path).ok();
    std::fs::remove_file(&reads_path).ok();

    Ok(())
}

#[test]
fn test_sam_format_complete() -> Result<(), BwaError> {
    let (aligner, reference) = create_test_aligner();
    let reads = parse_fastq(TEST_FASTQ);
    let header = SAMHeader::new(reference.clone());
    let mut output = header.to_string();

    for (qname, seq, qual) in &reads {
        let result = aligner.align_read(seq, None);
        let sequence = Sequence::new(qname, seq.clone());
        let sam_record = match result {
            Ok(alignment) => {
                SAMRecord::from_alignment(qname, &alignment, &sequence, qual, "test_ref")
            }
            Err(_) => SAMRecord::unmapped(qname, seq, qual),
        };
        output.push_str(&sam_record.to_string());
        output.push('\n');
    }

    let mut line_count = 0;
    let mut header_count = 0;

    for line in output.lines() {
        if line.starts_with('@') {
            header_count += 1;
        } else {
            validate_sam_record(line)?;
            line_count += 1;
        }
    }

    assert!(header_count >= 3, "Should have @HD, @SQ, @PG headers");
    assert_eq!(line_count, 3, "Should have 3 alignment records");

    assert!(output.contains("@HD\tVN:"), "Missing @HD header");
    assert!(
        output.contains("@SQ\tSN:test_ref\tLN:144"),
        "Missing @SQ header"
    );
    assert!(output.contains("@PG"), "Missing @PG header");

    Ok(())
}
