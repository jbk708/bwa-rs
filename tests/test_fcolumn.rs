use bwa_mem::{FMIndex, Reference};

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

#[test]
fn test_simple_pattern_search() {
    // Simple reference with known pattern
    let ref_seq = "AAAAAGGGGAAAAACCCCAAAAA";
    let reference = Reference::parse_fasta(&format!(">test\n{}", ref_seq)).unwrap();
    let index = FMIndex::build(&reference);

    // Test searching for various patterns
    let patterns = ["GGGG", "AAAAA", "CCCC", "ACCCC", "GGGGAAAAA", "AAAAACCCC"];

    for pattern in patterns {
        let bytes = seq_to_bytes(pattern);
        let (left, right) = index.search(&bytes);
        let positions = index.find_all(&bytes);

        println!(
            "Pattern '{}': range=[{}, {}), positions={:?}",
            pattern, left, right, positions
        );

        // Verify by checking actual locations in reference
        let mut expected: Vec<u32> = Vec::new();
        for (i, _) in ref_seq.match_indices(pattern) {
            expected.push(i as u32);
        }

        // Sort both since FM-index returns positions in SA order, not position order
        let mut sorted_positions = positions.clone();
        sorted_positions.sort();
        expected.sort();

        assert_eq!(
            sorted_positions, expected,
            "Pattern '{}' positions mismatch",
            pattern
        );
    }
}

#[test]
fn test_large_pattern_search() {
    // Create a reference with pattern at known position
    let ref_len = 50_000;
    let pattern_pos = 10000u32;
    let pattern = "GGGGAAAAACCCC";

    let mut ref_bytes = Vec::with_capacity(ref_len);
    for i in 0..ref_len {
        ref_bytes.push([b'A', b'C', b'G', b'T'][(i * 17 + 13) % 4]);
    }
    ref_bytes.splice(
        pattern_pos as usize..pattern_pos as usize + pattern.len(),
        pattern.bytes(),
    );

    let fasta = format!(">chr1\n{}", std::str::from_utf8(&ref_bytes).unwrap());
    let reference = Reference::parse_fasta(&fasta).unwrap();
    let index = FMIndex::build(&reference);

    eprintln!("F-column: {:?}", index.f_column);

    let pattern_bytes = seq_to_bytes(pattern);
    let (left, right) = index.search(&pattern_bytes);
    let positions = index.find_all(&pattern_bytes);

    eprintln!(
        "Pattern '{}' at pos {}: range=[{}, {}), positions={:?}",
        pattern, pattern_pos, left, right, positions
    );

    assert!(!positions.is_empty(), "Should find pattern");
    assert!(
        positions.contains(&pattern_pos),
        "Should find at position {}",
        pattern_pos
    );
}

#[test]
fn test_debug_fm_index() {
    let ref_seq = "AAAAAGGGGAAAAACCCCAAAAA";
    let reference = Reference::parse_fasta(&format!(">test\n{}", ref_seq)).unwrap();
    let index = FMIndex::build(&reference);

    eprintln!("Reference: {}", ref_seq);
    eprintln!("len (FM-index): {}", index.len);
    eprintln!("F-column: {:?}", index.f_column);
    eprintln!("BWT len: {}", index.bwt.len());
    eprintln!("SA len: {}", index.sa.len());
    eprintln!("BWT: {:?}", index.bwt.as_slice());
    eprintln!("SA: {:?}", index.sa.as_slice());

    // Manual trace G search
    let c = 2u8;
    let occ_left = index.occ.occ(c, 0);
    let occ_right = index.occ.occ(c, index.len);
    let f_c = (index.f_column[c as usize] as usize).saturating_sub(1);
    eprintln!("\nManual trace for G (c=2):");
    eprintln!("  F[2] - 1 = {}", f_c);
    eprintln!("  occ(2, 0) = {}", occ_left);
    eprintln!("  occ(2, {}) = {}", index.len, occ_right);
    eprintln!(
        "  new range: [{} + {}, {} + {}) = [{}, {})",
        f_c,
        occ_left,
        f_c,
        occ_right,
        f_c + occ_left as usize,
        f_c + occ_right as usize
    );

    // Show suffixes in SA that start with G
    eprintln!("\nSuffixes starting with G:");
    for (i, &pos) in index.sa.as_slice().iter().enumerate() {
        if pos < ref_seq.len() as u32 {
            let suffix = &ref_seq[pos as usize..];
            if suffix.starts_with("G") {
                eprintln!(
                    "  SA[{}] = {} -> '{}'",
                    i,
                    pos,
                    &suffix[..10.min(suffix.len())]
                );
            }
        } else {
            eprintln!("  SA[{}] = {} -> 'SENTINEL'", i, pos);
        }
    }
}

#[test]
fn test_debug_gggg_search() {
    let ref_seq = "AAAAAGGGGAAAAACCCCAAAAA";
    let reference = Reference::parse_fasta(&format!(">test\n{}", ref_seq)).unwrap();
    let index = FMIndex::build(&reference);

    // Search for GGGG
    let pattern = seq_to_bytes("GGGG");
    let (left, right) = index.search(&pattern);
    let positions = index.find_all(&pattern);

    eprintln!("\nSearch for GGGG:");
    eprintln!("Pattern: {:?}", pattern);
    eprintln!("Range: [{}, {})", left, right);
    eprintln!("Positions: {:?}", positions);

    // Verify each position
    for &pos in &positions {
        let end = (pos as usize + 4).min(ref_seq.len());
        let found = &ref_seq[pos as usize..end];
        eprintln!("  Position {}: substring = '{}'", pos, found);
    }

    // Expected position is 5 (GGGGAAAAACCCC starts at position 5)
    eprintln!("\nExpected position: 5");
    eprintln!("Expected substring: '{}'", &ref_seq[5..9]);
}

#[test]
fn test_debug_wavelet_rank() {
    let ref_seq = "AAAAAGGGGAAAAACCCCAAAAA";
    let reference = Reference::parse_fasta(&format!(">test\n{}", ref_seq)).unwrap();
    let index = FMIndex::build(&reference);

    let bwt = index.bwt.as_slice();
    eprintln!("BWT len: {}", bwt.len());
    eprintln!("BWT: {:?}", bwt);

    // Count G manually
    let mut manual_count = 0u32;
    for i in 0..bwt.len() {
        if bwt[i] == 2 {
            manual_count += 1;
        }
    }
    eprintln!("Manual G count in BWT: {}", manual_count);

    // Check wavelet tree rank at different positions
    for i in 0..=bwt.len() {
        let wt_rank = index.occ.occ(2, i);
        let mut expected = 0u32;
        for j in 0..i {
            if bwt[j] == 2 {
                expected += 1;
            }
        }
        if wt_rank != expected {
            eprintln!(
                "MISMATCH at i={}: wt_rank={}, expected={}",
                i, wt_rank, expected
            );
        }
    }
    eprintln!("Done checking G ranks");

    // Check A rank too
    eprintln!("\nChecking A ranks:");
    for i in 0..=bwt.len() {
        let wt_rank = index.occ.occ(0, i);
        let mut expected = 0u32;
        for j in 0..i {
            if bwt[j] == 0 {
                expected += 1;
            }
        }
        if wt_rank != expected {
            eprintln!(
                "MISMATCH at i={}: wt_rank={}, expected={}",
                i, wt_rank, expected
            );
        }
    }
    eprintln!("Done checking A ranks");
}

#[test]
fn test_debug_search_trace() {
    let ref_seq = "AAAAAGGGGAAAAACCCCAAAAA";
    let reference = Reference::parse_fasta(&format!(">test\n{}", ref_seq)).unwrap();
    let index = FMIndex::build(&reference);

    let pattern = seq_to_bytes("GGGG");
    let _n = pattern.len();

    eprintln!("Pattern: {:?}", pattern);
    eprintln!("F-column: {:?}", index.f_column);
    eprintln!("BWT: {:?}", index.bwt.as_slice());
    eprintln!("SA: {:?}", index.sa.as_slice());

    let mut left = 0usize;
    let mut right = index.len;

    eprintln!("\nInitial range: [{}, {})", left, right);

    for (step, &c) in pattern.iter().rev().enumerate() {
        let char_name = match c {
            0 => "A",
            1 => "C",
            2 => "G",
            3 => "T",
            4 => "$",
            _ => "?",
        };

        let occ_left = index.occ.occ(c, left);
        let occ_right = index.occ.occ(c, right);
        let f_c = (index.f_column[c as usize] as usize).saturating_sub(1);

        eprintln!("\nStep {}: c={} ({})", step, c, char_name);
        eprintln!("  Before: left={}, right={}", left, right);
        eprintln!(
            "  F[{}]={}, F[{}]-1={}",
            c, index.f_column[c as usize], c, f_c
        );
        eprintln!("  occ({}, {}) = {}", c, left, occ_left);
        eprintln!("  occ({}, {}) = {}", c, right, occ_right);

        left = f_c + occ_left as usize;
        right = f_c + occ_right as usize;

        eprintln!("  After: left={}, right={}", left, right);

        // Show suffixes in range
        if left < right && right <= index.sa.len() {
            eprintln!("  Suffixes in range [{}, {}):", left, right);
            for i in left..right.min(left + 5) {
                if let Some(pos) = index.sa.get(i) {
                    let suffix = if pos < ref_seq.len() as u32 {
                        &ref_seq[pos as usize..(pos as usize + 10).min(ref_seq.len())]
                    } else {
                        "SENTINEL"
                    };
                    eprintln!("    SA[{}] = {} -> '{}'", i, pos, suffix);
                }
            }
        }

        if left >= right {
            eprintln!("  RANGE COLLAPSED!");
            break;
        }
    }

    let positions = index.find_all(&pattern);
    eprintln!("\nFinal range: [{}, {})", left, right);
    eprintln!("Final positions: {:?}", positions);
}

#[test]
fn test_debug_padded_sa() {
    let ref_seq = "AAAAAGGGGAAAAACCCCAAAAA";
    let reference = Reference::parse_fasta(&format!(">test\n{}", ref_seq)).unwrap();
    let sequence = reference.as_slice();
    let len = sequence.len();

    // Build padded sequence (what the FM-index does)
    let mut padded = Vec::with_capacity(len + 1);
    padded.extend_from_slice(&sequence);
    padded.push(4); // Sentinel

    eprintln!("Original sequence: {:?}", sequence);
    eprintln!("Padded sequence: {:?}", padded);
    eprintln!("Padded length: {}", padded.len());

    // Build SA of padded
    let sa_padded = bwa_mem::sa::SuffixArray::build(&padded);
    eprintln!("\nSA of padded sequence:");
    for (i, &pos) in sa_padded.as_slice().iter().enumerate() {
        let suffix = if pos < padded.len() as u32 {
            let end = (pos as usize + 15).min(padded.len());
            format!("{:?}", &padded[pos as usize..end])
        } else {
            "OUT OF BOUNDS".to_string()
        };
        eprintln!("  SA[{}] = {} -> {}", i, pos, suffix);
    }

    // Now check the suffix at position 5 of original
    eprintln!("\nSuffix at position 5 of original:");
    eprintln!("  Original[5:] = {:?}", &sequence[5..]);
    eprintln!("  Padded[5:6] = {:?}", &padded[5..6]);

    // Find where position 5 (of original) ended up in the SA
    // It should be at a position where SA value is 5
    for (i, &pos) in sa_padded.as_slice().iter().enumerate() {
        if pos == 5 {
            eprintln!("  Position 5 found at SA[{}]", i);
            // Check what suffix that is
            let suffix = &padded[5..];
            eprintln!("  Suffix at position 5: {:?}", suffix);
        }
    }
}

#[test]
fn test_debug_acac_occ() {
    let ref_seq = Reference::parse_fasta(">test\nACAC").unwrap();
    let index = FMIndex::build(&ref_seq);

    eprintln!("Index len: {}", index.len);
    eprintln!("BWT: {:?}", index.bwt.as_slice());
    eprintln!("SA: {:?}", index.sa.as_slice());

    let bwt = index.bwt.as_slice();

    // Manual occ for C
    for i in 0..=5 {
        let mut count = 0u32;
        for j in 0..i {
            if bwt[j] == 1 {
                count += 1;
            }
        }
        let wt_occ = index.occ.occ(1, i);
        eprintln!("occ(C, {}) - manual: {}, wavelet: {}", i, count, wt_occ);
    }

    // Manual occ for A
    for i in 0..=5 {
        let mut count = 0u32;
        for j in 0..i {
            if bwt[j] == 0 {
                count += 1;
            }
        }
        let wt_occ = index.occ.occ(0, i);
        eprintln!("occ(A, {}) - manual: {}, wavelet: {}", i, count, wt_occ);
    }
}

#[test]
fn test_debug_sa_acac() {
    // Test SA construction for ACAC
    let seq = vec![0u8, 1, 0, 1, 4]; // ACAC$
    let sa = bwa_mem::sa::SuffixArray::build(&seq);

    eprintln!("Sequence: {:?}", seq);
    eprintln!("SA: {:?}", sa.as_slice());

    // Verify suffixes
    for (i, &pos) in sa.as_slice().iter().enumerate() {
        let end = (pos as usize + 10).min(seq.len());
        let suffix = &seq[pos as usize..end];
        eprintln!("SA[{}] = {} -> suffix: {:?}", i, pos, suffix);
    }
}

#[test]
fn test_debug_sa_ordering() {
    // Manually verify suffix ordering
    let seq = vec![0u8, 1, 0, 1, 4]; // ACAC$

    // All suffixes
    let suffixes: Vec<(&[u8], usize)> = (0..seq.len()).map(|i| (&seq[i..], i)).collect();

    // Sort manually
    let mut sorted = suffixes.clone();
    sorted.sort_by(|a, b| a.0.cmp(b.0));

    println!("Sorted suffixes:");
    for (i, (suffix, pos)) in sorted.iter().enumerate() {
        println!("  SA[{}] = {} -> {:?}", i, pos, suffix);
    }

    // Check if SA matches
    let sa = bwa_mem::sa::SuffixArray::build(&seq);
    println!("\nlibsais SA: {:?}", sa.as_slice());

    let expected: Vec<usize> = sorted.iter().map(|(_, pos)| *pos).collect();
    let actual: Vec<u32> = sa.as_slice().to_vec();

    if expected != actual.iter().map(|&x| x as usize).collect::<Vec<_>>() {
        println!("\nSA MISMATCH!");
        println!("Expected: {:?}", expected);
        println!("Actual: {:?}", actual);
    } else {
        println!("\nSA matches expected!");
    }
}

#[test]
fn test_debug_f_column_with_4() {
    // With sentinel = 4 being larger than A, C, G, T:
    // The F-column should reflect: A=1, C=?, G=?, T=?, $=?

    let ref_seq = "ACAC";
    let sequence: Vec<u8> = ref_seq
        .bytes()
        .map(|b| match b {
            b'A' => 0,
            b'C' => 1,
            b'G' => 2,
            b'T' => 3,
            _ => 4,
        })
        .collect();

    println!("Sequence: {:?}", sequence);
    println!(
        "Counts: A={}, C={}, G={}, T={}",
        sequence.iter().filter(|&&x| x == 0).count(),
        sequence.iter().filter(|&&x| x == 1).count(),
        sequence.iter().filter(|&&x| x == 2).count(),
        sequence.iter().filter(|&&x| x == 3).count()
    );

    // F[c] = number of characters < c + 1
    // A (0): 0 chars < 0, F[0] = 1
    // C (1): 1 char < 1 (A=1), F[1] = 2
    // G (2): 2 chars < 2 (A+C=2), F[2] = 3
    // T (3): 3 chars < 3 (A+C+G=3), F[3] = 4
    // $ (4): 4 chars < 4 (A+C+G+T=4), F[4] = 5 = n+1

    println!("Expected F-column: [1, 2, 3, 4, 5]");

    // But in our code we have:
    // F[0] = 1
    // F[1] = 1 + total[0] = 1 + 2 = 3
    // F[2] = 1 + total[0] + total[1] = 1 + 2 + 2 = 5
    // etc.

    // The issue is the formula: F[i] = 1 + sum(total[0..i])
    // This gives F[0] = 1, F[1] = 1 + A, F[2] = 1 + A + C, etc.
    // But if sentinel (4) is at the end, F[4] = 1 + A + C + G + T = n + 1 = 5

    // Actually this seems correct for sentinel being last!
    // The issue must be elsewhere.
}

#[test]
fn test_debug_fm_for_acac() {
    let ref_seq = Reference::parse_fasta(">test\nACAC").unwrap();
    let index = FMIndex::build(&ref_seq);

    println!("Index len: {}", index.len);
    println!("F-column: {:?}", index.f_column);
    println!("BWT: {:?}", index.bwt.as_slice());
    println!("SA: {:?}", index.sa.as_slice());

    // Expected F for ACAC:
    // A count = 2, C count = 2
    // F[0] (A) = 1
    // F[1] (C) = A + 1 = 3
    // F[2] (G) = A + C + 1 = 5
    // F[3] (T) = A + C + G + 1 = 5
    // F[4] ($) = n + 1 = 5
    println!("\nExpected F: [1, 3, 5, 5, 5]");

    // Now trace search for C
    // Search for "C" (c=1)
    // left = 0, right = 5
    // occ(1, 0) = 0, occ(1, 5) = 2
    // F[1] - 1 = 3 - 1 = 2
    // new range = [2, 2 + 2) = [2, 4)

    let (left, right) = index.search(&[1]); // Search for C
    println!("\nSearch for C: range = [{}, {})", left, right);

    // Check suffixes in range
    for i in left..right {
        if let Some(pos) = index.sa.get(i) {
            println!("  SA[{}] = {}", i, pos);
        }
    }
}

#[test]
fn test_debug_bwt_acac() {
    let ref_seq = Reference::parse_fasta(">test\nACAC").unwrap();
    let index = FMIndex::build(&ref_seq);

    let bwt = index.bwt.as_slice();

    println!("BWT len: {}", bwt.len());
    println!("BWT: {:?}", bwt);

    // Manual BWT verification
    // SA = [0, 2, 1, 3, 4]
    // BWT[i] = character before SA[i]
    // SA[0] = 0: before pos 0 = sentinel = 4
    // SA[1] = 2: before pos 2 = seq[1] = 1 (C)
    // SA[2] = 1: before pos 1 = seq[0] = 0 (A)
    // SA[3] = 3: before pos 3 = seq[2] = 0 (A)
    // SA[4] = 4: before pos 4 = seq[3] = 1 (C)

    println!("\nExpected BWT: [4, 1, 0, 0, 1]");
    println!("Actual BWT:   {:?}", bwt);

    // Check each position
    for (i, &c) in bwt.iter().enumerate() {
        let char_name = match c {
            0 => "A",
            1 => "C",
            2 => "G",
            3 => "T",
            4 => "$",
            _ => "?",
        };
        let suffix_pos = index.sa.get(i).unwrap_or(999);
        println!(
            "BWT[{}] = {} ({}) for SA[{}] = {}",
            i, c, char_name, i, suffix_pos
        );
    }

    // Manual count of C in BWT
    let mut c_count = 0u32;
    for i in 0..bwt.len() {
        if bwt[i] == 1 {
            c_count += 1;
        }
    }
    println!("\nManual C count in BWT: {}", c_count);

    // Wavelet tree occ
    for i in 0..=bwt.len() {
        println!("occ(C, {}) = {}", i, index.occ.occ(1, i));
    }
}

#[test]
fn test_debug_fm_search_acgt() {
    let reference = Reference::parse_fasta(">test\nACGTACGTACGT").unwrap();
    let fm_index = FMIndex::build(&reference);

    let pattern = [0, 1, 2, 3]; // ACGT

    eprintln!("FMIndex Pattern: {:?}", pattern);
    eprintln!("FMIndex F-column: {:?}", fm_index.f_column);
    eprintln!("FMIndex len: {}", fm_index.len);
    eprintln!("FMIndex BWT: {:?}", fm_index.bwt.as_slice());

    // Search backwards
    let mut left = 0usize;
    let mut right = fm_index.len;

    eprintln!("\nInitial range: [{}, {})", left, right);

    for (step, &c) in pattern.iter().rev().enumerate() {
        let occ_left = fm_index.occ.occ(c, left);
        let occ_right = fm_index.occ.occ(c, right);
        // Note: FMIndex uses F[c] - 1
        let f_c = (fm_index.f_column[c as usize] as usize).saturating_sub(1);

        eprintln!("\nStep {}: c={}", step, c);
        eprintln!(
            "  F[{}] = {}, F[{}]-1 = {}",
            c, fm_index.f_column[c as usize], c, f_c
        );
        eprintln!("  occ({}, {}) = {}", c, left, occ_left);
        eprintln!("  occ({}, {}) = {}", c, right, occ_right);

        left = f_c + occ_left as usize;
        right = f_c + occ_right as usize;

        eprintln!("  New range: [{}, {})", left, right);

        if left >= right {
            eprintln!("  RANGE COLLAPSED!");
            break;
        }
    }

    let positions = fm_index.find_all(&pattern);
    eprintln!("\nFinal positions: {:?}", positions);
}
