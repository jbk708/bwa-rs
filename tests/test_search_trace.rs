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
fn test_search_trace_gggg() {
    let ref_seq = "AAAAAGGGGAAAAACCCCAAAAA";
    let reference = Reference::parse_fasta(&format!(">test\n{}", ref_seq)).unwrap();
    let index = FMIndex::build(&reference);

    eprintln!("=== Reference: {} ===", ref_seq);
    eprintln!("Reference slice: {:?}", reference.as_slice());
    eprintln!("Length: {}", reference.total_len());
    eprintln!("");

    eprintln!("=== F-column ===");
    eprintln!("F = {:?}", index.f_column);
    eprintln!("");

    eprintln!("=== SA (suffix array) ===");
    eprintln!("SA = {:?}", index.sa.as_slice());
    for (i, &pos) in index.sa.as_slice().iter().enumerate() {
        let suffix = &ref_seq[pos as usize..];
        eprintln!("SA[{:2}] = {:2} -> suffix '{}'", i, pos, suffix);
    }
    eprintln!("");

    eprintln!("=== BWT (burrows-wheeler transform) ===");
    eprintln!("BWT = {:?}", index.bwt.as_slice());
    for (i, &c) in index.bwt.as_slice().iter().enumerate() {
        let sa_pos = index.sa.as_slice()[i];
        let char_at_sa = ref_seq.chars().nth(sa_pos as usize).unwrap_or('?');
        eprintln!("BWT[{:2}] = {} (preceding '{}')", i, c, char_at_sa);
    }
    eprintln!("");

    // Now trace the search for GGGG
    let pattern = "GGGG";
    let pattern_bytes = seq_to_bytes(pattern);

    eprintln!("=== Searching for '{}' ===", pattern);
    eprintln!("Pattern bytes: {:?}", pattern_bytes);
    eprintln!("");

    // Manual search trace
    let mut left = 0usize;
    let mut right = index.len;

    for (step, &c) in pattern_bytes.iter().rev().enumerate() {
        let char_name = match c {
            0 => "A",
            1 => "C",
            2 => "G",
            3 => "T",
            _ => "N",
        };
        let occ_left = index.occ.occ(c, left) as usize;
        let occ_right = index.occ.occ(c, right) as usize;
        let f_c = index.f_column[c as usize] as usize;

        eprintln!("Step {}: Search for '{}' (c={})", step, char_name, c);
        eprintln!("  F[{}] = {}", c, f_c);
        eprintln!("  occ({}, {}) = {}", c, left, occ_left);
        eprintln!("  occ({}, {}) = {}", c, right, occ_right);
        eprintln!("  left = {} + {} = {}", f_c, occ_left, f_c + occ_left);
        eprintln!("  right = {} + {} = {}", f_c, occ_right, f_c + occ_right);

        left = f_c + occ_left;
        right = f_c + occ_right;
        eprintln!("  New range: [{}, {})", left, right);
        eprintln!("");
    }

    // Compare with built-in search
    let (result_left, result_right) = index.search(&pattern_bytes);
    eprintln!(
        "Built-in search result: [{}, {})",
        result_left, result_right
    );

    // Get positions
    let positions = index.find_all(&pattern_bytes);
    eprintln!("Found positions: {:?}", positions);

    // Verify by checking actual locations in reference
    eprintln!("\n=== Verification ===");
    for (i, _c) in ref_seq.char_indices() {
        let chars: String = ref_seq.chars().skip(i).take(4).collect();
        if chars == pattern {
            eprintln!("Pattern '{}' found at position {} in reference", pattern, i);
        }
    }
}

#[test]
fn test_full_bwt_analysis() {
    let ref_seq = "AAAAAGGGGAAAAACCCCAAAAA";
    let reference = Reference::parse_fasta(&format!(">test\n{}", ref_seq)).unwrap();
    let index = FMIndex::build(&reference);

    // Count actual G's in BWT
    let bwt = index.bwt.as_slice();
    let mut g_count = 0;
    for (i, &c) in bwt.iter().enumerate() {
        if c == 2 {
            g_count += 1;
            let sa_pos = index.sa.as_slice()[i];
            eprintln!(
                "BWT[{}] = G, SA[{}] = {}, preceding char = {:?}",
                i,
                i,
                sa_pos,
                ref_seq.chars().nth(sa_pos as usize)
            );
        }
    }
    eprintln!("\nTotal G's in BWT: {}", g_count);
    eprintln!("occ(2, 23) from function: {}", index.occ.occ(2, 23));
    eprintln!("occ(2, 22) from function: {}", index.occ.occ(2, 22));

    // Verify the expected SA indices for G suffixes
    eprintln!("\nExpected G suffix SA indices:");
    for (i, &pos) in index.sa.as_slice().iter().enumerate() {
        if ref_seq.chars().nth(pos as usize) == Some('G') {
            eprintln!("  SA[{}] = {}", i, pos);
        }
    }

    // Check what the range should be after searching for single G
    eprintln!("\nRange for 'G' pattern:");
    eprintln!("  F[2] = {}", index.f_column[2]);
    eprintln!("  Expected: [F[2], F[2] + count_of_G_in_BWT)");
    eprintln!(
        "  = [{}, {} + {})",
        index.f_column[2], index.f_column[2], g_count
    );
}

#[test]
fn test_wavelet_rank_trace() {
    use wavelet_matrix::WaveletMatrix;

    let bwt: Vec<u64> = vec![
        0, 0, 0, 0, 1, 2, 4, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 2, 2, 2, 0,
    ];
    let wm = WaveletMatrix::new(&bwt);

    eprintln!("BWT length: {}", bwt.len());
    eprintln!("Testing wm.rank() for value 2 (G):");
    for k in [0, 5, 6, 19, 20, 22, 23] {
        let rank = wm.rank(k, 2);
        eprintln!("  wm.rank({}, 2) = {}", k, rank);
    }
}
