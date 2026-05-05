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
fn test_gggg_search() {
    let ref_seq = "AAAAAGGGGAAAAACCCCAAAAA";
    let reference = Reference::parse_fasta(&format!(">test\n{}", ref_seq)).unwrap();
    let index = FMIndex::build(&reference);
    
    println!("Reference: {}", ref_seq);
    println!("F-column: {:?}", index.f_column);
    
    // Find all occurrences of "GGGG" in reference
    let mut expected = Vec::new();
    for (i, _) in ref_seq.match_indices("GGGG") {
        expected.push(i as u32);
    }
    println!("Expected positions for GGGG: {:?}", expected);
    
    // Search for GGGG
    let pattern = seq_to_bytes("GGGG");
    let (left, right) = index.search(&pattern);
    let positions = index.find_all(&pattern);
    
    println!("Search for GGGG: range=[{}, {}), positions={:?}", left, right, positions);
    
    // Verify positions are in expected
    for pos in &positions {
        let suffix = &ref_seq[*pos as usize..];
        println!("  Position {}: suffix starts with '{}'", pos, &suffix[..10.min(suffix.len())]);
        assert!(suffix.starts_with("GGGG"), "Position {} should start with GGGG", pos);
    }
    
    assert_eq!(positions, expected, "GGGG positions mismatch");
}

#[test]
fn test_ggggaaaaacccc_search() {
    let ref_seq = "AAAAAGGGGAAAAACCCCAAAAA";
    let reference = Reference::parse_fasta(&format!(">test\n{}", ref_seq)).unwrap();
    let index = FMIndex::build(&reference);
    
    // Find all occurrences of "GGGGAAAAACCCC"
    let mut expected = Vec::new();
    for (i, _) in ref_seq.match_indices("GGGGAAAAACCCC") {
        expected.push(i as u32);
    }
    println!("\nExpected positions for GGGGAAAAACCCC: {:?}", expected);
    
    // Search for GGGGAAAAACCCC
    let pattern = seq_to_bytes("GGGGAAAAACCCC");
    let (left, right) = index.search(&pattern);
    let positions = index.find_all(&pattern);
    
    println!("Search for GGGGAAAAACCCC: range=[{}, {}), positions={:?}", left, right, positions);
    
    // Verify
    for pos in &positions {
        let suffix = &ref_seq[*pos as usize..];
        assert!(suffix.starts_with("GGGGAAAAACCCC"), "Position {} should start with GGGGAAAAACCCC", pos);
    }
    
    assert_eq!(positions, expected, "GGGGAAAAACCCC positions mismatch");
}

#[test]
fn test_search_trace_long_pattern() {
    let ref_seq = "AAAAAGGGGAAAAACCCCAAAAA";
    let reference = Reference::parse_fasta(&format!(">test\n{}", ref_seq)).unwrap();
    let index = FMIndex::build(&reference);
    
    let pattern = seq_to_bytes("GGGGAAAAACCCC");
    println!("Pattern: {:?}", pattern);
    println!("F-column: {:?}", index.f_column);
    println!("BWT: {:?}", index.bwt.as_slice());
    println!("SA: {:?}", index.sa.as_slice());
    
    // Search backwards through pattern
    let mut left = 0usize;
    let mut right = index.len;
    
    println!("\nStarting search with range [{}, {})", left, right);
    
    for (step, &c) in pattern.iter().rev().enumerate() {
        let char_name = match c { 0 => "A", 1 => "C", 2 => "G", 3 => "T", _ => "N" };
        let occ_left = index.occ.occ(c, left);
        let occ_right = index.occ.occ(c, right);
        let f_c = (index.f_column[c as usize] as usize).saturating_sub(1);
        
        println!("\nStep {}: c={} ({})", step, c, char_name);
        println!("  F[{}] - 1 = {}", c, f_c);
        println!("  occ({}, {}) = {}", c, left, occ_left);
        println!("  occ({}, {}) = {}", c, right, occ_right);
        println!("  new range: [{} + {}, {} + {}) = [{}, {})", 
            f_c, occ_left, f_c, occ_right, f_c + occ_left as usize, f_c + occ_right as usize);
        
        left = f_c + occ_left as usize;
        right = f_c + occ_right as usize;
        
        println!("  -> Range [{}, {})", left, right);
        
        if left >= right {
            println!("  RANGE COLLAPSED!");
            break;
        }
        
        // Show what suffixes are in range
        println!("  Suffixes in range:");
        for i in left..right.min(index.sa.len()) {
            if let Some(pos) = index.sa.get(i) {
                let suffix = &ref_seq[pos as usize..];
                println!("    SA[{}] = {} -> '{}'", i, pos, &suffix[..10.min(suffix.len())]);
            }
        }
    }
    
    let positions = index.find_all(&pattern);
    println!("\nFinal positions: {:?}", positions);
}
