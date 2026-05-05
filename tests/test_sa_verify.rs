use bwa_mem::{FMIndex, Reference};

#[test]
fn test_verify_sa_and_fcolumn() {
    let ref_seq = "AAAAAGGGGAAAAACCCCAAAAA";
    let reference = Reference::parse_fasta(&format!(">test\n{}", ref_seq)).unwrap();
    let index = FMIndex::build(&reference);
    
    println!("Reference: {}", ref_seq);
    println!("Length: {}", ref_seq.len());
    println!("Index len (n+1): {}", index.len);
    println!("F-column: {:?}", index.f_column);
    println!("\nSA (index -> position -> suffix):");
    
    let mut a_positions = Vec::new();
    let mut c_positions = Vec::new();
    let mut g_positions = Vec::new();
    
    for i in 0..index.sa.len() {
        if let Some(pos) = index.sa.get(i) {
            if pos as usize >= ref_seq.len() {
                // Sentinel position
                println!("  SA[{:2}] = {:2} -> 'SENTINEL'", i, pos);
                continue;
            }
            let suffix = &ref_seq[pos as usize..];
            let first_char = ref_seq.chars().nth(pos as usize).unwrap();
            println!("  SA[{:2}] = {:2} -> '{}' (starts with {})", 
                i, pos, suffix.chars().take(10).collect::<String>(), first_char);
            
            match first_char {
                'A' => a_positions.push((i, pos)),
                'C' => c_positions.push((i, pos)),
                'G' => g_positions.push((i, pos)),
                _ => {}
            }
        }
    }
    
    println!("\nA-suffixes: {:?}", a_positions);
    println!("C-suffixes: {:?}", c_positions);
    println!("G-suffixes: {:?}", g_positions);
    
    // Verify F-column
    println!("\nF-column verification:");
    if !a_positions.is_empty() {
        println!("  First A at SA index {}, so F[0] should be {} (1-indexed)", 
            a_positions[0].0, a_positions[0].0 + 1);
    }
    if !c_positions.is_empty() {
        println!("  First C at SA index {}, so F[1] should be {} (1-indexed)", 
            c_positions[0].0, c_positions[0].0 + 1);
    }
    if !g_positions.is_empty() {
        println!("  First G at SA index {}, so F[2] should be {} (1-indexed)", 
            g_positions[0].0, g_positions[0].0 + 1);
    }
}
