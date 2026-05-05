//! Test for T50: FM-index search not finding patterns in long references
//!
//! This test reproduces the issue where FM-index search returns incorrect results
//! for longer patterns in large references.

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
fn test_large_reference_search() {
    // T50: Reproduce FM-index search issue on large references
    // Create a ~200KB reference with a known pattern at position 100000
    let ref_len = 200_013;

    // Build reference: random-ish sequence with GGGGAAAAACCCC at position 100000
    let mut ref_bytes = Vec::with_capacity(ref_len);
    for i in 0..ref_len {
        ref_bytes.push([b'A', b'C', b'G', b'T'][(i * 17 + 13) % 4]);
    }
    // Insert known pattern at position 100000
    ref_bytes.splice(100000..100013, b"GGGGAAAAACCCC".iter().copied());

    let fasta = format!(">chr1\n{}", std::str::from_utf8(&ref_bytes).unwrap());
    let reference = Reference::parse_fasta(&fasta).unwrap();

    assert_eq!(reference.total_len(), ref_len);

    // Build FM-index
    let index = FMIndex::build(&reference);
    assert_eq!(index.len, ref_len + 1);  // n+1 with sentinel

    // Verify pattern exists in reference at position 100000
    let ref_slice = reference.as_slice();
    assert_eq!(
        &ref_slice[100000..100013],
        seq_to_bytes("GGGGAAAAACCCC").as_slice(),
        "Pattern should be at position 100000"
    );

    // Search for pattern - should find at least one match
    let pattern_bytes = seq_to_bytes("GGGGAAAAACCCC");
    let (left, right) = index.search(&pattern_bytes);

    println!("Search for GGGGAAAAACCCC: range = ({}, {})", left, right);

    // Get actual positions
    let positions = index.find_all(&pattern_bytes);
    println!("Found positions: {:?}", positions);

    // The search should return a valid range
    assert!(left < right, "Search should find matches: ({}, {})", left, right);

    // Should find positions
    assert!(
        !positions.is_empty(),
        "FM-index should find GGGGAAAAACCCC in large reference"
    );

    // Verify at least one position is correct (100000)
    assert!(
        positions.contains(&100000),
        "Should find pattern at position 100000, got: {:?}",
        positions
    );
}

#[test]
fn test_simple_vs_large_consistency() {
    // Verify that simple and large references behave consistently
    let simple_ref = "AAAAAGGGGAAAAACCCCAAAAA";
    let simple_pattern = "GGGGAAAAACCCC";

    let simple_fasta = format!(">test\n{}", simple_ref);
    let simple_ref_obj = Reference::parse_fasta(&simple_fasta).unwrap();
    let simple_index = FMIndex::build(&simple_ref_obj);
    let simple_positions = simple_index.find_all(&seq_to_bytes(simple_pattern));

    println!("Simple ref positions: {:?}", simple_positions);
    assert!(!simple_positions.is_empty(), "Simple reference should find pattern");

    // Now test large reference
    let mut large_ref_bytes = Vec::new();
    for i in 0..50_000 {
        large_ref_bytes.push([b'A', b'C', b'G', b'T'][(i * 17 + 13) % 4]);
    }
    large_ref_bytes.extend_from_slice(simple_ref.as_bytes());

    let large_fasta = format!(">chr1\n{}", std::str::from_utf8(&large_ref_bytes).unwrap());
    let large_ref_obj = Reference::parse_fasta(&large_fasta).unwrap();
    let large_index = FMIndex::build(&large_ref_obj);
    let large_positions = large_index.find_all(&seq_to_bytes(simple_pattern));

    println!("Large ref positions: {:?}", large_positions);
    assert!(!large_positions.is_empty(), "Large reference should also find pattern");
}

#[test]
fn test_ggg_pattern_large_ref() {
    // Test the simpler GGGG pattern from the issue
    let pattern = b"GGGG";
    let ref_len = 200_013;

    // Build reference with known GGGG pattern at position 100000
    let mut ref_bytes = Vec::with_capacity(ref_len);
    for i in 0..ref_len {
        ref_bytes.push([b'A', b'C', b'G', b'T'][(i * 7) % 4]); // More Gs
    }
    ref_bytes.splice(100000..100004, pattern.iter().copied());

    let fasta = format!(">chr1\n{}", std::str::from_utf8(&ref_bytes).unwrap());
    let reference = Reference::parse_fasta(&fasta).unwrap();
    let index = FMIndex::build(&reference);

    // Search for GGGG
    let (left, right) = index.search(&seq_to_bytes("GGGG"));
    let positions = index.find_all(&seq_to_bytes("GGGG"));

    println!("GGGG search: range = ({}, {}), positions = {:?}", left, right, positions);

    assert!(!positions.is_empty(), "Should find GGGG in large reference");

    assert!(
        positions.contains(&100000),
        "Should find GGGG at position 100000"
    );
}

#[test]
fn test_f_column_correctness() {
    // Test that F-column calculation is correct for various reference sizes
    for size in [100, 1000, 10000, 50000, 200000] {
        let ref_bytes: Vec<u8> = (0..size)
            .map(|i| [b'A', b'C', b'G', b'T'][(i * 17 + 13) % 4])
            .collect();

        let fasta = format!(">test\n{}", std::str::from_utf8(&ref_bytes).unwrap());
        let reference = Reference::parse_fasta(&fasta).unwrap();
        let index = FMIndex::build(&reference);

        // F column should be valid
        for c in 0..4u8 {
            if c < 4 {
                assert!(
                    index.f_column[(c + 1) as usize] >= index.f_column[c as usize],
                    "F-column should be non-decreasing"
                );
            }
        }
        // F[4] - 1 should be >= reference length
        assert!(
            index.f_column[4] >= index.len as u32,
            "F[4]={} should be >= len={} for size {}",
            index.f_column[4],
            index.len,
            size
        );

        // Test search for single character
        let (left, right) = index.search(&[2]); // G
        assert!(
            left <= right && right <= index.len,
            "Search range ({}, {}) should be within [0, {}] for size {}",
            left,
            right,
            index.len,
            size
        );
    }
}

#[test]
fn test_large_ref_debug() {
    // Test with a simpler large reference
    let ref_len = 50_000;

    // Build reference with GGGGAAAAACCCC at position 10000
    let mut ref_bytes = Vec::with_capacity(ref_len);
    for i in 0..ref_len {
        ref_bytes.push([b'A', b'C', b'G', b'T'][(i * 17 + 13) % 4]);
    }
    // Insert pattern at position 10000
    ref_bytes.splice(10000..10013, b"GGGGAAAAACCCC".iter().copied());

    let fasta = format!(">chr1\n{}", std::str::from_utf8(&ref_bytes).unwrap());
    let reference = Reference::parse_fasta(&fasta).unwrap();
    let index = FMIndex::build(&reference);

    eprintln!("Reference length: {}", index.len);
    eprintln!("F-column: {:?}", index.f_column);

    // Check if pattern exists
    let pattern_bytes = seq_to_bytes("GGGGAAAAACCCC");
    eprintln!("Pattern bytes: {:?}", pattern_bytes);

    // Verify pattern in reference
    let ref_slice = reference.as_slice();
    let pattern_at_10000 = &ref_slice[10000..10013];
    eprintln!(
        "Pattern at position 10000: {:?}",
        pattern_at_10000
    );
    eprintln!("Expected pattern: {:?}", pattern_bytes);

    // Search step by step
    let mut left = 0usize;
    let mut right = index.len;

    eprintln!("\nSearch trace:");
    for (i, &c) in pattern_bytes.iter().rev().enumerate() {
        let occ_left = index.occ.occ(c, left) as usize;
        let occ_right = index.occ.occ(c, right) as usize;
        let f_c = (index.f_column[c as usize] as usize).saturating_sub(1);

        eprintln!(
            "Step {}: c={}, F-1={}, occ={}, {}, range=[{}, {})",
            i,
            c,
            f_c,
            occ_left,
            occ_right,
            f_c + occ_left,
            f_c + occ_right
        );

        left = f_c + occ_left;
        right = f_c + occ_right;

        if left >= right {
            eprintln!("  Range collapsed!");
            break;
        }
    }

    let positions = index.find_all(&pattern_bytes);
    eprintln!("\nFound positions: {:?}", positions);

    assert!(!positions.is_empty(), "Should find pattern");
}