//! Supermaximal MEM finding using binary search.
//!
//! This implementation provides O(n log m) expected time MEM finding,
//! where n is query length and m is reference length. `find_supermaximal_mems`
//! is the primary seeding path; `collect_short_seeds` adds bwa's third seeding
//! round, used only to derive the suboptimal-alignment score for MAPQ.

use crate::fm_index::FMIndex;
use crate::types::MEM;

pub const DEFAULT_MIN_MEM_LEN: usize = 19;
pub const DEFAULT_MAX_OCC: usize = 500;

/// bwa's `max_mem_intv` (third seeding round). Short exact seeds whose SA
/// interval is at most this are collected at every read position so that
/// inexact secondary loci — where only a portion of the read matches exactly —
/// surface as chainable seeds.
pub const MAX_MEM_INTV: usize = 20;

/// Re-seeding splits a long SMEM at its midpoint when the SMEM itself has very
/// few occurrences but is long enough to harbour a shorter repeated core. The
/// factor and width mirror bwa's defaults.
const SPLIT_FACTOR: f32 = 1.5;
const SPLIT_WIDTH: usize = 10;

/// Binary-search for the longest prefix of `query[start..]` (length `>= min_len`)
/// whose occurrence count is `>= min_occ`. Occurrence count is non-increasing in
/// length, so the longest qualifying prefix is the boundary. Returns 0 when even
/// the shortest prefix has fewer than `min_occ` hits.
fn longest_match_min_occ(
    index: &FMIndex,
    query: &[u8],
    start: usize,
    min_len: usize,
    min_occ: usize,
) -> usize {
    let max_possible = query.len().saturating_sub(start);
    if max_possible < min_len {
        return 0;
    }
    if index.count(&query[start..start + min_len]) < min_occ {
        return 0;
    }
    let mut left = min_len;
    let mut right = max_possible;
    while left < right {
        let mid = (left + right).div_ceil(2);
        if index.count(&query[start..start + mid]) >= min_occ {
            left = mid;
        } else {
            right = mid - 1;
        }
    }
    left
}

/// Longest prefix of `query[start..]` (length `>= min_len`) that occurs at all.
fn longest_match(index: &FMIndex, query: &[u8], start: usize, min_len: usize) -> usize {
    longest_match_min_occ(index, query, start, min_len, 1)
}

/// Shortest prefix of `query[start..]` (length `>= min_len`, `<= remaining`)
/// whose occurrence count is `<= max_intv`. Returns 0 when none qualifies.
fn shortest_match_max_occ(
    index: &FMIndex,
    query: &[u8],
    start: usize,
    min_len: usize,
    max_intv: usize,
) -> usize {
    let max_possible = query.len().saturating_sub(start);
    if max_possible < min_len {
        return 0;
    }
    // If the longest possible prefix still occurs too often, no L qualifies.
    if index.count(&query[start..start + max_possible]) > max_intv {
        return 0;
    }
    // count(min_len) is the largest count among candidates; if already small, done.
    if index.count(&query[start..start + min_len]) <= max_intv {
        return min_len;
    }
    let mut left = min_len; // count(left) > max_intv
    let mut right = max_possible; // count(right) <= max_intv
    while left + 1 < right {
        let mid = (left + right) / 2;
        if index.count(&query[start..start + mid]) <= max_intv {
            right = mid;
        } else {
            left = mid;
        }
    }
    right
}

/// bwa's third seeding round (`bwt_seed_strategy1` with `max_mem_intv`):
/// scan the read left-to-right, collecting at each position the shortest exact
/// seed (length `>= min_len`) whose occurrence count is `<= max_intv`, then
/// advancing past it. Surfaces short low-occurrence seeds (including at inexact
/// secondary loci) that the SMEM/reseed rounds miss.
///
/// These seeds are used only to derive the suboptimal alignment score for MAPQ
/// (bwa's `sub`/XS); they are intentionally kept out of the primary seeding
/// path (`find_supermaximal_mems`) so that primary placement is unaffected.
pub fn collect_short_seeds(
    index: &FMIndex,
    query: &[u8],
    min_len: usize,
    max_intv: usize,
) -> Vec<MEM> {
    let n = query.len();
    let ref_len = index.len;
    let mut seeds = Vec::new();
    let mut x = 0usize;
    while x + min_len <= n {
        let len = shortest_match_max_occ(index, query, x, min_len, max_intv);
        if len == 0 {
            x += 1;
            continue;
        }
        for &ref_start in &index.find_all_capped(&query[x..x + len], max_intv) {
            if ref_start as usize + len <= ref_len {
                seeds.push(MEM::new(x, ref_start as usize, len));
            }
        }
        x += len;
    }
    seeds
}

/// Find supermaximal MEMs using binary search, with an occurrence cap and a
/// bwa-style re-seeding pass for long unique SMEMs that span a repeated core.
///
/// `max_occ` is the maximum number of reference hits a SMEM may have; SMEMs
/// exceeding it are not enumerated for chaining but mark a repetitive span of the
/// read. Re-seeding explores the midpoint of each long, rare SMEM and appends the
/// resulting shorter MEMs to the output. Returns the MEMs together with bwa's
/// `frac_rep`: the fraction of the read covered by over-frequent (repetitive)
/// SMEMs, used to down-weight MAPQ.
pub fn find_supermaximal_mems(
    index: &FMIndex,
    query: &[u8],
    min_len: usize,
    max_occ: usize,
) -> (Vec<MEM>, f32) {
    if query.is_empty() || min_len == 0 {
        return (Vec::new(), 0.0);
    }

    let n = query.len();
    let ref_len = index.len;
    let mut raw_mems: Vec<MEM> = Vec::new();
    let mut rep_spans: Vec<(usize, usize)> = Vec::new();

    for query_start in 0..n {
        let max_possible = n - query_start;
        if max_possible < min_len {
            break;
        }

        let len = longest_match(index, query, query_start, min_len);
        if len == 0 {
            continue;
        }

        let occ = index.count(&query[query_start..query_start + len]);
        if occ > max_occ {
            rep_spans.push((query_start, query_start + len));
            continue;
        }
        for &ref_start in &index.find_all(&query[query_start..query_start + len]) {
            if ref_start as usize + len <= ref_len {
                raw_mems.push(MEM::new(query_start, ref_start as usize, len));
            }
        }
    }

    let frac_rep = repetitive_fraction(&rep_spans, n);

    // Apply supermaximal filter to the raw SMEMs.
    let filtered_smems = filter_supermaximal(raw_mems);

    // Reseeds are intentionally NOT passed through filter_supermaximal: they are
    // shorter than their parent SMEM and would be dropped as "contained". Their
    // purpose is to expose repeat loci that the long SMEM masked, so they are
    // appended to the supermaximal set rather than competing with it.
    let mut result = filtered_smems;

    // Re-seed pass: for each surviving SMEM that is long enough and rare enough,
    // split at the midpoint and look for shorter seeds covering the repeat loci.
    let n_smems = result.len();
    for i in 0..n_smems {
        let smem = result[i];
        let qs = smem.query_start;
        let len = smem.length;

        // Only re-seed when the SMEM is long relative to min_len.
        if (len as f32) < min_len as f32 * SPLIT_FACTOR {
            continue;
        }

        // Only re-seed when the SMEM itself has very few occurrences — meaning
        // the flanking context is unique but the core may be repeated.
        let smem_occ = index.count(&query[qs..qs + len]);
        if smem_occ > SPLIT_WIDTH {
            continue;
        }

        // Split at the midpoint and look for a shorter seed that is MORE frequent
        // than the parent (bwa's min_intv = parent_occ + 1): the longest match
        // spanning the midpoint whose occurrence count exceeds the SMEM's. This is
        // what exposes a repeated core buried inside an otherwise-unique SMEM.
        let mid = (smem.query_start + smem.query_end()).div_ceil(2);
        let sub_len = longest_match_min_occ(index, query, mid, min_len, smem_occ + 1);
        if sub_len == 0 {
            continue;
        }

        let sub_positions = index.find_all_capped(&query[mid..mid + sub_len], max_occ);
        for &ref_start in &sub_positions {
            if ref_start as usize + sub_len <= ref_len {
                result.push(MEM::new(mid, ref_start as usize, sub_len));
            }
        }
    }

    (result, frac_rep)
}

/// bwa's `frac_rep`: the merged query-coordinate length covered by repetitive
/// SMEMs (those exceeding `max_occ`), as a fraction of the read length. `spans`
/// must be in non-decreasing start order, matching the seeding scan.
fn repetitive_fraction(spans: &[(usize, usize)], read_len: usize) -> f32 {
    if read_len == 0 {
        return 0.0;
    }
    let (mut b, mut e, mut l_rep) = (0usize, 0usize, 0usize);
    for &(sb, se) in spans {
        if sb > e {
            l_rep += e - b;
            b = sb;
            e = se;
        } else {
            e = e.max(se);
        }
    }
    l_rep += e - b;
    (l_rep as f32 / read_len as f32).min(1.0)
}

/// Filter MEMs to keep only supermaximal ones.
/// A MEM is supermaximal if no other MEM extends further in both query and reference.
fn filter_supermaximal(mems: Vec<MEM>) -> Vec<MEM> {
    if mems.is_empty() {
        return mems;
    }

    let mut mems = mems;
    mems.sort_by_key(|b| std::cmp::Reverse(b.length));

    let mut result = Vec::new();
    let mut max_ref_end = 0usize;
    let mut max_query_end = 0usize;

    for mem in mems {
        let query_end = mem.query_end();
        let ref_end = mem.ref_end();

        // Keep if not contained in a longer MEM
        if query_end <= max_query_end && ref_end <= max_ref_end {
            continue;
        }

        result.push(mem);
        max_query_end = max_query_end.max(query_end);
        max_ref_end = max_ref_end.max(ref_end);
    }

    result
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_supermaximal_basic() {
        use crate::fm_index::FMIndex;
        use crate::reference::Reference;

        let ref_seq = Reference::parse_fasta(">test\nACGTACGT").unwrap();
        let index = FMIndex::build(&ref_seq);

        // ACGT in 2-bit encoding: A=0, C=1, G=2, T=3
        let query = [0u8, 1, 2, 3]; // ACGT
        let (mems, _) = find_supermaximal_mems(&index, &query, 2, DEFAULT_MAX_OCC);

        assert!(!mems.is_empty(), "Should find MEMs for ACGT");
        // Due to FM-index limitations, we might find shorter MEMs
        let max_len = mems.iter().map(|m| m.length).max().unwrap();
        assert!(max_len >= 2, "Should find MEMs of at least length 2");
    }

    #[test]
    fn test_supermaximal_no_match() {
        use crate::fm_index::FMIndex;
        use crate::reference::Reference;

        let ref_seq = Reference::parse_fasta(">test\nACGTACGT").unwrap();
        let index = FMIndex::build(&ref_seq);

        // TTTT (2-bit encoded as 3) doesn't exist as a contiguous run in reference
        let query = [3u8, 3, 3, 3];
        let (mems, _) = find_supermaximal_mems(&index, &query, 2, DEFAULT_MAX_OCC);

        // The algorithm will find single T matches but not a supermaximal MEM of length >= 2
        // that doesn't overlap with other MEMs
        assert!(
            mems.is_empty(),
            "Should not find supermaximal MEMs for TTTT"
        );
    }

    #[test]
    fn test_supermaximal_min_len() {
        use crate::fm_index::FMIndex;
        use crate::reference::Reference;

        let ref_seq = Reference::parse_fasta(">test\nACGTACGT").unwrap();
        let index = FMIndex::build(&ref_seq);

        let query = [0u8, 1, 2, 3]; // ACGT
        let (mems, _) = find_supermaximal_mems(&index, &query, 10, DEFAULT_MAX_OCC);

        assert!(
            mems.is_empty(),
            "Should not find MEMs with min_len > query length"
        );
    }

    #[test]
    fn test_supermaximal_empty_query() {
        use crate::fm_index::FMIndex;
        use crate::reference::Reference;

        let ref_seq = Reference::parse_fasta(">test\nACGT").unwrap();
        let index = FMIndex::build(&ref_seq);

        let (mems, _) = find_supermaximal_mems(&index, &[], 2, DEFAULT_MAX_OCC);
        assert!(mems.is_empty());
    }

    #[test]
    fn test_filter_supermaximal_contained() {
        // MEM at (query=0, ref=0, len=10) contains MEM at (query=2, ref=2, len=8)
        let mems = vec![MEM::new(0, 0, 10), MEM::new(2, 2, 8), MEM::new(5, 5, 3)];

        let filtered = filter_supermaximal(mems);
        assert_eq!(filtered.len(), 1);
        assert_eq!(filtered[0].length, 10);
    }

    #[test]
    fn test_filter_supermaximal_disjoint() {
        // Two disjoint MEMs
        let mems = vec![MEM::new(0, 0, 5), MEM::new(10, 10, 5)];

        let filtered = filter_supermaximal(mems);
        assert_eq!(filtered.len(), 2);
    }

    #[test]
    fn test_supermaximal_multiple_positions() {
        use crate::fm_index::FMIndex;
        use crate::reference::Reference;

        let ref_seq = Reference::parse_fasta(">test\nACGTACGTACGT").unwrap();
        let index = FMIndex::build(&ref_seq);

        // Query contains substring that appears in reference
        let query = [0u8, 1, 2, 3]; // ACGT
        let (mems, _) = find_supermaximal_mems(&index, &query, 2, DEFAULT_MAX_OCC);

        // Should find at least some MEMs
        assert!(!mems.is_empty(), "Should find MEMs for ACGT");
        let positions: Vec<_> = mems.iter().map(|m| m.ref_start).collect();
        // Verify at least one valid position was found
        assert!(!positions.is_empty(), "Should have found valid positions");
    }

    /// Verify that re-seeding surfaces repeat loci hidden by a long unique SMEM.
    ///
    /// Reference layout (2-bit encoded, A=0 C=1 G=2 T=3):
    ///   [unique flank A×6] [repeat core ACGT×3 = 12 bp] [unique flank T×6]
    ///   [gap CCCC = 4 bp]
    ///   [unique flank G×6] [repeat core ACGT×3 = 12 bp] [unique flank A×6]
    ///
    /// The repeat core (ACGTACGTACGT, 12 bp) appears at two loci. The unique
    /// flanks on either side extend each copy into a long SMEM covering the
    /// full segment.  Without re-seeding only one chain exists; with re-seeding
    /// the midpoint of the SMEM yields seeds from both loci.
    #[test]
    fn test_reseed_surfaces_repeat_loci() {
        use crate::fm_index::FMIndex;
        use crate::reference::Reference;

        // Build reference in 2-bit encoding as FASTA (A C G T characters)
        // Locus 1: AAAAAAA + ACGTACGTACGT + TTTTTTT  (6+12+6 = 24 bp)
        // Separator: CCCC (4 bp)
        // Locus 2: GGGGGGG + ACGTACGTACGT + AAAAAAA  (6+12+6 = 24 bp)
        let fasta = ">t\n\
            AAAAAACGTACGTACGTTTTTT\
            CCCC\
            GGGGGGACGTACGTACGTAAAAAA";
        let ref_seq = Reference::parse_fasta(fasta).unwrap();
        let index = FMIndex::build(&ref_seq);

        // The shared core in 2-bit: ACGTACGTACGT
        // Embed it in a query that also contains flanking unique context so the
        // full query forms one SMEM over the first locus:
        //   AAAAAA (6) + ACGTACGTACGT (12) = 18 bp query
        // With min_len=4 and SPLIT_FACTOR=1.5, a 12-bp SMEM triggers a split.
        // Encode manually: A=0, C=1, G=2, T=3
        let query: Vec<u8> = "AAAAAACGTACGTACGT"
            .bytes()
            .map(|b| match b {
                b'A' => 0,
                b'C' => 1,
                b'G' => 2,
                b'T' => 3,
                _ => 4,
            })
            .collect();

        // min_len=4, SPLIT_FACTOR=1.5: a SMEM of length >= 6 with <= SPLIT_WIDTH
        // occurrences is split at its midpoint, surfacing the repeated core.
        let min_len = 4;
        let (mems, _) = find_supermaximal_mems(&index, &query, min_len, DEFAULT_MAX_OCC);

        // The repeat core occurs at two reference loci; seeding must surface both.
        let distinct_ref_starts: std::collections::HashSet<usize> =
            mems.iter().map(|m| m.ref_start).collect();
        assert!(
            distinct_ref_starts.len() >= 2,
            "Seeding should expose >= 2 distinct ref loci; got {:?} from MEMs {:?}",
            distinct_ref_starts,
            mems
        );
    }

    /// Verify `shortest_match_max_occ` returns the correct boundary length and 0
    /// when nothing qualifies.
    ///
    /// Reference: "ACGTACGT" (8 bp, 2-bit encoded).
    ///   - "A"       occurs 2 times  (> max_intv=1)
    ///   - "AC"      occurs 2 times  (> max_intv=1)
    ///   - "ACG"     occurs 2 times  (> max_intv=1)
    ///   - "ACGT"    occurs 2 times  (> max_intv=1)
    ///   - "ACGTA"   occurs 1 time   (<= max_intv=1)  ← boundary
    ///
    /// With min_len=1 and max_intv=1, the shortest qualifying prefix from pos 0
    /// should be length 5 ("ACGTA").
    ///
    /// With max_intv=0 nothing qualifies, so 0 is returned.
    #[test]
    fn test_shortest_match_max_occ_boundary() {
        use crate::fm_index::FMIndex;
        use crate::reference::Reference;

        // ACGTACGT: "ACGT" appears twice, "ACGTA" appears once.
        let ref_seq = Reference::parse_fasta(">test\nACGTACGT").unwrap();
        let index = FMIndex::build(&ref_seq);

        // With max_intv=1, the shortest prefix from position 0 of length >= 1
        // whose count is <= 1 is "ACGTA" (length 5).
        let len = shortest_match_max_occ(&index, &[0u8, 1, 2, 3, 0, 1, 2, 3], 0, 1, 1);
        assert_eq!(
            len, 5,
            "shortest prefix with count<=1 should be length 5 (ACGTA)"
        );

        // With max_intv=0 nothing qualifies (every non-empty prefix occurs >= 1 time).
        let len_none = shortest_match_max_occ(&index, &[0u8, 1, 2, 3, 0, 1, 2, 3], 0, 1, 0);
        assert_eq!(len_none, 0, "nothing should qualify when max_intv=0");

        // With min_len larger than what remains, should also return 0.
        let len_short = shortest_match_max_occ(&index, &[0u8, 1, 2, 3], 0, 10, 1);
        assert_eq!(
            len_short, 0,
            "should return 0 when fewer than min_len bases remain"
        );
    }

    /// Verify that the third seeding round surfaces a secondary locus that shares
    /// only a short exact core with the query.
    ///
    /// Reference (2-bit encoded, A=0 C=1 G=2 T=3):
    ///   `TTTTTTACGTACGT` `CCCC` `GGGGGGACGTACGT`
    ///    locus 1 (0..14)  sep    locus 2 (18..32)
    ///
    /// The query `TTTTTTACGTACGT` matches locus 1 in full but shares only the
    /// `ACGTACGT` core (8 bp, occurs twice) with locus 2. The third round scans
    /// every position for short low-occurrence matches and must surface a seed in
    /// locus 2 (ref offset >= 18) covering that core.
    #[test]
    fn test_third_round_surfaces_inexact_secondary() {
        use crate::fm_index::FMIndex;
        use crate::reference::Reference;

        // Reference: TTTTTTACGTACGT + CCCC + GGGGGGACGTACGT
        //            0..13           14..17  18..31
        let fasta = ">t\nTTTTTTACGTACGTCCCCGGGGGGACGTACGT";
        let ref_seq = Reference::parse_fasta(fasta).unwrap();
        let index = FMIndex::build(&ref_seq);

        // Query = TTTTTTACGTACGT (14 bp): matches locus 1 exactly.
        // Locus 2 shares only "ACGTACGT" (8 bp) starting at ref offset 24 (18+6).
        let encode = |s: &str| -> Vec<u8> {
            s.bytes()
                .map(|b| match b {
                    b'A' => 0,
                    b'C' => 1,
                    b'G' => 2,
                    b'T' => 3,
                    _ => 4,
                })
                .collect()
        };
        let query = encode("TTTTTTACGTACGT");

        // Round 1 (SMEM) emits the unique full-query match at locus 1; the
        // midpoint reseed (round 2) is not guaranteed to land on locus 2's
        // shared core. The third round (`collect_short_seeds`) scans every
        // position and finds the 8-mer "ACGTACGT" (count=2 <= MAX_MEM_INTV) at
        // BOTH loci, so it must surface a seed at locus 2.
        let min_len = 4;
        let (primary_mems, _) = find_supermaximal_mems(&index, &query, min_len, DEFAULT_MAX_OCC);
        let short_seeds = collect_short_seeds(&index, &query, min_len, MAX_MEM_INTV);

        // Locus 2 occupies ref offsets 18..32 (GGGGGG at 18, ACGTACGT at 24).
        // The third round scans short low-occurrence matches across the read and
        // must surface at least one seed inside locus 2 covering the shared core.
        assert!(
            short_seeds.iter().any(|m| m.ref_start >= 18),
            "third round must surface a seed in locus 2 (ref_start >= 18); got: {:?}",
            short_seeds
        );

        // Combined, the candidate set covers >= 2 distinct ref loci.
        let distinct: std::collections::HashSet<usize> = primary_mems
            .iter()
            .chain(short_seeds.iter())
            .map(|m| m.ref_start)
            .collect();
        assert!(
            distinct.len() >= 2,
            "Should expose >= 2 distinct ref loci; got {:?}",
            distinct
        );
    }
}
