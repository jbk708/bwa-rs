//! Paired-end read handling and pairing logic.

use crate::types::{AlignmentResult, Orientation};

/// Std-dev multiplier for the proper-pair insert-size window (bwa `mem_pestat` MAX_STDDEV).
const MAX_STDDEV: f64 = 4.0;

/// Fully-assembled per-read mate fields for SAM output.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct MateFields {
    /// Fully assembled SAM flag for this read.
    pub flag: u16,
    /// RNEXT field: "=" when read and mate share reference coordinates, "*" otherwise.
    pub rnext: &'static str,
    /// 1-based mate POS to print; 0 when no coordinate is available.
    pub pnext: i64,
    /// Signed template length per bwa convention.
    pub tlen: i64,
    /// When this read is unmapped and the mate is mapped, the mate's 0-based position
    /// is placed here so the CLI can use it as this read's POS.
    pub placed_pos: Option<usize>,
    /// Mate's CIGAR in SAM M-form (MC:Z), present only when the mate is mapped.
    pub mc: Option<String>,
}

#[derive(Clone, Debug, Default)]
pub struct InsertSizeDistribution {
    pub mean: f64,
    pub std_dev: f64,
    count: u64,
    sum: f64,
    sum_sq: f64,
}

impl InsertSizeDistribution {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn with_params(mean: f64, std_dev: f64) -> Self {
        Self {
            mean,
            std_dev,
            count: 0,
            sum: 0.0,
            sum_sq: 0.0,
        }
    }

    pub fn add(&mut self, insert_size: u32) {
        let size = insert_size as f64;
        self.count += 1;
        self.sum += size;
        self.sum_sq += size * size;

        if self.count == 1 {
            self.mean = size;
            self.std_dev = 0.0;
        } else {
            let n = self.count as f64;
            let new_mean = self.mean + (size - self.mean) / n;
            self.std_dev = ((n - 1.0) * self.std_dev * self.std_dev
                + (size - self.mean) * (size - new_mean))
                / n;
            self.std_dev = self.std_dev.sqrt();
            self.mean = new_mean;
        }
    }

    pub fn is_anomalous(&self, insert_size: u32) -> bool {
        let expected_max = self.mean + 3.0 * self.std_dev;
        let expected_min = (self.mean - 3.0 * self.std_dev).max(0.0);
        (insert_size as f64) > expected_max || (insert_size as f64) < expected_min
    }

    pub fn lower_bound(&self) -> u32 {
        (self.mean - MAX_STDDEV * self.std_dev).max(0.0) as u32
    }

    pub fn upper_bound(&self) -> u32 {
        (self.mean + MAX_STDDEV * self.std_dev) as u32
    }
}

#[derive(Clone, Debug)]
pub struct PairedResult {
    pub read1: AlignmentResult,
    pub read2: AlignmentResult,
    pub orientation: Orientation,
    pub insert_size: i32,
    pub proper_pair: bool,
}

impl PairedResult {
    pub fn both_mapped(&self) -> bool {
        (self.read1.flag & 0x4) == 0 && (self.read2.flag & 0x4) == 0
    }

    pub fn opposite_strands(&self) -> bool {
        self.read1.reverse_strand != self.read2.reverse_strand
    }
}

pub fn pair_reads(
    read1: AlignmentResult,
    read2: AlignmentResult,
    insert_dist: &InsertSizeDistribution,
) -> PairedResult {
    let mut paired = PairedResult {
        read1: read1.clone(),
        read2: read2.clone(),
        orientation: Orientation::FR,
        insert_size: 0,
        proper_pair: false,
    };

    if !read1.reverse_strand && read2.reverse_strand {
        paired.orientation = Orientation::FR;
    } else if read1.reverse_strand && !read2.reverse_strand {
        paired.orientation = Orientation::RF;
    } else if !read1.reverse_strand && !read2.reverse_strand {
        paired.orientation = Orientation::FF;
    } else {
        paired.orientation = Orientation::RR;
    }

    if paired.both_mapped() {
        let r1_pos = read1.position as i32;
        let r2_pos = read2.position as i32;

        let min_pos = r1_pos.min(r2_pos);
        let max_pos = r1_pos.max(r2_pos);

        paired.insert_size = max_pos - min_pos + 1;
        paired.proper_pair = paired.opposite_strands()
            && (paired.insert_size as u32) >= insert_dist.lower_bound()
            && (paired.insert_size as u32) <= insert_dist.upper_bound();
    }

    paired.read1.flag |= 0x1;
    paired.read2.flag |= 0x1;

    if paired.proper_pair {
        paired.read1.flag |= 0x2;
        paired.read2.flag |= 0x2;
    }

    paired
}

/// Insert size of an FR-oriented pair measured between the two reads' 5′ ends, or
/// `None` when the pair is not FR-concordant.
///
/// Each read is anchored at its 5′-most reference coordinate (the forward read's start,
/// the reverse read's end), matching the [`template_length`] convention. The pair is FR
/// when the forward read's 5′ end lies at or before the reverse read's 5′ end, which
/// tolerates the small dovetail overlaps bwa still treats as FR. Returns `None` if either
/// read is unmapped, the strands match, or the orientation is RF/outward.
fn fr_insert_size(read1: &AlignmentResult, read2: &AlignmentResult) -> Option<u32> {
    if read1.is_unmapped() || read2.is_unmapped() {
        return None;
    }
    if read1.reverse_strand == read2.reverse_strand {
        return None;
    }
    let (fwd, rev) = if !read1.reverse_strand {
        (read1, read2)
    } else {
        (read2, read1)
    };
    let fwd_5p = fwd.position as i64;
    let rev_5p = rev.position as i64 + rev.cigar.reference_length() as i64 - 1;
    if fwd_5p > rev_5p {
        return None;
    }
    Some((rev_5p - fwd_5p + 1) as u32)
}

/// Returns true when both reads form a proper FR pair within the insert-size bounds.
///
/// Requires FR orientation (see [`fr_insert_size`]) and a 5′-to-5′ insert size within
/// `[dist.lower_bound(), dist.upper_bound()]`.
pub fn is_proper_pair(
    read1: &AlignmentResult,
    read2: &AlignmentResult,
    dist: &InsertSizeDistribution,
) -> bool {
    match fr_insert_size(read1, read2) {
        Some(isize) => isize >= dist.lower_bound() && isize <= dist.upper_bound(),
        None => false,
    }
}

/// Computes the signed template length for `read` given `mate`, using the bwa convention.
///
/// Returns 0 when either read is unmapped.
pub fn template_length(read: &AlignmentResult, mate: &AlignmentResult) -> i64 {
    if read.is_unmapped() || mate.is_unmapped() {
        return 0;
    }
    let p0 = read.position as i64
        + if read.reverse_strand {
            (read.cigar.reference_length() as i64) - 1
        } else {
            0
        };
    let p1 = mate.position as i64
        + if mate.reverse_strand {
            (mate.cigar.reference_length() as i64) - 1
        } else {
            0
        };
    if p0 < p1 {
        p1 - p0 + 1
    } else if p0 > p1 {
        -(p0 - p1 + 1)
    } else {
        0
    }
}

/// Returns the FR-oriented 5′-to-5′ insert size for feeding the insert-size distribution.
///
/// Delegates to [`fr_insert_size`]: `Some` only for mapped, opposite-strand, FR-concordant
/// pairs, so the distribution is estimated from concordant pairs alone (matching bwa's
/// per-orientation binning).
pub fn pair_insert_size(read1: &AlignmentResult, read2: &AlignmentResult) -> Option<u32> {
    fr_insert_size(read1, read2)
}

/// Assembles the SAM mate fields for one read of a pair.
///
/// `is_read1` distinguishes read1 (0x40) from read2 (0x80).
/// `proper_pair` should be pre-computed via [`is_proper_pair`].
pub fn mate_fields(
    read: &AlignmentResult,
    mate: &AlignmentResult,
    is_read1: bool,
    proper_pair: bool,
) -> MateFields {
    let read_unmapped = read.is_unmapped();
    let mate_unmapped = mate.is_unmapped();

    let mut flag: u16 = 0x1;
    if read_unmapped {
        flag |= 0x4;
    } else if read.reverse_strand {
        flag |= 0x10;
    }
    if is_read1 {
        flag |= 0x40;
    } else {
        flag |= 0x80;
    }
    if mate_unmapped {
        flag |= 0x8;
    } else if mate.reverse_strand {
        flag |= 0x20;
    }
    if proper_pair {
        flag |= 0x2;
    }

    let tlen = template_length(read, mate);

    let (rnext, pnext) = if !mate_unmapped {
        ("=", mate.position as i64 + 1)
    } else if !read_unmapped {
        ("=", read.position as i64 + 1)
    } else {
        ("*", 0)
    };

    let placed_pos = if read_unmapped && !mate_unmapped {
        Some(mate.position)
    } else {
        None
    };

    let mc = if mate_unmapped {
        None
    } else {
        Some(mate.cigar.to_sam_string())
    };

    MateFields {
        flag,
        rnext,
        pnext,
        tlen,
        placed_pos,
        mc,
    }
}

pub fn rescue_orphan(
    orphan: &AlignmentResult,
    mate_pos: usize,
    insert_dist: &InsertSizeDistribution,
) -> Option<AlignmentResult> {
    let distance = (orphan.position as i32 - mate_pos as i32).unsigned_abs();

    if distance <= insert_dist.upper_bound() {
        let mut rescued = orphan.clone();
        rescued.flag |= 0x8;
        return Some(rescued);
    }

    None
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::types::Cigar;

    #[test]
    fn test_insert_distribution() {
        let mut dist = InsertSizeDistribution::new();

        dist.add(200);
        dist.add(210);
        dist.add(190);
        dist.add(205);

        assert!((dist.mean - 201.25).abs() < 0.1);
        assert!(dist.std_dev > 0.0);
    }

    #[test]
    fn test_anomalous_pair() {
        let dist = InsertSizeDistribution::with_params(200.0, 30.0);

        assert!(dist.is_anomalous(1000));
        assert!(!dist.is_anomalous(200));
    }

    #[test]
    fn test_orientation() {
        let mut read1 = AlignmentResult::new(100, Cigar::new());
        read1.reverse_strand = false;
        let mut read2 = AlignmentResult::new(200, Cigar::new());
        read2.reverse_strand = true;

        let dist = InsertSizeDistribution::with_params(200.0, 30.0);
        let paired = pair_reads(read1, read2, &dist);

        assert_eq!(paired.orientation, Orientation::FR);
    }

    fn mapped_read(position: usize, ref_len: u32, reverse: bool) -> AlignmentResult {
        use crate::types::CigarOp;
        let mut cigar = Cigar::new();
        cigar.push(CigarOp::M, ref_len);
        let mut r = AlignmentResult::new(position, cigar);
        r.reverse_strand = reverse;
        r
    }

    fn unmapped_read() -> AlignmentResult {
        let mut r = AlignmentResult::new(0, Cigar::new());
        r.flag = 0x4;
        r
    }

    // --- is_proper_pair ---

    #[test]
    fn proper_pair_fr_within_bounds() {
        let r1 = mapped_read(100, 100, false);
        let r2 = mapped_read(200, 100, true);
        let dist = InsertSizeDistribution::with_params(250.0, 50.0);
        assert!(is_proper_pair(&r1, &r2, &dist));
    }

    #[test]
    fn proper_pair_fr_out_of_bounds() {
        let r1 = mapped_read(100, 100, false);
        let r2 = mapped_read(5000, 100, true);
        let dist = InsertSizeDistribution::with_params(250.0, 50.0);
        assert!(!is_proper_pair(&r1, &r2, &dist));
    }

    #[test]
    fn proper_pair_same_strand_not_proper() {
        let r1 = mapped_read(100, 100, false);
        let r2 = mapped_read(200, 100, false);
        let dist = InsertSizeDistribution::with_params(250.0, 50.0);
        assert!(!is_proper_pair(&r1, &r2, &dist));
    }

    #[test]
    fn proper_pair_rf_not_proper() {
        // forward read is to the RIGHT of the reverse read → not FR orientation
        let r1 = mapped_read(100, 100, true); // reverse at 100
        let r2 = mapped_read(300, 100, false); // forward at 300 (forward > reverse → RF)
        let dist = InsertSizeDistribution::with_params(400.0, 100.0);
        assert!(!is_proper_pair(&r1, &r2, &dist));
    }

    #[test]
    fn proper_pair_one_unmapped() {
        let r1 = mapped_read(100, 100, false);
        let r2 = unmapped_read();
        let dist = InsertSizeDistribution::with_params(250.0, 50.0);
        assert!(!is_proper_pair(&r1, &r2, &dist));
    }

    // --- template_length ---

    #[test]
    fn template_length_fr_pair() {
        // read1 at 100 (fwd), read2 at 300 (rev, ref_len=100)
        // p0 = 100 (fwd, so +0), p1 = 300+99 = 399 (rev)
        // tlen(r1,r2) = 399-100+1 = 300
        let r1 = mapped_read(100, 100, false);
        let r2 = mapped_read(300, 100, true);
        assert_eq!(template_length(&r1, &r2), 300);
        assert_eq!(template_length(&r2, &r1), -300);
    }

    #[test]
    fn template_length_antisymmetry() {
        let r1 = mapped_read(50, 75, false);
        let r2 = mapped_read(200, 75, true);
        let tlen_r1 = template_length(&r1, &r2);
        let tlen_r2 = template_length(&r2, &r1);
        assert_eq!(tlen_r1, -tlen_r2);
        assert!(tlen_r1 > 0);
    }

    #[test]
    fn template_length_equal_anchor_zero() {
        // Both reads start at same position, same ref_len, both forward
        // p0 = 0, p1 = 0 → tlen = 0
        let r1 = mapped_read(0, 50, false);
        let r2 = mapped_read(0, 50, false);
        assert_eq!(template_length(&r1, &r2), 0);
    }

    #[test]
    fn template_length_either_unmapped_zero() {
        let r1 = mapped_read(100, 100, false);
        let r2 = unmapped_read();
        assert_eq!(template_length(&r1, &r2), 0);
        assert_eq!(template_length(&r2, &r1), 0);
    }

    // --- pair_insert_size ---

    #[test]
    fn pair_insert_size_fr_pair() {
        // r1: [100, 200), r2: [200, 300) → span = 300 - 100 = 200
        let r1 = mapped_read(100, 100, false);
        let r2 = mapped_read(200, 100, true);
        assert_eq!(pair_insert_size(&r1, &r2), Some(200));
    }

    #[test]
    fn pair_insert_size_none_when_unmapped() {
        let r1 = mapped_read(100, 100, false);
        let r2 = unmapped_read();
        assert_eq!(pair_insert_size(&r1, &r2), None);
        assert_eq!(pair_insert_size(&r2, &r1), None);
    }

    #[test]
    fn proper_pair_dovetail_fr() {
        // Forward read starts a few bp AFTER the reverse read's start (dovetail overlap),
        // so the alignment-start ordering is inverted, but the 5' ends are still FR.
        let fwd = mapped_read(102, 100, false); // 5' = 102
        let rev = mapped_read(100, 100, true); // 5' = 100 + 99 = 199
        assert_eq!(fr_insert_size(&fwd, &rev), Some(98));
        let dist = InsertSizeDistribution::with_params(250.0, 50.0);
        assert!(
            is_proper_pair(&fwd, &rev, &dist),
            "small dovetail is still a proper FR pair"
        );
    }

    #[test]
    fn proper_pair_large_insert_within_4sigma() {
        // Insert 800 is beyond mean+3*std (710) but within mean+4*std (830): bwa MAX_STDDEV.
        let r1 = mapped_read(100, 100, false);
        let r2 = mapped_read(800, 100, true); // 5' insert = (800+99) - 100 + 1 = 800
        assert_eq!(fr_insert_size(&r1, &r2), Some(800));
        let dist = InsertSizeDistribution::with_params(350.0, 120.0);
        assert!(
            is_proper_pair(&r1, &r2, &dist),
            "insert within 4 std-devs is proper"
        );
    }

    // --- mate_fields flag bits ---

    #[test]
    fn mate_fields_both_mapped_proper_fr() {
        let r1 = mapped_read(100, 100, false);
        let r2 = mapped_read(200, 100, true);
        let dist = InsertSizeDistribution::with_params(250.0, 50.0);
        let proper = is_proper_pair(&r1, &r2, &dist);
        assert!(proper);

        let mf1 = mate_fields(&r1, &r2, true, proper);
        assert_eq!(mf1.flag & 0x1, 0x1, "paired bit always set");
        assert_eq!(mf1.flag & 0x40, 0x40, "read1 bit set");
        assert_eq!(mf1.flag & 0x80, 0, "read2 bit not set for read1");
        assert_eq!(mf1.flag & 0x2, 0x2, "proper pair set");
        assert_eq!(mf1.flag & 0x4, 0, "read1 is mapped");
        assert_eq!(mf1.flag & 0x8, 0, "mate (r2) is mapped");
        assert_eq!(mf1.flag & 0x20, 0x20, "mate reverse set (r2 is reverse)");
        assert_eq!(mf1.flag & 0x10, 0, "read1 not reverse");
        assert_eq!(mf1.rnext, "=");
        assert_eq!(mf1.pnext, 201);
        assert_eq!(mf1.placed_pos, None);
        assert_eq!(
            mf1.mc,
            Some("100M".to_string()),
            "r2 is mapped so mc is Some"
        );

        let mf2 = mate_fields(&r2, &r1, false, proper);
        assert_eq!(mf2.flag & 0x1, 0x1, "paired bit always set");
        assert_eq!(mf2.flag & 0x80, 0x80, "read2 bit set");
        assert_eq!(mf2.flag & 0x40, 0, "read1 bit not set for read2");
        assert_eq!(mf2.flag & 0x2, 0x2, "proper pair set");
        assert_eq!(mf2.flag & 0x4, 0, "read2 is mapped");
        assert_eq!(mf2.flag & 0x8, 0, "mate (r1) is mapped");
        assert_eq!(mf2.flag & 0x20, 0, "mate (r1) not reverse");
        assert_eq!(mf2.flag & 0x10, 0x10, "read2 is reverse");
        assert_eq!(mf2.rnext, "=");
        assert_eq!(mf2.pnext, 101);
        assert_eq!(mf2.placed_pos, None);
        assert_eq!(
            mf2.mc,
            Some("100M".to_string()),
            "r1 is mapped so mc is Some"
        );
    }

    #[test]
    fn mate_fields_no_proper_when_out_of_bounds() {
        let r1 = mapped_read(100, 100, false);
        let r2 = mapped_read(5000, 100, true);
        let dist = InsertSizeDistribution::with_params(250.0, 50.0);
        let proper = is_proper_pair(&r1, &r2, &dist);
        assert!(!proper);

        let mf1 = mate_fields(&r1, &r2, true, proper);
        assert_eq!(mf1.flag & 0x2, 0, "proper pair bit should not be set");
    }

    #[test]
    fn mate_fields_no_proper_same_strand() {
        let r1 = mapped_read(100, 100, false);
        let r2 = mapped_read(200, 100, false);
        let dist = InsertSizeDistribution::with_params(250.0, 50.0);
        let proper = is_proper_pair(&r1, &r2, &dist);
        assert!(!proper);

        let mf1 = mate_fields(&r1, &r2, true, proper);
        assert_eq!(mf1.flag & 0x2, 0, "proper pair bit not set for same-strand");
        assert_eq!(mf1.flag & 0x20, 0, "mate forward so 0x20 not set");
    }

    #[test]
    fn mate_fields_mate_unmapped() {
        let r1 = mapped_read(100, 100, false);
        let r2 = unmapped_read();

        let mf1 = mate_fields(&r1, &r2, true, false);
        assert_eq!(mf1.flag & 0x8, 0x8, "mate unmapped bit set");
        assert_eq!(mf1.flag & 0x20, 0, "0x20 not set when mate unmapped");
        assert_eq!(mf1.rnext, "=", "unmapped mate placed at read's coord");
        assert_eq!(mf1.pnext, 101, "pnext = read's 1-based pos");
        assert_eq!(mf1.placed_pos, None, "mapped read has no placed_pos");
        assert_eq!(mf1.mc, None, "mate is unmapped so mc is None");

        let mf2 = mate_fields(&r2, &r1, false, false);
        assert_eq!(mf2.flag & 0x4, 0x4, "this read is unmapped");
        assert_eq!(mf2.flag & 0x8, 0, "mate (r1) is mapped");
        assert_eq!(mf2.rnext, "=", "placed at mate's coord");
        assert_eq!(mf2.pnext, 101, "pnext = mate's 1-based pos");
        assert_eq!(
            mf2.placed_pos,
            Some(100),
            "unmapped read inherits mate's 0-based pos"
        );
        assert_eq!(
            mf2.mc,
            Some("100M".to_string()),
            "mate (r1) is mapped so mc is Some"
        );
    }

    #[test]
    fn mate_fields_both_unmapped() {
        let r1 = unmapped_read();
        let r2 = unmapped_read();

        let mf1 = mate_fields(&r1, &r2, true, false);
        assert_eq!(mf1.flag & 0x4, 0x4, "this read unmapped");
        assert_eq!(mf1.flag & 0x8, 0x8, "mate unmapped");
        assert_eq!(mf1.rnext, "*");
        assert_eq!(mf1.pnext, 0);
        assert_eq!(mf1.tlen, 0);
        assert_eq!(mf1.placed_pos, None);
        assert_eq!(mf1.mc, None, "both unmapped so mc is None");
    }

    #[test]
    fn mate_fields_tlen_sign_fr_pair() {
        // r1 fwd at 100 (ref_len 100), r2 rev at 300 (ref_len 100)
        // template_length(r1,r2) > 0 (r1 is left anchor)
        // template_length(r2,r1) < 0
        let r1 = mapped_read(100, 100, false);
        let r2 = mapped_read(300, 100, true);
        let proper = true;

        let mf1 = mate_fields(&r1, &r2, true, proper);
        let mf2 = mate_fields(&r2, &r1, false, proper);
        assert!(mf1.tlen > 0, "left read tlen positive");
        assert!(mf2.tlen < 0, "right read tlen negative");
        assert_eq!(mf1.tlen, -mf2.tlen);
    }

    #[test]
    fn mate_fields_read1_read2_bits_exclusive() {
        let r1 = mapped_read(100, 100, false);
        let r2 = mapped_read(200, 100, true);

        let mf_as_r1 = mate_fields(&r1, &r2, true, false);
        let mf_as_r2 = mate_fields(&r1, &r2, false, false);

        assert_ne!(mf_as_r1.flag & 0x40, 0, "0x40 set when is_read1=true");
        assert_eq!(mf_as_r1.flag & 0x80, 0, "0x80 not set when is_read1=true");
        assert_eq!(mf_as_r2.flag & 0x40, 0, "0x40 not set when is_read1=false");
        assert_ne!(mf_as_r2.flag & 0x80, 0, "0x80 set when is_read1=false");
    }
}
