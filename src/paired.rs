//! Paired-end read handling and pairing logic.

use crate::types::{AlignmentResult, Orientation};

/// Std-dev multiplier for the proper-pair insert-size window (bwa `mem_pestat` MAX_STDDEV).
const MAX_STDDEV: f64 = 4.0;
/// IQR multiplier used to reject outliers when computing mean/std_dev.
const OUTLIER_BOUND: f64 = 2.0;
/// IQR multiplier used to compute the proper-pair bounds before widening toward mean±MAX_STDDEV*σ.
const MAPPING_BOUND: f64 = 3.0;
/// Minimum number of concordant insert sizes required to produce a valid distribution.
const MIN_DIR_CNT: usize = 10;

/// bwa `-U` default unpaired-read penalty subtracted from the unpaired score when
/// choosing between paired and unpaired placements.
const PEN_UNPAIRED: i32 = 17;

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
    low: u32,
    high: u32,
}

impl InsertSizeDistribution {
    pub fn new() -> Self {
        Self::default()
    }

    /// Construct directly from a known mean/std_dev, deriving the proper-pair bounds
    /// from `mean ± MAX_STDDEV * std_dev` alone. Unlike [`Self::from_insert_sizes`] this
    /// has no percentile data, so it omits the IQR widening; it is intended for tests and
    /// callers that supply distribution parameters rather than raw observations.
    pub fn with_params(mean: f64, std_dev: f64) -> Self {
        let low = (mean - MAX_STDDEV * std_dev).max(0.0) as u32;
        let high = (mean + MAX_STDDEV * std_dev) as u32;
        Self {
            mean,
            std_dev,
            low,
            high,
        }
    }

    /// Estimate the insert-size distribution from a slice of FR-concordant insert sizes
    /// using bwa's `mem_pestat` robust algorithm: IQR-based outlier filtering for
    /// mean/std_dev, then IQR-based proper-pair bounds widened toward mean±MAX_STDDEV*σ.
    pub fn from_insert_sizes(sizes: &[u32]) -> Self {
        if sizes.len() < MIN_DIR_CNT {
            return Self::default();
        }

        let n = sizes.len();
        let mut sorted = sizes.to_vec();
        sorted.sort_unstable();

        let percentile = |q: f64| sorted[((q * n as f64 + 0.499) as usize).min(n - 1)];
        let p25 = percentile(0.25);
        let p75 = percentile(0.75);
        let iqr = p75 as f64 - p25 as f64;

        let filter_low = ((p25 as f64 - OUTLIER_BOUND * iqr + 0.499) as i64).max(1);
        let filter_high = (p75 as f64 + OUTLIER_BOUND * iqr + 0.499) as i64;

        let inliers: Vec<f64> = sizes
            .iter()
            .filter(|&&x| {
                let v = x as i64;
                v >= filter_low && v <= filter_high
            })
            .map(|&x| x as f64)
            .collect();

        if inliers.is_empty() {
            return Self::default();
        }

        let mean = inliers.iter().sum::<f64>() / inliers.len() as f64;
        let variance = inliers
            .iter()
            .map(|&x| (x - mean) * (x - mean))
            .sum::<f64>()
            / inliers.len() as f64;
        let std_dev = variance.sqrt();

        let sigma_low = mean - MAX_STDDEV * std_dev;
        let sigma_high = mean + MAX_STDDEV * std_dev;

        let mut low = ((p25 as f64 - MAPPING_BOUND * iqr + 0.499) as i64).max(1);
        let mut high = (p75 as f64 + MAPPING_BOUND * iqr + 0.499) as i64;

        if (low as f64) > sigma_low {
            low = (sigma_low + 0.499) as i64;
        }
        if (high as f64) < sigma_high {
            high = (sigma_high + 0.499) as i64;
        }

        let low = low.max(1) as u32;
        let high = high.max(1) as u32;

        Self {
            mean,
            std_dev,
            low,
            high,
        }
    }

    pub fn is_anomalous(&self, insert_size: u32) -> bool {
        let expected_max = self.mean + 3.0 * self.std_dev;
        let expected_min = (self.mean - 3.0 * self.std_dev).max(0.0);
        (insert_size as f64) > expected_max || (insert_size as f64) < expected_min
    }

    pub fn lower_bound(&self) -> u32 {
        self.low
    }

    pub fn upper_bound(&self) -> u32 {
        self.high
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

/// bwa's raw_mapq: linear in score difference, scaled by the match score `a`.
fn raw_mapq(diff: i32, a: i32) -> i32 {
    (6.02 * diff as f64 / a as f64 + 0.499) as i32
}

/// Abramowitz & Stegun 7.1.26 rational approximation for erfc.
///
/// Max absolute error ~1.5e-7. Handles negative inputs via erfc(-x) = 2 - erfc(x).
fn erfc(x: f64) -> f64 {
    let sign = if x < 0.0 { -1.0_f64 } else { 1.0_f64 };
    let x = x.abs();
    let t = 1.0 / (1.0 + 0.3275911 * x);
    let y = 1.0
        - (((((1.061405429 * t - 1.453152027) * t) + 1.421413741) * t - 0.284496736) * t
            + 0.254829592)
            * t
            * (-x * x).exp();
    1.0 - sign * y
}

/// bwa's insert-size log-likelihood bonus (mem_pair). Zero when the
/// distribution has no spread.
fn insert_bonus(insert_size: u32, dist: &InsertSizeDistribution, a: i32) -> i32 {
    if dist.std_dev <= 0.0 {
        return 0;
    }
    let ns = (insert_size as f64 - dist.mean) / dist.std_dev;
    let arg = ns.abs() * std::f64::consts::FRAC_1_SQRT_2;
    let bonus = 0.721 * (2.0 * erfc(arg)).ln() * a as f64 + 0.499;
    bonus as i32
}

/// Port of bwa-mem2 `mem_sam_pe`'s paired-MAPQ recalculation for a pair with one
/// candidate placement per end. Returns the adjusted (read1_mapq, read2_mapq).
/// `match_score` is bwa's `opt->a`. Inputs come from each mate's AlignmentResult
/// (score, mapq as the SE estimate, frac_rep). The pair must be FR-concordant
/// (call only when is_proper_pair && both mapped); insert size is taken from
/// fr_insert_size.
///
/// bwa's actual gate is the pairing score `o > 0`; with both mates filtered above
/// the `-T` threshold a proper FR pair always clears it, so callers gate on
/// [`is_proper_pair`]. With a single candidate per end there are no alternative
/// pairings, so the suboptimal-pairing count `n_sub` is 0 and bwa's
/// `4.343 * ln(n_sub + 1)` multi-hit penalty drops out.
pub fn pair_mapq(
    read1: &AlignmentResult,
    read2: &AlignmentResult,
    dist: &InsertSizeDistribution,
    match_score: i32,
) -> (u8, u8) {
    let insert_size = match fr_insert_size(read1, read2) {
        Some(s) => s,
        None => return (read1.mapq, read2.mapq),
    };

    let paired_score = read1.score + read2.score;
    let o = paired_score + insert_bonus(insert_size, dist, match_score);
    let subo = (paired_score - PEN_UNPAIRED).max(0);

    let mut q_pe = raw_mapq(o - subo, match_score).clamp(0, 60);
    q_pe = (q_pe as f64 * (1.0 - 0.5 * (read1.frac_rep + read2.frac_rep) as f64) + 0.499) as i32;

    let adjust = |q_se: i32| -> u8 {
        let q = if q_se > q_pe {
            q_se
        } else if q_pe < q_se + 40 {
            q_pe
        } else {
            q_se + 40
        };
        q as u8
    };

    (adjust(read1.mapq as i32), adjust(read2.mapq as i32))
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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::types::Cigar;

    #[test]
    fn test_insert_distribution() {
        let sizes: Vec<u32> = vec![200, 198, 202, 199, 201, 200, 197, 203, 200, 201, 199, 202];
        let dist = InsertSizeDistribution::from_insert_sizes(&sizes);

        assert!(
            (dist.mean - 200.0).abs() < 2.0,
            "mean near 200, got {}",
            dist.mean
        );
        assert!(dist.std_dev >= 0.0, "std_dev non-negative");
        assert!(
            dist.upper_bound() < 1000,
            "upper_bound reasonable, got {}",
            dist.upper_bound()
        );
    }

    #[test]
    fn test_pestat_filters_outliers() {
        let mut sizes: Vec<u32> = (0..100).map(|i| 335 + (i % 11)).collect();
        sizes.extend_from_slice(&[50_000_000, 50_000_001, 49_999_999, 50_000_002, 50_000_003]);

        let dist = InsertSizeDistribution::from_insert_sizes(&sizes);

        assert!(
            dist.std_dev < 10_000.0,
            "outliers must be filtered: std_dev={}",
            dist.std_dev
        );
        assert!(
            dist.upper_bound() < 100_000,
            "upper_bound must not be in millions: {}",
            dist.upper_bound()
        );
    }

    #[test]
    fn test_pestat_too_few_pairs() {
        let dist = InsertSizeDistribution::from_insert_sizes(&[300, 320, 310]);
        assert_eq!(dist.mean, 0.0, "too few pairs → mean=0");
        assert_eq!(dist.upper_bound(), 0, "too few pairs → upper_bound=0");
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

    fn scored_read(
        position: usize,
        ref_len: u32,
        reverse: bool,
        score: i32,
        mapq: u8,
        frac_rep: f32,
    ) -> AlignmentResult {
        let mut r = mapped_read(position, ref_len, reverse);
        r.score = score;
        r.mapq = mapq;
        r.frac_rep = frac_rep;
        r
    }

    // --- erfc ---

    #[test]
    fn erfc_matches_known_values() {
        assert!((erfc(0.0) - 1.0).abs() < 1e-6, "erfc(0) ≈ 1");
        assert!((erfc(1.0) - 0.1573).abs() < 1e-3, "erfc(1) ≈ 0.1573");
        assert!((erfc(-1.0) - 1.8427).abs() < 1e-3, "erfc(-1) ≈ 1.8427");
    }

    // --- pair_mapq ---

    #[test]
    fn pair_mapq_unique_proper_pair_stays_60() {
        // Near-mean insert + high scores → q_pe >> 60, so both SE-60 mates stay 60.
        let dist = InsertSizeDistribution::with_params(300.0, 30.0);
        let r1 = scored_read(100, 150, false, 150, 60, 0.0);
        let r2 = scored_read(350, 150, true, 150, 60, 0.0);

        let (m1, m2) = pair_mapq(&r1, &r2, &dist, 1);
        assert_eq!(m1, 60, "high unique pair: read1 mapq stays 60");
        assert_eq!(m2, 60, "high unique pair: read2 mapq stays 60");
    }

    #[test]
    fn pair_mapq_low_se_bumped_toward_pe() {
        // A high-scoring pair lifts a mate whose SE mapq is low via min(q_pe, se+40).
        let dist = InsertSizeDistribution::with_params(300.0, 30.0);
        let r1 = scored_read(100, 150, false, 200, 60, 0.0);
        let r2 = scored_read(350, 150, true, 200, 27, 0.0);

        let (_, m2) = pair_mapq(&r1, &r2, &dist, 1);
        assert!(
            m2 > 27,
            "paired evidence should raise low-SE mate: got {m2}"
        );
    }

    #[test]
    fn pair_mapq_off_mean_insert_lowers_qpe() {
        // Off-mean insert → negative insert bonus → lower q_pe than the near-mean pair.
        let dist = InsertSizeDistribution::with_params(300.0, 30.0);

        let r1_near = scored_read(100, 150, false, 100, 60, 0.0);
        let r2_near = scored_read(250, 150, true, 100, 60, 0.0);
        let r1_far = scored_read(100, 150, false, 100, 60, 0.0);
        let r2_far = scored_read(950, 150, true, 100, 60, 0.0);

        let (m1_near, _) = pair_mapq(&r1_near, &r2_near, &dist, 1);
        let (m1_far, _) = pair_mapq(&r1_far, &r2_far, &dist, 1);

        assert!(
            m1_far <= m1_near,
            "off-mean insert should not raise mapq vs near-mean: far={m1_far}, near={m1_near}"
        );
    }

    #[test]
    fn pair_mapq_frac_rep_downweights() {
        // With a low SE mapq the paired q_pe drives the result, so high frac_rep lowers it.
        let dist = InsertSizeDistribution::with_params(300.0, 30.0);
        let make_pair = |frac: f32| {
            (
                scored_read(100, 150, false, 200, 5, frac),
                scored_read(250, 150, true, 200, 5, frac),
            )
        };

        let (r1_clean, r2_clean) = make_pair(0.0);
        let (r1_rep, r2_rep) = make_pair(0.8);

        let (m1_clean, _) = pair_mapq(&r1_clean, &r2_clean, &dist, 1);
        let (m1_rep, _) = pair_mapq(&r1_rep, &r2_rep, &dist, 1);

        assert!(
            m1_rep < m1_clean,
            "high frac_rep should lower mapq: rep={m1_rep}, clean={m1_clean}"
        );
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
