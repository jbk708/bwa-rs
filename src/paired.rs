//! Paired-end read handling and pairing logic.

use crate::types::{AlignmentResult, Orientation};

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
        (self.mean - 3.0 * self.std_dev).max(0.0) as u32
    }

    pub fn upper_bound(&self) -> u32 {
        (self.mean + 3.0 * self.std_dev) as u32
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
}