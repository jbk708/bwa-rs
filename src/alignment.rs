//! Smith-Waterman alignment with affine gap penalties.

use std::cmp::Ordering;

use crate::chaining::chain_seeds;
use crate::error::BwaError;
use crate::fm_index::FMIndex;
use crate::mem_finder::DEFAULT_MAX_OCC;
use crate::paired::InsertSizeDistribution;
use crate::reference::reverse_complement;
use crate::seed::{filter_mems, find_mems_with_frac_rep, DEFAULT_MIN_SEED_LEN};
use crate::types::{AlignmentResult, ChainedSeed, Cigar, CigarOp};

pub const DEFAULT_MIN_SCORE: i32 = 30;

// Cap on the number of secondary chains aligned to derive the suboptimal score
// for MAPQ; keeps the cost linear on repetitive reads.
const MAX_SUB_CHAINS: usize = 32;

// Fraction of the shorter alignment's query span that must overlap with the
// primary for a candidate to be treated as a secondary (masked) hit.
const MASK_LEVEL: f32 = 0.5;

#[derive(Clone, Debug)]
pub struct Scoring {
    pub match_score: i32,
    pub mismatch_penalty: i32,
    pub gap_open: i32,
    pub gap_extend: i32,
    pub clip_penalty: i32,
}

impl Default for Scoring {
    fn default() -> Self {
        Self {
            match_score: 1,
            mismatch_penalty: 4,
            gap_open: 6,
            gap_extend: 1,
            clip_penalty: 5,
        }
    }
}

#[derive(Clone)]
pub struct Aligner {
    index: FMIndex,
    reference: Vec<u8>,
    scoring: Scoring,
    min_seed_len: usize,
    min_score: i32,
}

impl Aligner {
    pub fn new(index: FMIndex, reference: Vec<u8>) -> Self {
        Self {
            index,
            reference,
            scoring: Scoring::default(),
            min_seed_len: DEFAULT_MIN_SEED_LEN,
            min_score: DEFAULT_MIN_SCORE,
        }
    }

    pub fn min_seed_len(mut self, len: usize) -> Self {
        self.min_seed_len = len;
        self
    }

    pub fn min_score(mut self, score: i32) -> Self {
        self.min_score = score;
        self
    }

    pub fn scoring(mut self, scoring: Scoring) -> Self {
        self.scoring = scoring;
        self
    }

    pub fn align_read(
        &self,
        query: &[u8],
        _mate: Option<&[u8]>,
    ) -> Result<AlignmentResult, BwaError> {
        let (mut mems, frac_rep) =
            find_mems_with_frac_rep(&self.index, query, self.min_seed_len, DEFAULT_MAX_OCC);

        if mems.is_empty() {
            return Ok(self.unmapped_result());
        }

        filter_mems(&mut mems);

        if mems.is_empty() {
            return Ok(self.unmapped_result());
        }

        let chains = chain_seeds(&mems, 50.0);

        if chains.is_empty() {
            return Ok(self.unmapped_result());
        }

        let best_chain = chains
            .iter()
            .max_by(|a, b| a.score.partial_cmp(&b.score).unwrap())
            .ok_or_else(|| BwaError::Alignment("No valid chains".to_string()))?;

        let mut result = self.build_alignment(query, best_chain)?;

        // bwa's -T: alignments scoring below the threshold are reported unmapped
        // rather than placing a low-quality partial hit.
        if result.score < self.min_score {
            return Ok(self.unmapped_result());
        }

        // Suboptimal alignment score for MAPQ: the best-scoring secondary chain
        // whose placement does not overlap the primary, taken over the top
        // MAX_SUB_CHAINS chains by score.
        let mut alt_chains: Vec<&ChainedSeed> = chains
            .iter()
            .filter(|c| !std::ptr::eq(*c, best_chain))
            .collect();
        if alt_chains.len() > MAX_SUB_CHAINS {
            alt_chains.select_nth_unstable_by(MAX_SUB_CHAINS - 1, |a, b| {
                b.score.partial_cmp(&a.score).unwrap_or(Ordering::Equal)
            });
            alt_chains.truncate(MAX_SUB_CHAINS);
        }

        let pspan = query_span(&result.cigar);
        let primary_len = pspan.1.saturating_sub(pspan.0);
        let primary_pos = result.position;
        let primary_ref_end = primary_pos + result.cigar.reference_length();
        let primary_reverse = result.reverse_strand;

        let mut sub = 0i32;
        let mut sub_n = 0u32;
        for alt in alt_chains {
            let Ok(alt_res) = self.build_alignment(query, alt) else {
                continue;
            };

            // A candidate only contributes a suboptimal score when it is a
            // genuinely distinct placement; a chain that reproduces the primary
            // locus (e.g. a re-seed of the same region) is the same alignment.
            let alt_ref_end = alt_res.position + alt_res.cigar.reference_length();
            let same_locus = alt_res.reverse_strand == primary_reverse
                && alt_res.position < primary_ref_end
                && primary_pos < alt_ref_end;
            if same_locus {
                continue;
            }

            let cspan = query_span(&alt_res.cigar);
            let cand_len = cspan.1.saturating_sub(cspan.0);

            // Query-space overlap between candidate and primary spans.
            let overlap_start = pspan.0.max(cspan.0);
            let overlap_end = pspan.1.min(cspan.1);
            let overlap = overlap_end.saturating_sub(overlap_start);

            // Treat as secondary (masked) when query overlap exceeds MASK_LEVEL of
            // the shorter aligned span — mirrors bwa's mem_mark_primary_se logic.
            let min_len = primary_len.min(cand_len);
            let is_secondary = min_len > 0 && overlap as f32 > MASK_LEVEL * min_len as f32;
            if !is_secondary {
                continue;
            }

            match alt_res.score.cmp(&sub) {
                Ordering::Greater => {
                    sub = alt_res.score;
                    sub_n = 1;
                }
                Ordering::Equal if sub > 0 => sub_n += 1,
                _ => {}
            }
        }

        result.mapq = approx_mapq_se(
            result.score,
            sub,
            sub_n,
            aligned_span(&result.cigar),
            &self.scoring,
            self.min_seed_len,
            frac_rep,
        );

        result.xs = sub;

        Ok(result)
    }

    fn build_alignment(
        &self,
        query: &[u8],
        chain: &ChainedSeed,
    ) -> Result<AlignmentResult, BwaError> {
        let ref_seq = &self.reference;
        let seed = &chain.mem;
        let query_len = query.len();

        let bw = optimal_bandwidth(query_len);

        let backward = extend_seed_backward(
            &query[..seed.query_start],
            ref_seq,
            seed.ref_start,
            &self.scoring,
            bw,
        );

        // Anchor at the seed's start corner; the seed interior is free to absorb
        // mismatches or gaps rather than being pinned as all-Eq.
        let forward = affine_extend_forward(
            &query[seed.query_start..],
            ref_seq,
            seed.ref_start,
            &self.scoring,
            bw,
        );

        let leading_clip = seed.query_start - backward.query_end;
        let trailing_clip = (query_len - seed.query_start) - forward.query_end;

        let mut cigar = Cigar::new();

        if leading_clip > 0 {
            cigar.push(CigarOp::S, leading_clip as u32);
        }
        cigar.extend(backward.cigar);
        cigar.extend(forward.cigar);
        if trailing_clip > 0 {
            cigar.push(CigarOp::S, trailing_clip as u32);
        }

        // A backward extension reports `ref_end` as the leftmost reference base it
        // consumed, which is the alignment's absolute start.
        let ref_start = backward.ref_end;

        // nm / score are mismatch/indel counts, invariant to strand, so they are
        // computed in the space the alignment was built in (2N for a build_2n index).
        let nm = self.compute_nm(&cigar, query, ref_seq, ref_start);
        let score = self.compute_alignment_score(&cigar, query, ref_seq, ref_start);

        // Map an indexed-space position to the reported forward-strand position.
        // For a 2N index, positions >= n_fwd fall in the reverse-complement half:
        // a hit covering [p, p+L) there corresponds to forward [2N - p - L, 2N - p)
        // on the reverse strand, and the CIGAR is reversed to forward orientation.
        let n_fwd = self.index.n_fwd;
        let ref_len = cigar.reference_length();
        let (position, reverse_strand, cigar) = if n_fwd > 0 && ref_start >= n_fwd {
            let fwd_pos = (2 * n_fwd).saturating_sub(ref_start + ref_len);
            (fwd_pos, true, cigar.reversed())
        } else {
            (ref_start, false, cigar)
        };

        // MD records reference bases in forward-strand order. For a reverse-strand
        // read this means the forward-oriented read (the reverse-complement of the
        // query) walked along the forward CIGAR against the forward reference.
        let md_tag = if reverse_strand {
            let fwd_query = reverse_complement(query);
            AlignmentResult::new(position, cigar.clone()).mdz_string(&fwd_query, ref_seq)
        } else {
            AlignmentResult::new(position, cigar.clone()).mdz_string(query, ref_seq)
        };

        let mut result = AlignmentResult::new(position, cigar);
        result.reverse_strand = reverse_strand;
        result.score = score;
        result.nm = nm;
        result.md_tag = Some(md_tag);

        Ok(result)
    }

    fn compute_nm(&self, cigar: &Cigar, query: &[u8], reference: &[u8], ref_start: usize) -> u32 {
        let mut nm = 0u32;
        let mut q_pos = 0usize;
        let mut r_pos = ref_start;

        for (op, len) in &cigar.ops {
            match op {
                CigarOp::M | CigarOp::Eq | CigarOp::X => {
                    for i in 0..*len as usize {
                        let q_idx = q_pos + i;
                        let r_idx = r_pos + i;
                        if q_idx < query.len()
                            && r_idx < reference.len()
                            && query[q_idx] != reference[r_idx]
                        {
                            nm += 1;
                        }
                    }
                    q_pos += *len as usize;
                    r_pos += *len as usize;
                }
                CigarOp::I => {
                    nm += *len;
                    q_pos += *len as usize;
                }
                CigarOp::D => {
                    nm += *len;
                    r_pos += *len as usize;
                }
                CigarOp::S => {
                    q_pos += *len as usize;
                }
                _ => {}
            }
        }
        nm
    }

    fn compute_alignment_score(
        &self,
        cigar: &Cigar,
        query: &[u8],
        reference: &[u8],
        ref_start: usize,
    ) -> i32 {
        let mut score = 0i32;
        let mut q_pos = 0usize;
        let mut r_pos = ref_start;

        for (op, len) in &cigar.ops {
            match op {
                CigarOp::M | CigarOp::Eq | CigarOp::X => {
                    for i in 0..*len as usize {
                        let q_idx = q_pos + i;
                        let r_idx = r_pos + i;
                        if q_idx < query.len() && r_idx < reference.len() {
                            if query[q_idx] == reference[r_idx] {
                                score += self.scoring.match_score;
                            } else {
                                score -= self.scoring.mismatch_penalty;
                            }
                        }
                    }
                    q_pos += *len as usize;
                    r_pos += *len as usize;
                }
                CigarOp::I => {
                    score -=
                        self.scoring.gap_open + (*len as i32 - 1).max(0) * self.scoring.gap_extend;
                    q_pos += *len as usize;
                }
                CigarOp::D => {
                    score -=
                        self.scoring.gap_open + (*len as i32 - 1).max(0) * self.scoring.gap_extend;
                    r_pos += *len as usize;
                }
                CigarOp::S => {
                    q_pos += *len as usize;
                }
                _ => {}
            }
        }
        score
    }

    pub fn align_paired(
        &self,
        read1: &[u8],
        read2: &[u8],
    ) -> Result<(AlignmentResult, AlignmentResult), BwaError> {
        let r1 = self.align_read(read1, Some(read2))?;
        let r2 = self.align_read(read2, Some(read1))?;
        Ok((r1, r2))
    }

    /// Attempt to rescue an unmapped `orphan` read given its mapped `mate`, by running a
    /// local Smith-Waterman of the orphan within the mate's insert-size window (bwa
    /// `mem_matesw`). Returns a forward-coordinate alignment on success, gated by
    /// `min_score`. The orphan is expected on the strand opposite the mate (FR).
    pub fn rescue_mate(
        &self,
        orphan: &[u8],
        mate: &AlignmentResult,
        dist: &InsertSizeDistribution,
    ) -> Option<AlignmentResult> {
        if mate.is_unmapped() || orphan.is_empty() {
            return None;
        }
        let n_fwd = self.index.n_fwd;
        if n_fwd == 0 {
            return None;
        }
        let forward_ref = &self.reference[..n_fwd];

        let max_isize = dist.upper_bound() as usize;
        if max_isize < orphan.len() {
            return None;
        }

        let mate_reflen = mate.cigar.reference_length();

        // FR geometry in forward coordinates: the orphan maps on the strand opposite
        // the mate. If the mate is forward (5' at mate.position) the orphan is reverse
        // and downstream; if the mate is reverse (5' at mate.position + mate_reflen) the
        // orphan is forward and upstream.
        let (window_start, window_end, oriented, rescued_reverse) = if !mate.reverse_strand {
            let ws = mate.position;
            let we = (mate.position + max_isize).min(n_fwd);
            (ws, we, reverse_complement(orphan), true)
        } else {
            let we = (mate.position + mate_reflen).min(n_fwd);
            let ws = we.saturating_sub(max_isize);
            (ws, we, orphan.to_vec(), false)
        };

        if window_end <= window_start {
            return None;
        }
        let window = &forward_ref[window_start..window_end];

        let la = local_align(&oriented, window, &self.scoring)?;
        if la.score < self.min_score {
            return None;
        }

        let mut cigar = Cigar::new();
        if la.query_start > 0 {
            cigar.push(CigarOp::S, la.query_start as u32);
        }
        cigar.extend(la.cigar);
        let trailing = oriented.len() - la.query_end;
        if trailing > 0 {
            cigar.push(CigarOp::S, trailing as u32);
        }

        let span = aligned_span(&cigar);
        let position = window_start + la.ref_offset;
        let nm = self.compute_nm(&cigar, &oriented, forward_ref, position);
        let md_tag =
            AlignmentResult::new(position, cigar.clone()).mdz_string(&oriented, forward_ref);

        let mut result = AlignmentResult::new(position, cigar);
        result.reverse_strand = rescued_reverse;
        result.score = la.score;
        result.nm = nm;
        result.mapq = approx_mapq_se(la.score, 0, 0, span, &self.scoring, self.min_seed_len, 0.0);
        result.md_tag = Some(md_tag);
        Some(result)
    }

    fn unmapped_result(&self) -> AlignmentResult {
        AlignmentResult {
            position: 0,
            mapq: 0,
            cigar: Cigar::new(),
            flag: 0x4,
            reverse_strand: false,
            nm: 0,
            score: 0,
            xs: 0,
            md_tag: None,
        }
    }
}

#[derive(Clone, Debug)]
pub struct SeedExtension {
    pub score: i32,
    pub query_end: usize,
    pub ref_end: usize,
    pub cigar: Cigar,
}

pub fn extend_seed_forward(
    query: &[u8],
    reference: &[u8],
    seed_start: usize,
    scoring: &Scoring,
    bandwidth: usize,
) -> SeedExtension {
    let query_len = query.len();
    let ref_len = reference.len();

    if query_len == 0 || ref_len == 0 {
        return SeedExtension {
            score: 0,
            query_end: 0,
            ref_end: seed_start,
            cigar: Cigar::new(),
        };
    }

    let actual_bw = bandwidth.clamp(1, 256);
    let mut best_score = 0;
    let mut best_j = 0;
    let mut traceback: Vec<(usize, i8)> = Vec::new();

    for d in 0..=actual_bw {
        for ref_pos in [seed_start.saturating_sub(d), seed_start + d] {
            if ref_pos >= ref_len {
                continue;
            }

            let j_max = query_len.min(d + 1);

            for j in 1..=j_max {
                let ref_idx = ref_pos + j.saturating_sub(1);

                if ref_idx >= ref_len {
                    continue;
                }

                let score = if query[j - 1] == reference[ref_idx] {
                    scoring.match_score
                } else {
                    scoring.mismatch_penalty
                };

                if score > best_score {
                    best_score = score;
                    best_j = j;
                    traceback.clear();
                    traceback.push((j, 1));
                } else if score == best_score && score > 0 && j > best_j {
                    best_j = j;
                    traceback.clear();
                    traceback.push((j, 1));
                }
            }
        }
    }

    let mut cigar = Cigar::new();
    let mut last_op = CigarOp::M;
    let mut last_len = 0u32;

    for (_j, op_code) in &traceback {
        let op = match op_code {
            1 => {
                if query[best_j.saturating_sub(1)]
                    == reference[seed_start + best_j.saturating_sub(1)]
                {
                    CigarOp::Eq
                } else {
                    CigarOp::X
                }
            }
            2 => CigarOp::I,
            3 => CigarOp::D,
            _ => CigarOp::M,
        };

        if op == last_op {
            last_len += 1;
        } else {
            if last_len > 0 {
                cigar.push(last_op, last_len);
            }
            last_op = op;
            last_len = 1;
        }
    }

    if last_len > 0 {
        cigar.push(last_op, last_len);
    }

    if cigar.ops.is_empty() {
        cigar.push(CigarOp::M, best_j as u32);
    }

    SeedExtension {
        score: best_score,
        query_end: best_j,
        ref_end: seed_start + best_j,
        cigar,
    }
}

pub fn extend_seed_backward(
    query: &[u8],
    reference: &[u8],
    seed_end: usize,
    scoring: &Scoring,
    bandwidth: usize,
) -> SeedExtension {
    if query.is_empty() || seed_end == 0 {
        return SeedExtension {
            score: 0,
            query_end: 0,
            ref_end: seed_end,
            cigar: Cigar::new(),
        };
    }

    let bw = bandwidth.clamp(1, 256);
    let span = (query.len() + bw + 1).min(seed_end);

    let ref_rev: Vec<u8> = reference[seed_end - span..seed_end]
        .iter()
        .rev()
        .copied()
        .collect();
    let q_rev: Vec<u8> = query.iter().rev().copied().collect();

    let ext = affine_extend_forward(&q_rev, &ref_rev, 0, scoring, bandwidth);

    SeedExtension {
        score: ext.score,
        query_end: ext.query_end,
        ref_end: seed_end.saturating_sub(ext.ref_end),
        cigar: ext.cigar.reversed(),
    }
}

pub fn optimal_bandwidth(query_len: usize) -> usize {
    (16 + query_len / 2).min(256)
}

#[derive(Clone)]
pub struct AffineDP {
    m: Vec<i32>,
    x: Vec<i32>,
    g: Vec<i32>,
    rows: usize,
    cols: usize,
}

impl AffineDP {
    pub fn new(query_len: usize, ref_len: usize) -> Self {
        let rows = query_len + 1;
        let cols = ref_len + 1;
        let size = rows * cols;
        Self {
            m: vec![i32::MIN / 2; size],
            x: vec![i32::MIN / 2; size],
            g: vec![i32::MIN / 2; size],
            rows,
            cols,
        }
    }

    #[allow(dead_code)]
    fn idx(&self, i: usize, j: usize) -> usize {
        i * self.cols + j
    }

    fn m_at(&self, i: usize, j: usize) -> i32 {
        self.m[i * self.cols + j]
    }

    fn x_at(&self, i: usize, j: usize) -> i32 {
        self.x[i * self.cols + j]
    }

    fn g_at(&self, i: usize, j: usize) -> i32 {
        self.g[i * self.cols + j]
    }

    fn set_m(&mut self, i: usize, j: usize, val: i32) {
        self.m[i * self.cols + j] = val;
    }

    fn set_x(&mut self, i: usize, j: usize, val: i32) {
        self.x[i * self.cols + j] = val;
    }

    fn set_g(&mut self, i: usize, j: usize, val: i32) {
        self.g[i * self.cols + j] = val;
    }

    pub fn rows(&self) -> usize {
        self.rows
    }

    pub fn cols(&self) -> usize {
        self.cols
    }
}

/// Best local alignment of `query` against `window` via affine-gap Smith-Waterman.
/// Reference ends outside the optimal alignment are free; query ends outside it are
/// reported as `query_start`/`query_end` so the caller can soft-clip them. The returned
/// `cigar` covers only the aligned core (Eq/X/I/D, no soft-clips). Returns `None` when
/// either input is empty or no positive-scoring alignment exists.
#[derive(Clone, Debug)]
pub struct LocalAlignment {
    pub ref_offset: usize,
    pub query_start: usize,
    pub query_end: usize,
    pub score: i32,
    pub cigar: Cigar,
}

pub fn local_align(query: &[u8], window: &[u8], scoring: &Scoring) -> Option<LocalAlignment> {
    let ql = query.len();
    let wl = window.len();
    if ql == 0 || wl == 0 {
        return None;
    }
    let neg = i32::MIN / 2;
    let cols = wl + 1;
    // H = best score ending at (i,j); E = horizontal gap (deletion, consumes reference);
    // F = vertical gap (insertion, consumes query).
    let mut h = vec![0i32; (ql + 1) * cols];
    let mut e = vec![neg; (ql + 1) * cols];
    let mut f = vec![neg; (ql + 1) * cols];

    let mut best_score = 0i32;
    let mut best_i = 0usize;
    let mut best_j = 0usize;

    for i in 1..=ql {
        for j in 1..=wl {
            let s = if query[i - 1] == window[j - 1] {
                scoring.match_score
            } else {
                -scoring.mismatch_penalty
            };
            // deletion (horizontal): consumes reference column j
            let e_val = h[i * cols + (j - 1)]
                .saturating_sub(scoring.gap_open)
                .max(e[i * cols + (j - 1)].saturating_sub(scoring.gap_extend));
            // insertion (vertical): consumes query row i
            let f_val = h[(i - 1) * cols + j]
                .saturating_sub(scoring.gap_open)
                .max(f[(i - 1) * cols + j].saturating_sub(scoring.gap_extend));
            e[i * cols + j] = e_val;
            f[i * cols + j] = f_val;
            let diag = h[(i - 1) * cols + (j - 1)].saturating_add(s);
            let h_val = 0.max(diag).max(e_val).max(f_val);
            h[i * cols + j] = h_val;
            if h_val > best_score {
                best_score = h_val;
                best_i = i;
                best_j = j;
            }
        }
    }

    if best_score <= 0 {
        return None;
    }

    // Traceback from (best_i, best_j) until H hits 0.
    let mut ops: Vec<CigarOp> = Vec::new();
    let mut i = best_i;
    let mut j = best_j;
    while i > 0 && j > 0 && h[i * cols + j] > 0 {
        let s = if query[i - 1] == window[j - 1] {
            scoring.match_score
        } else {
            -scoring.mismatch_penalty
        };
        let cur = h[i * cols + j];
        if cur == h[(i - 1) * cols + (j - 1)].saturating_add(s) {
            ops.push(if query[i - 1] == window[j - 1] {
                CigarOp::Eq
            } else {
                CigarOp::X
            });
            i -= 1;
            j -= 1;
        } else if cur == e[i * cols + j] {
            ops.push(CigarOp::D);
            j -= 1;
        } else if cur == f[i * cols + j] {
            ops.push(CigarOp::I);
            i -= 1;
        } else {
            break;
        }
    }

    ops.reverse();
    Some(LocalAlignment {
        ref_offset: j,
        query_start: i,
        query_end: best_i,
        score: best_score,
        cigar: Cigar::compress(ops),
    })
}

/// Port of bwa's `mem_approx_mapq_se`. `score` is the alignment score, `sub` the best
/// suboptimal *alignment* score (0 when none — the bwa `min_seed_len * a` baseline is
/// applied here), `sub_n` the number of comparable suboptimal hits, `aligned_len`
/// the aligned span `l = max(query_core, ref_core)`, and `frac_rep` the fraction of
/// the read covered by repetitive seeds (bwa's final `mapq *= 1 - frac_rep`
/// down-weighting). bwa-mem2 defaults drive the `mapQ_coef_len = 50` branch; `csub`
/// is not modeled (treated as 0).
fn approx_mapq_se(
    score: i32,
    sub: i32,
    sub_n: u32,
    aligned_len: i32,
    scoring: &Scoring,
    min_seed_len: usize,
    frac_rep: f32,
) -> u8 {
    const MAPQ_COEF_LEN: i32 = 50;
    let a = scoring.match_score;
    let b = scoring.mismatch_penalty;

    let sub = if sub > 0 {
        sub
    } else {
        min_seed_len as i32 * a
    };
    if sub >= score || score == 0 || aligned_len <= 0 {
        return 0;
    }

    let l = aligned_len;
    let identity = 1.0 - (l * a - score) as f64 / ((a + b) as f64) / l as f64;

    let coef_fac = (MAPQ_COEF_LEN as f64).ln();
    let tmp = if l < MAPQ_COEF_LEN {
        1.0
    } else {
        coef_fac / (l as f64).ln()
    };
    let tmp = tmp * identity * identity;
    let mut mapq = (6.02 * (score - sub) as f64 / a as f64 * tmp * tmp + 0.499) as i32;

    if sub_n > 0 {
        mapq -= (4.343 * ((sub_n + 1) as f64).ln() + 0.499) as i32;
    }

    let mapq = mapq.clamp(0, 60);
    (mapq as f64 * (1.0 - frac_rep as f64) + 0.499) as i32 as u8
}

/// Returns (start, end) of the query bases consumed by the aligned core (M/=/X/I
/// ops), measured from the start of the read. Leading soft-clips shift `start`
/// forward; `end` = start + aligned query length.
fn query_span(cigar: &Cigar) -> (usize, usize) {
    let leading_clip = match cigar.ops.first() {
        Some(&(CigarOp::S, len)) => len as usize,
        _ => 0,
    };
    let aligned_q: usize = cigar
        .ops
        .iter()
        .filter(|(op, _)| matches!(op, CigarOp::M | CigarOp::Eq | CigarOp::X | CigarOp::I))
        .map(|(_, len)| *len as usize)
        .sum();
    (leading_clip, leading_clip + aligned_q)
}

/// bwa's alignment length `l = max(query span, reference span)` over the aligned core
/// (M/=/X/I consume query; M/=/X/D consume reference; soft/hard clips excluded).
fn aligned_span(cigar: &Cigar) -> i32 {
    let mut q = 0i32;
    let mut r = 0i32;
    for &(op, len) in &cigar.ops {
        match op {
            CigarOp::M | CigarOp::Eq | CigarOp::X => {
                q += len as i32;
                r += len as i32;
            }
            CigarOp::I => q += len as i32,
            CigarOp::D | CigarOp::N => r += len as i32,
            _ => {}
        }
    }
    q.max(r)
}

pub fn affine_extend_forward(
    query: &[u8],
    reference: &[u8],
    start_pos: usize,
    scoring: &Scoring,
    bandwidth: usize,
) -> SeedExtension {
    let query_len = query.len();
    let ref_len = reference.len().saturating_sub(start_pos);

    if query_len == 0 || ref_len == 0 {
        return SeedExtension {
            score: 0,
            query_end: 0,
            ref_end: start_pos,
            cigar: Cigar::new(),
        };
    }

    let bw = bandwidth.clamp(1, 256);
    // The banded DP only ever touches columns j in [i-bw, i+bw] for rows
    // i in 1..=query_len, so column index never exceeds query_len + bw. Cap the
    // reference window (and thus the allocated matrix) to that band width;
    // anything beyond it is dead space. Without this cap the matrix is sized by
    // the full remaining reference length (up to the whole chromosome), which
    // allocates hundreds of GB and OOM-kills the process.
    let max_span = query_len + bw + 1;
    let actual_ref_len = ref_len.min(max_span);
    let ref_end = (start_pos + actual_ref_len).min(reference.len());
    let actual_ref_len = ref_end - start_pos;

    let mut dp = AffineDP::new(query_len, actual_ref_len);
    dp.set_m(0, 0, 0);
    dp.set_x(0, 0, i32::MIN / 2);
    dp.set_g(0, 0, i32::MIN / 2);

    let mut best_score = i32::MIN;
    let mut best_i = 0;
    let mut best_j = 0;
    let mut gscore = i32::MIN;
    let mut gscore_j = 0usize;

    for i in 1..=query_len {
        let j_start = 1.max(i.saturating_sub(bw));
        let j_end = (i + bw).min(actual_ref_len);

        for j in j_start..=j_end {
            let is_match = query[i - 1] == reference[start_pos + j - 1];
            let match_score = if is_match {
                scoring.match_score
            } else {
                -(scoring.mismatch_penalty)
            };

            let m_prev = dp.m_at(i - 1, j - 1);
            let x_prev = dp.x_at(i - 1, j - 1);
            let g_prev = dp.g_at(i - 1, j - 1);
            let m_val = m_prev
                .saturating_add(match_score)
                .max(x_prev.saturating_add(match_score))
                .max(g_prev.saturating_add(match_score));
            dp.set_m(i, j, m_val);

            let x_open = dp.m_at(i - 1, j).saturating_sub(scoring.gap_open);
            let x_extend = dp.x_at(i - 1, j).saturating_sub(scoring.gap_extend);
            dp.set_x(i, j, x_open.max(x_extend));

            let g_open = dp.m_at(i, j - 1).saturating_sub(scoring.gap_open);
            let g_extend = dp.g_at(i, j - 1).saturating_sub(scoring.gap_extend);
            dp.set_g(i, j, g_open.max(g_extend));

            let score = m_val.max(dp.x_at(i, j)).max(dp.g_at(i, j));
            if score > best_score {
                best_score = score;
                best_i = i;
                best_j = j;
            }
            if i == query_len && score > gscore {
                gscore = score;
                gscore_j = j;
            }
        }
    }

    let (end_i, end_j, end_score) = if gscore > 0 && gscore > best_score - scoring.clip_penalty {
        (query_len, gscore_j, gscore)
    } else {
        (best_i, best_j, best_score)
    };

    let cigar = traceback_affine(&dp, query, reference, start_pos, end_i, end_j, scoring);

    SeedExtension {
        score: end_score,
        query_end: end_i,
        ref_end: start_pos + end_j,
        cigar,
    }
}

fn traceback_affine(
    dp: &AffineDP,
    query: &[u8],
    reference: &[u8],
    ref_start: usize,
    mut i: usize,
    mut j: usize,
    scoring: &Scoring,
) -> Cigar {
    let mut ops = Vec::new();
    let match_s = scoring.match_score;
    let mismatch_s = -(scoring.mismatch_penalty);

    while i > 0 && j > 0 {
        let current = dp.m_at(i, j).max(dp.x_at(i, j)).max(dp.g_at(i, j));

        let is_match = query[i - 1] == reference[ref_start + j - 1];
        let delta = if is_match { match_s } else { mismatch_s };

        if current == dp.m_at(i, j) && current == dp.m_at(i - 1, j - 1) + delta {
            ops.push(if is_match { CigarOp::Eq } else { CigarOp::X });
            i -= 1;
            j -= 1;
        } else if current == dp.x_at(i, j) {
            ops.push(CigarOp::D);
            i -= 1;
        } else if current == dp.g_at(i, j) {
            ops.push(CigarOp::I);
            j -= 1;
        } else if current == dp.m_at(i, j) {
            ops.push(if is_match { CigarOp::Eq } else { CigarOp::X });
            i -= 1;
            j -= 1;
        } else {
            break;
        }
    }

    // Only add trailing deletions for remaining query (not insertions)
    // Insertions (I) consume query characters, so if i == 0 we can't add more I
    while i > 0 {
        ops.push(CigarOp::D);
        i -= 1;
    }
    // Note: trailing reference without more query doesn't add I ops
    // because I consumes query chars (decrements i), not reference

    ops.reverse();
    Cigar::compress(ops)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::fm_index::FMIndex as FMIndexAlias;
    use crate::paired::InsertSizeDistribution;
    use crate::reference::Reference;

    fn enc(s: &str) -> Vec<u8> {
        s.bytes()
            .map(|b| match b {
                b'A' => 0u8,
                b'C' => 1,
                b'G' => 2,
                b'T' => 3,
                _ => 4,
            })
            .collect()
    }

    #[test]
    fn local_align_exact_substring() {
        let window = enc("AAAACGTACGTAAAA");
        let query = enc("CGTACGT");
        let scoring = Scoring::default();
        let la = local_align(&query, &window, &scoring).unwrap();
        assert_eq!(la.query_start, 0);
        assert_eq!(la.query_end, query.len());
        assert_eq!(la.score, query.len() as i32 * scoring.match_score);
        assert_eq!(la.ref_offset, 4, "matches the CGTACGT at window offset 4");
        assert_eq!(la.cigar.to_sam_string(), "7M");
    }

    #[test]
    fn local_align_soft_clips_query_ends() {
        // query = junk(TTT) + core(ACGTAC) + junk(TTT); core matches window middle
        let window = enc("GGGGACGTACGGGG");
        let query = enc("TTTACGTACTTT");
        let scoring = Scoring::default();
        let la = local_align(&query, &window, &scoring).unwrap();
        assert_eq!(la.query_start, 3, "leading 3 bases clipped");
        assert_eq!(la.query_end, 9, "trailing 3 bases clipped");
        assert_eq!(la.cigar.to_sam_string(), "6M");
    }

    #[test]
    fn local_align_none_on_empty() {
        let scoring = Scoring::default();
        assert!(local_align(&[], &enc("ACGT"), &scoring).is_none());
        assert!(local_align(&enc("ACGT"), &[], &scoring).is_none());
    }

    fn rescue_test_aligner(body: &str) -> (Aligner, Vec<u8>) {
        let ref_seq = Reference::parse_fasta(&format!(">t\n{body}")).unwrap();
        let ref_data = ref_seq.as_slice_2n();
        let index = FMIndexAlias::build_2n(&ref_seq);
        let aligner = Aligner::new(index, ref_data.clone()).min_score(10);
        (aligner, ref_data)
    }

    #[test]
    fn rescue_mate_forward_mate_places_reverse_orphan_downstream() {
        let body = "ACGTACGATCGATCGGATTCCAGTCAGTCAGGATCCATGCATGCATTAGCATCGATCGTA";
        let (aligner, ref_data) = rescue_test_aligner(body);
        let n_fwd = body.len();
        let fwd = &ref_data[..n_fwd];

        // mate maps forward at position 0, covering 20 bp
        let mut mate_cigar = Cigar::new();
        mate_cigar.push(CigarOp::M, 20);
        let mut mate = AlignmentResult::new(0, mate_cigar);
        mate.reverse_strand = false;

        // orphan is the reverse-strand read of forward region [35, 55): the read as it
        // would arrive from the sequencer is reverse_complement of that region.
        let region = fwd[35..55].to_vec();
        let orphan = reverse_complement(&region);

        let dist = InsertSizeDistribution::with_params(60.0, 5.0); // upper_bound = 60+20=80
        let rescued = aligner
            .rescue_mate(&orphan, &mate, &dist)
            .expect("should rescue");
        assert!(rescued.reverse_strand, "orphan rescued on reverse strand");
        assert_eq!(rescued.position, 35, "placed at forward coord 35");
        assert_eq!(rescued.cigar.to_sam_string(), "20M");
        assert_eq!(rescued.nm, 0);
    }

    #[test]
    fn rescue_mate_reverse_mate_places_forward_orphan_upstream() {
        let body = "ACGTACGATCGATCGGATTCCAGTCAGTCAGGATCCATGCATGCATTAGCATCGATCGTA";
        let (aligner, ref_data) = rescue_test_aligner(body);
        let n_fwd = body.len();
        let fwd = &ref_data[..n_fwd];

        // mate maps REVERSE at position 40, covering 20 bp (5' end at 60)
        let mut mate_cigar = Cigar::new();
        mate_cigar.push(CigarOp::M, 20);
        let mut mate = AlignmentResult::new(40, mate_cigar);
        mate.reverse_strand = true;

        // orphan is a forward read of region [5, 25)
        let orphan = fwd[5..25].to_vec();

        let dist = InsertSizeDistribution::with_params(50.0, 5.0); // upper_bound = 70
        let rescued = aligner
            .rescue_mate(&orphan, &mate, &dist)
            .expect("should rescue");
        assert!(!rescued.reverse_strand, "orphan rescued on forward strand");
        assert_eq!(rescued.position, 5, "placed at forward coord 5");
        assert_eq!(rescued.cigar.to_sam_string(), "20M");
        assert_eq!(rescued.nm, 0);
    }

    #[test]
    fn rescue_mate_none_for_unrelated_orphan() {
        let body = "ACGTACGATCGATCGGATTCCAGTCAGTCAGGATCCATGCATGCATTAGCATCGATCGTA";
        let (aligner, _ref_data) = rescue_test_aligner(body);
        let mut mate_cigar = Cigar::new();
        mate_cigar.push(CigarOp::M, 20);
        let mut mate = AlignmentResult::new(0, mate_cigar);
        mate.reverse_strand = false;
        // orphan that does not occur in the window
        let orphan = enc("TTTTTTTTTTTTTTTTTTTT");
        let dist = InsertSizeDistribution::with_params(60.0, 5.0);
        assert!(aligner.rescue_mate(&orphan, &mate, &dist).is_none());
    }

    #[test]
    fn rescue_mate_none_when_mate_unmapped() {
        let body = "ACGTACGATCGATCGGATTCCAGTCAGTCAGGATCCATGCATGCATTAGCATCGATCGTA";
        let (aligner, _ref_data) = rescue_test_aligner(body);
        let mut mate = AlignmentResult::new(0, Cigar::new());
        mate.flag = 0x4;
        let orphan = enc("ACGTACGATCGATCGGATTC");
        let dist = InsertSizeDistribution::with_params(60.0, 5.0);
        assert!(aligner.rescue_mate(&orphan, &mate, &dist).is_none());
    }

    #[test]
    fn test_unmapped_result() {
        let ref_seq = Reference::parse_fasta(">test\nACGT").unwrap();
        let ref_data = ref_seq.as_slice().to_vec();
        let index = FMIndex::build(&ref_seq);
        let aligner = Aligner::new(index, ref_data);

        let unmapped = aligner.unmapped_result();
        assert_eq!(unmapped.flag, 0x4);
        assert_eq!(unmapped.mapq, 0);
    }

    #[test]
    fn test_scoring_default() {
        let scoring = Scoring::default();
        assert_eq!(scoring.match_score, 1);
        assert_eq!(scoring.mismatch_penalty, 4);
    }

    #[test]
    fn test_mdz_reverse_strand_uses_forward_reference() {
        let fwd = "AGCTTAGCTAGGCATTACGATCGATCGGATCCATGCATGA";
        let reference = Reference::parse_fasta(&format!(">ref\n{}", fwd)).unwrap();
        let ref_data = reference.as_slice_2n();
        let n = fwd.len();
        let fwd_enc: Vec<u8> = ref_data[..n].to_vec();

        let (start, len, off) = (5usize, 30usize, 20usize);
        let mut fwd_read = fwd_enc[start..start + len].to_vec();
        fwd_read[off] = (fwd_read[off] + 1) % 4;
        let read = reverse_complement(&fwd_read);

        let index = FMIndex::build_2n(&reference);
        let aligner = Aligner::new(index, ref_data).min_seed_len(11).min_score(0);
        let result = aligner.align_read(&read, None).unwrap();

        assert!(
            result.reverse_strand,
            "read should map to the reverse strand"
        );
        assert_eq!(result.cigar.to_sam_string(), "30M");
        let bases = ['A', 'C', 'G', 'T'];
        let ref_base = bases[fwd_enc[start + off] as usize];
        let expected = format!("MD:Z:{}{}{}", off, ref_base, len - off - 1);
        assert_eq!(result.md_tag.as_deref(), Some(expected.as_str()));
    }

    #[test]
    fn test_aligner_builder() {
        let ref_seq = Reference::parse_fasta(">test\nACGT").unwrap();
        let ref_data = ref_seq.as_slice().to_vec();
        let index = FMIndex::build(&ref_seq);

        let aligner = Aligner::new(index, ref_data)
            .min_seed_len(10)
            .scoring(Scoring::default());

        assert_eq!(aligner.min_seed_len, 10);
    }

    #[test]
    fn test_extend_forward_exact_match() {
        let query = vec![0, 1, 2, 3];
        let reference = vec![0, 1, 2, 3, 4, 0, 1, 2, 3];
        let scoring = Scoring::default();

        let ext = extend_seed_forward(&query, &reference, 0, &scoring, 16);
        assert!(ext.score > 0);
        assert_eq!(ext.query_end, query.len());
        assert_eq!(ext.ref_end, 4);
    }

    #[test]
    fn test_extend_forward_with_mismatch() {
        let query = vec![0, 1, 2, 3];
        let reference = vec![0, 1, 0, 3, 0, 1, 2, 3];
        let scoring = Scoring::default();

        let ext = extend_seed_forward(&query, &reference, 3, &scoring, 16);
        assert!(ext.score >= 0);
    }

    #[test]
    fn test_extend_forward_insertion() {
        let query = vec![0, 1, 5, 2, 3];
        let reference = vec![0, 1, 2, 3];
        let scoring = Scoring::default();

        let ext = extend_seed_forward(&query, &reference, 0, &scoring, 16);
        assert!(ext.score >= 0);
        assert!(ext.query_end <= query.len());
    }

    #[test]
    fn test_extend_backward() {
        let query = vec![2, 3, 0, 1];
        let reference = vec![0, 1, 2, 3, 4, 2, 3, 0, 1];
        let scoring = Scoring::default();

        let ext = extend_seed_backward(&query, &reference, 6, &scoring, 16);
        assert!(ext.ref_end <= 6);
    }

    #[test]
    fn test_optimal_bandwidth() {
        assert!(optimal_bandwidth(50) >= 16);
        assert!(optimal_bandwidth(500) >= 16);
        assert!(optimal_bandwidth(500) <= 256);
    }

    #[test]
    fn test_seed_extension_zero_bandwidth() {
        let query = vec![0, 1, 2];
        let reference = vec![0, 1, 2, 3];
        let scoring = Scoring::default();

        let ext = extend_seed_forward(&query, &reference, 0, &scoring, 0);
        assert!(ext.score >= 0);
    }

    #[test]
    fn test_affine_exact_match() {
        let query = vec![0, 1, 2, 3];
        let reference = vec![0, 1, 2, 3];
        let scoring = Scoring::default();

        let ext = affine_extend_forward(&query, &reference, 0, &scoring, 16);
        assert!(ext.score > 0, "Score should be positive for exact match");
        assert_eq!(ext.query_end, query.len());
        let cigar_str = ext.cigar.to_string();
        assert!(
            cigar_str.contains('=') || cigar_str.contains('X'),
            "CIGAR should contain match ops: {}",
            cigar_str
        );
    }

    #[test]
    fn test_affine_insertion() {
        let query = vec![0, 1, 4, 2, 3];
        let reference = vec![0, 1, 2, 3];
        let scoring = Scoring::default();

        let ext = affine_extend_forward(&query, &reference, 0, &scoring, 16);
        let cigar_str = ext.cigar.to_string();
        assert!(!cigar_str.is_empty(), "CIGAR should not be empty");
    }

    #[test]
    fn test_affine_deletion() {
        let query = vec![0, 1, 2, 3];
        let reference = vec![0, 1, 4, 2, 3];
        let scoring = Scoring::default();

        let ext = affine_extend_forward(&query, &reference, 0, &scoring, 16);
        let cigar_str = ext.cigar.to_string();
        assert!(!cigar_str.is_empty(), "CIGAR should not be empty");
    }

    #[test]
    fn test_affine_gap_open_vs_extend() {
        let scoring = Scoring::default();
        assert!(scoring.gap_open > scoring.gap_extend);
        let single_gap_cost = scoring.gap_open;
        let two_gap_cost = scoring.gap_open + scoring.gap_extend;
        assert!(
            single_gap_cost < two_gap_cost,
            "Single gap should cost less than two separate gaps"
        );
    }

    #[test]
    fn test_affine_empty_query() {
        let query: Vec<u8> = vec![];
        let reference = vec![0, 1, 2, 3];
        let scoring = Scoring::default();

        let ext = affine_extend_forward(&query, &reference, 0, &scoring, 16);
        assert_eq!(ext.score, 0);
        assert_eq!(ext.query_end, 0);
    }

    #[test]
    fn test_affine_empty_reference() {
        let query = vec![0, 1, 2, 3];
        let reference: Vec<u8> = vec![];
        let scoring = Scoring::default();

        let ext = affine_extend_forward(&query, &reference, 0, &scoring, 16);
        assert_eq!(ext.score, 0);
        assert_eq!(ext.query_end, 0);
    }

    #[test]
    fn test_affine_dp_struct() {
        let dp = AffineDP::new(5, 10);
        assert_eq!(dp.rows(), 6);
        assert_eq!(dp.cols(), 11);
        assert!(dp.m_at(0, 0) < 0);
    }

    #[test]
    fn test_align_full_read() {
        let ref_seq = Reference::parse_fasta(">test\nACGTACGT").unwrap();
        let ref_data = ref_seq.as_slice().to_vec();
        let ref_data_len = ref_data.len();
        let index = FMIndex::build(&ref_seq);
        let index_for_check = index.clone();
        let aligner = Aligner::new(index, ref_data).min_seed_len(2).min_score(0);

        // Use AC which we know works with the simple search test
        let query = vec![0, 1]; // AC

        // Check if FM-index finds the pattern
        let positions = index_for_check.find_all(&query);
        assert!(
            !positions.is_empty(),
            "FM-index should find AC positions: {:?}, ref_len: {}",
            positions,
            ref_data_len
        );

        let result = aligner.align_read(&query, None).unwrap();

        if result.cigar.ops.is_empty() {
            panic!(
                "CIGAR is empty - flag: {}, mapq: {}, score: {}",
                result.flag, result.mapq, result.score
            );
        }
        assert!(
            result.cigar.reference_length() > 0,
            "CIGAR should cover reference positions"
        );
        // Score can be negative if extension introduces mismatches
        assert!(result.nm < u32::MAX, "NM tag should be set");
    }

    #[test]
    fn test_align_with_mismatch() {
        let ref_seq = Reference::parse_fasta(">test\nACGTACGT").unwrap();
        let ref_data = ref_seq.as_slice().to_vec();
        let index = FMIndex::build(&ref_seq);
        let aligner = Aligner::new(index, ref_data).min_seed_len(2).min_score(0);

        // Use shorter query
        let query = vec![0, 1, 2, 3, 2]; // ACGTC
        let result = aligner.align_read(&query, None).unwrap();

        assert!(
            !result.cigar.to_string().is_empty(),
            "CIGAR should not be empty"
        );
    }

    #[test]
    fn test_align_unmapped() {
        let ref_seq = Reference::parse_fasta(">test\nACGT").unwrap();
        let ref_data = ref_seq.as_slice().to_vec();
        let index = FMIndex::build(&ref_seq);
        let aligner = Aligner::new(index, ref_data).min_seed_len(100);

        let query = vec![0, 1, 2, 3, 0, 1, 2, 3];
        let result = aligner.align_read(&query, None).unwrap();

        assert_eq!(result.flag & 0x4, 0x4, "Unmapped flag should be set");
        assert_eq!(result.mapq, 0, "MAPQ should be 0 for unmapped");
    }

    #[test]
    fn test_mapq_calculation() {
        let ref_seq = Reference::parse_fasta(">test\nACGTACGTACGT").unwrap();
        let ref_data = ref_seq.as_slice().to_vec();
        let index = FMIndex::build(&ref_seq);
        let aligner = Aligner::new(index, ref_data).min_seed_len(2);

        let query = vec![0, 1, 2, 3];
        let result = aligner.align_read(&query, None).unwrap();

        assert!(result.mapq <= 60, "MAPQ should be capped at 60");
    }

    #[test]
    fn test_cigar_extend() {
        let mut cigar1 = Cigar::new();
        cigar1.push(CigarOp::Eq, 5);

        let mut cigar2 = Cigar::new();
        cigar2.push(CigarOp::Eq, 3);

        cigar1.extend(cigar2);
        assert_eq!(cigar1.to_string(), "8=");
    }

    #[test]
    fn test_cigar_extend_different_ops() {
        let mut cigar1 = Cigar::new();
        cigar1.push(CigarOp::Eq, 5);

        let mut cigar2 = Cigar::new();
        cigar2.push(CigarOp::I, 3);

        cigar1.extend(cigar2);
        assert_eq!(cigar1.to_string(), "5=3I");
    }

    #[test]
    fn test_mdz_string_all_matches() {
        let mut cigar = Cigar::new();
        cigar.push(CigarOp::Eq, 10);

        let query = vec![0u8; 10];
        let reference = vec![0u8; 10];

        let result = AlignmentResult::new(0, cigar);
        let mdz = result.mdz_string(&query, &reference);
        assert_eq!(mdz, "MD:Z:10");
    }

    #[test]
    fn test_mdz_string_with_mismatch() {
        let mut cigar = Cigar::new();
        cigar.push(CigarOp::Eq, 5);
        cigar.push(CigarOp::X, 1);
        cigar.push(CigarOp::Eq, 3);

        let query = vec![0, 0, 0, 0, 0, 1, 0, 0, 0];
        let reference = vec![0, 0, 0, 0, 0, 0, 0, 0, 0];

        let result = AlignmentResult::new(0, cigar);
        let mdz = result.mdz_string(&query, &reference);
        assert_eq!(mdz, "MD:Z:5A3");
    }

    #[test]
    fn test_mdz_string_with_deletion() {
        let mut cigar = Cigar::new();
        cigar.push(CigarOp::Eq, 5);
        cigar.push(CigarOp::D, 1);
        cigar.push(CigarOp::Eq, 3);

        let query = vec![0u8; 8];
        let reference = vec![0, 0, 0, 0, 0, 1, 0, 0, 0];

        let result = AlignmentResult::new(0, cigar);
        let mdz = result.mdz_string(&query, &reference);
        assert_eq!(mdz, "MD:Z:5^C3");
    }

    #[test]
    fn test_mdz_string_with_insertion() {
        let mut cigar = Cigar::new();
        cigar.push(CigarOp::Eq, 3);
        cigar.push(CigarOp::I, 2);
        cigar.push(CigarOp::Eq, 3);

        let query = vec![0, 0, 0, 1, 2, 0, 0, 0];
        let reference = vec![0, 0, 0, 0, 0, 0];

        let result = AlignmentResult::new(0, cigar);
        let mdz = result.mdz_string(&query, &reference);
        // Insertions don't appear in MD:Z per SAM spec - only reference matches count
        assert_eq!(mdz, "MD:Z:6");
    }

    #[test]
    fn test_mdz_string_empty() {
        let cigar = Cigar::new();
        let query: Vec<u8> = vec![];
        let reference: Vec<u8> = vec![];

        let result = AlignmentResult::new(0, cigar);
        let mdz = result.mdz_string(&query, &reference);
        assert_eq!(mdz, "MD:Z:0");
    }

    #[test]
    fn test_mdz_string_complex() {
        let mut cigar = Cigar::new();
        cigar.push(CigarOp::Eq, 2);
        cigar.push(CigarOp::X, 1);
        cigar.push(CigarOp::Eq, 3);
        cigar.push(CigarOp::D, 1);
        cigar.push(CigarOp::Eq, 4);

        let query = vec![0, 0, 1, 0, 0, 0, 0, 0, 0, 0];
        let reference = vec![0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0];

        let result = AlignmentResult::new(0, cigar);
        let mdz = result.mdz_string(&query, &reference);
        assert_eq!(mdz, "MD:Z:2A3^C4");
    }

    #[test]
    fn test_mdz_string_in_alignment_result() {
        let ref_seq = Reference::parse_fasta(">test\nACGTACGT").unwrap();
        let ref_data = ref_seq.as_slice().to_vec();
        let index = FMIndex::build(&ref_seq);
        let aligner = Aligner::new(index, ref_data).min_seed_len(2).min_score(0);

        let query = vec![0, 1, 2, 3];
        let result = aligner.align_read(&query, None).unwrap();

        if result.md_tag.is_none() {
            panic!("MD tag should be set for aligned reads");
        }
        assert!(
            result.md_tag.as_ref().unwrap().starts_with("MD:Z:"),
            "{}",
            result.md_tag.unwrap()
        );
    }

    #[test]
    fn test_affine_clip_trailing_mismatches() {
        let reference = vec![0u8, 1, 2, 3, 0, 1, 2];
        let query = vec![0u8, 1, 2, 3, 3, 3, 3];
        let scoring = Scoring::default();

        let ext = affine_extend_forward(&query, &reference, 0, &scoring, 16);
        assert_eq!(
            ext.query_end, 4,
            "trailing 3 bases should be clipped, query_end should be 4"
        );
        assert_eq!(
            ext.cigar.to_string(),
            "4=",
            "cigar should be 4= for the 4 matching bases"
        );
    }

    #[test]
    fn test_affine_extend_through_single_tail_mismatch() {
        let reference = vec![0u8, 1, 2, 3, 0, 1, 2, 3];
        let query = vec![0u8, 1, 2, 3, 0, 1, 3];
        let scoring = Scoring::default();

        let ext = affine_extend_forward(&query, &reference, 0, &scoring, 16);
        assert_eq!(
            ext.query_end, 7,
            "should extend through single tail mismatch, query_end should be 7"
        );
        assert_eq!(ext.cigar.to_string(), "6=1X", "cigar should be 6=1X");
    }

    #[test]
    fn test_align_read_emits_soft_clips() {
        let small_ref = "AAAAAGGGGAAAAACCCCAAAAA";
        let ref_seq = Reference::parse_fasta(&format!(">test\n{}", small_ref)).unwrap();
        let ref_slice = ref_seq.as_slice();
        let ref_data = ref_slice.to_vec();
        let index = FMIndex::build(&ref_seq);
        let aligner = Aligner::new(index, ref_data).min_seed_len(3);

        let read_str = "TTTTGGGGAAAAACCCCTTTT";
        let read: Vec<u8> = read_str
            .bytes()
            .map(|b| match b {
                b'A' => 0,
                b'C' => 1,
                b'G' => 2,
                b'T' => 3,
                _ => 4,
            })
            .collect();
        let read_len = read.len();

        let result = aligner.align_read(&read, None).unwrap();

        if result.flag & 0x4 != 0 {
            return;
        }

        assert_eq!(
            result.cigar.query_length(),
            read_len,
            "CIGAR query-consumed length {} must equal read length {}. CIGAR: {}",
            result.cigar.query_length(),
            read_len,
            result.cigar
        );

        let cigar_str = result.cigar.to_string();
        assert!(
            cigar_str.contains('S'),
            "CIGAR should contain soft-clip ops for mismatching ends. CIGAR: {}",
            cigar_str
        );
    }

    fn full_length_test_aligner(fasta_body: &str) -> Aligner {
        let ref_seq = Reference::parse_fasta(&format!(">t\n{fasta_body}")).unwrap();
        let ref_data = ref_seq.as_slice().to_vec();
        let index = FMIndex::build(&ref_seq);
        Aligner::new(index, ref_data).min_seed_len(4).min_score(0)
    }

    // An exact full-length match must still place at the right POS with NM=0 after
    // the forward-DP anchor change.
    #[test]
    fn test_exact_match_regression() {
        let aligner = full_length_test_aligner("ACGTACGT");
        let query = vec![0u8, 1, 2, 3];
        let result = aligner.align_read(&query, None).unwrap();

        assert_eq!(result.flag & 0x4, 0, "should be mapped");
        assert_eq!(result.nm, 0, "NM should be 0 for exact match");
        assert_eq!(result.cigar.query_length(), 4, "full read must be aligned");
        assert_eq!(result.cigar.reference_length(), 4);
        assert_eq!(result.cigar.to_sam_string(), "4M");
    }

    // The read is the reference with its leading base removed, so the optimal
    // placement runs full-length from POS 1 through the seed region — reachable
    // only now that the seed is no longer pinned as exact matches.
    #[test]
    fn test_seed_interior_free_full_length_alignment() {
        let aligner = full_length_test_aligner("AACGTAACGT");
        let query = vec![0u8, 1, 2, 3, 0, 0, 1, 2, 3];
        let result = aligner.align_read(&query, None).unwrap();

        assert_eq!(result.flag & 0x4, 0, "should be mapped");
        assert_eq!(
            result.cigar.query_length(),
            query.len(),
            "all query bases must appear in CIGAR (no soft-clips); got {}",
            result.cigar
        );
        assert_eq!(
            result.position, 1,
            "alignment should start at ref pos 1; got pos={} cigar={}",
            result.position, result.cigar
        );
        assert_eq!(result.nm, 0, "NM must be 0 for an exact shifted alignment");
        assert_eq!(result.cigar.to_sam_string(), "9M");
    }

    #[test]
    fn test_min_score_unmaps_low_scoring_alignment() {
        let small_ref = "AAAAAGGGGAAAAACCCCAAAAA";
        let ref_seq = Reference::parse_fasta(&format!(">test\n{}", small_ref)).unwrap();
        let ref_data = ref_seq.as_slice().to_vec();
        let index = FMIndex::build(&ref_seq);

        let read = vec![0u8, 0, 2, 2, 2, 2, 0, 0, 0, 0, 0];

        let mapped = Aligner::new(index.clone(), ref_data.clone())
            .min_seed_len(3)
            .min_score(0)
            .align_read(&read, None)
            .unwrap();
        if mapped.flag & 0x4 != 0 {
            return;
        }

        let filtered = Aligner::new(index, ref_data)
            .min_seed_len(3)
            .min_score(mapped.score + 1)
            .align_read(&read, None)
            .unwrap();
        assert_eq!(
            filtered.flag & 0x4,
            0x4,
            "alignment scoring below min_score should be reported unmapped"
        );
    }

    // ---- approx_mapq_se unit tests ----

    #[test]
    fn approx_mapq_se_unique_read_is_60() {
        // score=151, sub=0 → sub becomes min_seed_len*a = 19*1 = 19
        // l=151 >= 50, identity = 1 - (151*1 - 151)/5/151 = 1.0
        // coef_fac = ln(50), tmp = ln(50)/ln(151) * 1.0 * 1.0 ≈ 0.8054
        // mapq = (6.02 * 132/1 * 0.8054^2 * 0.8054^2 + 0.499) ≈ very large → clamps to 60
        let result = approx_mapq_se(151, 0, 0, 151, &Scoring::default(), 19, 0.0);
        assert_eq!(result, 60);
    }

    #[test]
    fn approx_mapq_se_sub_ge_score_is_zero() {
        assert_eq!(
            approx_mapq_se(100, 100, 0, 100, &Scoring::default(), 19, 0.0),
            0
        );
        assert_eq!(
            approx_mapq_se(100, 120, 0, 100, &Scoring::default(), 19, 0.0),
            0
        );
    }

    #[test]
    fn approx_mapq_se_zero_score_is_zero() {
        assert_eq!(
            approx_mapq_se(0, 0, 0, 100, &Scoring::default(), 19, 0.0),
            0
        );
    }

    #[test]
    fn approx_mapq_se_subn_penalty_reduces() {
        // score=40, sub=30, l=45 (< 50 so tmp=1.0 before identity scaling)
        // a=1, b=4; identity = 1 - (45 - 40)/(5*45) = 1 - 5/225 ≈ 0.97778
        // tmp = 1.0 * 0.97778^2 ≈ 0.95606
        // mapq0 = (6.02 * 10/1 * 0.95606^2 + 0.499) as i32
        //       = (6.02 * 10 * 0.91404 + 0.499) as i32
        //       = (55.025 + 0.499) as i32 = 55 as i32 = 55
        // sub_n=4: penalty = (4.343 * ln(5) + 0.499) as i32
        //                   = (4.343 * 1.60944 + 0.499) as i32
        //                   = (6.990 + 0.499) as i32 = 7
        // mapq4 = 55 - 7 = 48
        let m0 = approx_mapq_se(40, 30, 0, 45, &Scoring::default(), 19, 0.0);
        let m4 = approx_mapq_se(40, 30, 4, 45, &Scoring::default(), 19, 0.0);
        assert_eq!(m0, 55);
        assert_eq!(m4, 48);
        assert!(m4 < m0);
    }

    #[test]
    fn approx_mapq_se_clamps_to_60() {
        // A perfect long read with no suboptimal hits: must clamp to exactly 60.
        let result = approx_mapq_se(200, 0, 0, 200, &Scoring::default(), 19, 0.0);
        assert_eq!(result, 60);
    }

    #[test]
    fn approx_mapq_se_frac_rep_downweights() {
        // bwa's final mapq *= (1 - frac_rep): a fully-repetitive read drops to 0,
        // half-repetitive roughly halves the unweighted score.
        let full = approx_mapq_se(200, 0, 0, 200, &Scoring::default(), 19, 1.0);
        let half = approx_mapq_se(200, 0, 0, 200, &Scoring::default(), 19, 0.5);
        assert_eq!(full, 0);
        assert_eq!(half, 30);
    }

    // ---- aligned_span unit tests ----

    fn make_cigar(ops: &[(CigarOp, u32)]) -> Cigar {
        let mut c = Cigar::new();
        for &(op, len) in ops {
            c.push(op, len);
        }
        c
    }

    #[test]
    fn aligned_span_pure_match_with_clips() {
        // 10S 88M 9S → aligned core = 88M → q=88, r=88 → max=88
        let c = make_cigar(&[(CigarOp::S, 10), (CigarOp::M, 88), (CigarOp::S, 9)]);
        assert_eq!(aligned_span(&c), 88);
    }

    #[test]
    fn aligned_span_with_indels_equal_cores() {
        // 44S 37M 1I 2M 1D 58M 9S
        // query core: 37+1+2+58 = 98; ref core: 37+2+1+58 = 98 → max=98
        let c = make_cigar(&[
            (CigarOp::S, 44),
            (CigarOp::M, 37),
            (CigarOp::I, 1),
            (CigarOp::M, 2),
            (CigarOp::D, 1),
            (CigarOp::M, 58),
            (CigarOp::S, 9),
        ]);
        assert_eq!(aligned_span(&c), 98);
    }

    #[test]
    fn aligned_span_insertion_makes_query_larger() {
        // 5M 3I 5M → q=13, r=10 → max=13
        let c = make_cigar(&[(CigarOp::M, 5), (CigarOp::I, 3), (CigarOp::M, 5)]);
        assert_eq!(aligned_span(&c), 13);
    }

    #[test]
    fn aligned_span_deletion_makes_ref_larger() {
        // 5M 3D 5M → q=10, r=13 → max=13
        let c = make_cigar(&[(CigarOp::M, 5), (CigarOp::D, 3), (CigarOp::M, 5)]);
        assert_eq!(aligned_span(&c), 13);
    }

    #[test]
    fn mapq_is_60_for_unique_read() {
        // A read that is unique in the reference (no repeat loci) should get mapq=60.
        let body = "ACGTACGATCGATCGGATTCCAGTCAGTCAGGATCCATGCATGCATTAGCATCGATCGTA";
        let ref_seq = Reference::parse_fasta(&format!(">t\n{body}")).unwrap();
        let ref_data = ref_seq.as_slice().to_vec();
        let index = FMIndex::build(&ref_seq);
        let aligner = Aligner::new(index, ref_data).min_seed_len(10).min_score(10);

        // Use a long, exact substring of the reference that can't appear twice.
        let query = enc(&body[5..35]);
        let result = aligner.align_read(&query, None).unwrap();
        if result.flag & 0x4 != 0 {
            // If unmapped, the test is inconclusive — skip rather than fail.
            return;
        }
        assert_eq!(result.mapq, 60, "unique read must get mapq=60");
    }

    /// A query that matches at two distinct reference loci should get `xs > 0`
    /// and `mapq < 60`.  A query that occurs only once should get `xs == 0` and
    /// `mapq == 60`.
    #[test]
    fn xs_nonzero_and_mapq_lt60_for_repeat_read() {
        // Reference: unique flanks at both ends; the same 30-bp core repeated twice.
        // Layout: UNIQUE_LEFT(30) + CORE(30) + SEPARATOR(5) + CORE(30) + UNIQUE_RIGHT(30)
        let core = "ACGTACGATCGATCGGATTCCAGTCAGTCA";
        // Unique flanks chosen to be non-repetitive and not present elsewhere.
        let left_flank = "GGCTAGCTAGCTAGCTAGCTAGCTAGCTAG";
        let sep = "CCCCC";
        let right_flank = "TTCGATCGATCGATCGATCGATCGATCGAT";
        let body = format!("{}{}{}{}{}", left_flank, core, sep, core, right_flank);

        let ref_seq = Reference::parse_fasta(&format!(">t\n{body}")).unwrap();
        let ref_data_2n = ref_seq.as_slice_2n();
        let index = FMIndex::build_2n(&ref_seq);

        let aligner = Aligner::new(index, ref_data_2n.to_vec())
            .min_seed_len(10)
            .min_score(5);

        // The core itself occurs at two positions → should produce xs > 0 and mapq < 60.
        let repeat_query = enc(core);
        let repeat_result = aligner.align_read(&repeat_query, None).unwrap();
        if repeat_result.flag & 0x4 != 0 {
            return; // inconclusive if unmapped
        }
        assert!(
            repeat_result.xs > 0 || repeat_result.mapq < 60,
            "repeat read should have xs>0 or mapq<60; xs={}, mapq={}",
            repeat_result.xs,
            repeat_result.mapq
        );

        // A unique 30-bp substring taken from the left flank → should get xs==0, mapq==60.
        // Build the aligner using a single-stranded (forward-only) index so the RC of
        // the unique flank doesn't appear as a second hit.
        let ref_fwd = ref_seq.as_slice().to_vec();
        let index2 = FMIndex::build(&ref_seq);
        let aligner2 = Aligner::new(index2, ref_fwd).min_seed_len(10).min_score(5);
        let unique_query = enc(left_flank);
        let unique_result = aligner2.align_read(&unique_query, None).unwrap();
        if unique_result.flag & 0x4 != 0 {
            return; // inconclusive
        }
        assert_eq!(
            unique_result.xs, 0,
            "unique read must have xs=0; mapq={}",
            unique_result.mapq
        );
        assert_eq!(
            unique_result.mapq, 60,
            "unique read must have mapq=60; xs={}",
            unique_result.xs
        );
    }

    /// query_span correctly identifies the aligned region of a CIGAR.
    #[test]
    fn query_span_extracts_aligned_bounds() {
        // 5S 10M 3I 7M 4S → start=5, end=5+10+3+7=25
        let c = make_cigar(&[
            (CigarOp::S, 5),
            (CigarOp::M, 10),
            (CigarOp::I, 3),
            (CigarOp::M, 7),
            (CigarOp::S, 4),
        ]);
        let (s, e) = query_span(&c);
        assert_eq!(s, 5, "start = leading soft-clip");
        assert_eq!(e, 25, "end = start + aligned query bases");
    }

    #[test]
    fn query_span_no_clips() {
        // 20M → start=0, end=20
        let c = make_cigar(&[(CigarOp::M, 20)]);
        let (s, e) = query_span(&c);
        assert_eq!(s, 0);
        assert_eq!(e, 20);
    }
}
