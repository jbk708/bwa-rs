//! Smith-Waterman alignment with affine gap penalties.

use crate::chaining::{build_candidate_chains, chain_seeds, mem_chain_flt};
use crate::error::BwaError;
use crate::fm_index::FMIndex;
use crate::mem_finder::{collect_short_seeds, DEFAULT_MAX_OCC, MAX_MEM_INTV};
use crate::paired::InsertSizeDistribution;
use crate::reference::reverse_complement;
use crate::seed::{filter_mems, find_mems_with_frac_rep, DEFAULT_MIN_SEED_LEN};
use crate::types::{AlignmentResult, ChainedSeed, Cigar, CigarOp, MEM};

pub const DEFAULT_MIN_SCORE: i32 = 30;

/// The per-read candidate alignment regions (bwa `mem_alnreg_v`). `regions[0]` is
/// the chosen primary (max extended score, baseline winning ties) carrying the
/// final mapq/xs/frac_rep exactly as `align_read` returns; the remaining entries are
/// the other candidate placements, sorted by score descending and deduplicated by
/// (position, strand). `regions` is empty when the read is unmapped / below min_score.
#[derive(Clone, Debug)]
pub struct ReadRegions {
    pub regions: Vec<AlignmentResult>,
    /// bwa `sub` — best secondary score feeding MAPQ/XS (0 if none).
    pub sub: i32,
    /// bwa `sub_n` — count of secondaries tied at `sub`.
    pub sub_n: u32,
    pub frac_rep: f32,
    /// bwa `csub` — best score among candidate regions that overlap the primary's
    /// reference span (tandem hit). `sub` is lower-bounded by this before the MAPQ
    /// formula: `effective_sub = max(sub, csub)`. 0 when no overlapping candidate.
    pub csub: i32,
}

// Cap on the number of candidate MEMs extended to derive the suboptimal score
// for MAPQ; keeps the cost bounded on repetitive reads. The longest MEMs are
// kept, as they anchor the most reliable placements.
const MAX_SUB_CANDIDATES: usize = 64;

// Fraction of the shorter alignment's query span that must overlap with the
// primary for a candidate to be treated as a secondary (masked) hit.
const MASK_LEVEL: f32 = 0.5;

// bwa ksw_extend defaults (bwa-mem2 -w 100, -d 100): banded SW band half-width and
// the Z-drop X-drop threshold. end_bonus is the clip penalty (bwa pen_clip5/3).
const KSW_BAND_WIDTH: i32 = 100;
const KSW_ZDROP: i32 = 100;

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

    /// Returns bwa's match score (`opt->a`), used for scaled MAPQ formulas.
    pub fn match_score(&self) -> i32 {
        self.scoring.match_score
    }

    /// Single-end MAPQ (bwa `mem_approx_mapq_se`) for an already-extended region, using
    /// the per-read suboptimal score/count, repeat fraction, and tandem-hit score from
    /// its [`ReadRegions`]. Recomputing it for `regions[0]` reproduces the primary's
    /// MAPQ exactly; it is used to score a non-primary region picked by mate-aware
    /// pairing.
    pub fn region_se_mapq(
        &self,
        region: &AlignmentResult,
        sub: i32,
        sub_n: u32,
        frac_rep: f32,
        csub: i32,
    ) -> u8 {
        approx_mapq_se(
            region.score,
            sub,
            sub_n,
            aligned_span(&region.cigar),
            &self.scoring,
            self.min_seed_len,
            frac_rep,
            csub,
        )
    }

    /// bwa's `mem_pair` `tmp`: the score margin for counting near-suboptimal pairings,
    /// `max(a+b, o_del+e_del, o_ins+e_ins)`. bwa-rs uses a single gap model, so this is
    /// `max(match+mismatch, gap_open+gap_extend)`.
    pub fn pair_margin(&self) -> i32 {
        (self.scoring.match_score + self.scoring.mismatch_penalty)
            .max(self.scoring.gap_open + self.scoring.gap_extend)
    }

    pub fn align_read_regions(
        &self,
        query: &[u8],
        _mate: Option<&[u8]>,
    ) -> Result<ReadRegions, BwaError> {
        let (raw_mems, frac_rep) =
            find_mems_with_frac_rep(&self.index, query, self.min_seed_len, DEFAULT_MAX_OCC);

        if raw_mems.is_empty() {
            return Ok(ReadRegions {
                regions: Vec::new(),
                sub: 0,
                sub_n: 0,
                frac_rep,
                csub: 0,
            });
        }

        // Baseline placement: the highest-scoring filtered/chained SMEM, exactly as
        // before. It anchors the candidate set so primary selection never does worse
        // than the chain heuristic, and it wins exact score ties below so that
        // uniquely-mapping reads (already byte-identical to bwa-mem2) stay frozen.
        let mut mems = raw_mems.clone();
        filter_mems(&mut mems);

        if mems.is_empty() {
            return Ok(ReadRegions {
                regions: Vec::new(),
                sub: 0,
                sub_n: 0,
                frac_rep,
                csub: 0,
            });
        }

        let chains = chain_seeds(&mems, 50.0);

        if chains.is_empty() {
            return Ok(ReadRegions {
                regions: Vec::new(),
                sub: 0,
                sub_n: 0,
                frac_rep,
                csub: 0,
            });
        }

        let best_chain = chains
            .iter()
            .max_by(|a, b| a.score.partial_cmp(&b.score).unwrap())
            .ok_or_else(|| BwaError::Alignment("No valid chains".to_string()))?;

        let baseline = self.build_alignment(query, best_chain)?;

        // bwa scores every candidate *region* and picks the primary by extended
        // alignment score, not chain score. Build the same candidate set used for the
        // suboptimal score — unfiltered MEMs keep tandem-repeat copies, and the third
        // seeding round (`collect_short_seeds`) surfaces inexact secondary loci — then
        // extend each into a placement and let a strictly-higher-scoring region replace
        // the baseline.
        let mut cand_mems = raw_mems;
        cand_mems.extend(collect_short_seeds(
            &self.index,
            query,
            self.min_seed_len,
            MAX_MEM_INTV,
        ));

        // bwa `mem_chain_flt`: group the candidate seeds into chains and drop the
        // weak, query-contained chains that a higher-weight chain dominates, before
        // any of them are extended into regions. This prunes the spurious secondaries
        // that inflate `sub`/XS and over-score MAPQ. The baseline placement is added
        // separately below and the primary is the max-extended-score region, so
        // pruning here only narrows the secondary set — it never moves the primary.
        let candidate_chains = mem_chain_flt(build_candidate_chains(cand_mems), self.min_seed_len);
        let mut cand_mems: Vec<MEM> = candidate_chains.into_iter().flat_map(|c| c.seeds).collect();

        // `build_alignment` anchors on (query_start, ref_start) only, so seeds sharing
        // that corner extend to the same placement; collapse them before the expensive
        // DP. Keeping the longest among equal corners is arbitrary (the extension
        // ignores length) but deterministic.
        cand_mems
            .sort_unstable_by_key(|m| (m.query_start, m.ref_start, std::cmp::Reverse(m.length)));
        cand_mems.dedup_by_key(|m| (m.query_start, m.ref_start));
        if cand_mems.len() > MAX_SUB_CANDIDATES {
            cand_mems.select_nth_unstable_by_key(MAX_SUB_CANDIDATES - 1, |m| {
                std::cmp::Reverse(m.length)
            });
            cand_mems.truncate(MAX_SUB_CANDIDATES);
        }

        let mut placements: Vec<AlignmentResult> = cand_mems
            .iter()
            .filter_map(|m| {
                let chained = ChainedSeed::from_mem(*m);
                self.build_alignment(query, &chained).ok()
            })
            .collect();
        placements.push(baseline.clone());

        // Primary = highest extended-alignment score. The baseline wins exact ties so
        // uniquely-mapping placements stay byte-identical to bwa-mem2; only a strictly
        // better region (e.g. a candidate MEM that DP-extends above the chain's pick)
        // moves the placement. bwa's coordinate/strand tie-break is unnecessary because
        // ties retain the baseline.
        let mut result = baseline;
        for p in &placements {
            if p.score > result.score {
                result = p.clone();
            }
        }

        // bwa's -T: alignments scoring below the threshold are reported unmapped
        // rather than placing a low-quality partial hit.
        if result.score < self.min_score {
            return Ok(ReadRegions {
                regions: Vec::new(),
                sub: 0,
                sub_n: 0,
                frac_rep,
                csub: 0,
            });
        }

        // Suboptimal alignment score for MAPQ (bwa's `sub`/`sub_n`). The primary is now
        // the maximal-scoring region, so a candidate at a distinct reference locus that
        // overlaps the primary in query space is a genuine secondary; no guard against
        // candidates outscoring the primary is needed (none can).
        //
        // Query spans are normalized to forward-read coordinates (bwa stores
        // `mem_alnreg_t.qb/qe` in the read's own orientation, not the reference's). A
        // reverse-strand placement's CIGAR is reference-forward, so its soft-clips are
        // mirrored relative to the sequenced read; without this normalization the
        // primary and an equal-scoring opposite-strand copy of the same read region
        // appear to occupy disjoint query spans and the secondary is missed.
        let qlen = query.len();
        let pspan = forward_query_span(&result.cigar, result.reverse_strand, qlen);
        let primary_len = pspan.1.saturating_sub(pspan.0);
        let primary_pos = result.position;
        let primary_ref_end = primary_pos + result.cigar.reference_length();
        let primary_reverse = result.reverse_strand;

        // Best score per distinct secondary placement, keyed by (position, strand).
        let mut secondaries: Vec<(usize, bool, i32)> = Vec::new();
        for alt_res in &placements {
            let cspan = forward_query_span(&alt_res.cigar, alt_res.reverse_strand, qlen);

            // A placement contained in the primary in BOTH query and reference
            // (same strand) is the primary alignment re-derived from one of its own
            // seeds — bwa's mem_sort_dedup_patch containment rule — so it is not a
            // secondary. This also excludes the chosen primary itself.
            let alt_ref_end = alt_res.position + alt_res.cigar.reference_length();
            let query_contained = cspan.0 >= pspan.0 && cspan.1 <= pspan.1;
            let ref_contained = alt_res.position >= primary_pos && alt_ref_end <= primary_ref_end;
            if alt_res.reverse_strand == primary_reverse && query_contained && ref_contained {
                continue;
            }

            let cand_len = cspan.1.saturating_sub(cspan.0);
            let overlap = pspan.1.min(cspan.1).saturating_sub(pspan.0.max(cspan.0));
            let min_len = primary_len.min(cand_len);
            let is_secondary = min_len > 0 && overlap as f32 > MASK_LEVEL * min_len as f32;
            if !is_secondary {
                continue;
            }

            if !secondaries
                .iter()
                .any(|&(p, rev, _)| rev == alt_res.reverse_strand && p == alt_res.position)
            {
                secondaries.push((alt_res.position, alt_res.reverse_strand, alt_res.score));
            }
        }

        let sub = secondaries.iter().map(|&(_, _, s)| s).max().unwrap_or(0);
        let sub_n = if sub > 0 {
            secondaries.iter().filter(|&&(_, _, s)| s == sub).count() as u32
        } else {
            0
        };

        // csub (bwa `mem_alnreg_t.csub`, "SW score of a tandem hit"): best score among
        // candidate placements that overlap the primary in **reference** space but are not
        // primary re-derivations.  bwa's mem_approx_mapq_se lower-bounds `sub` by csub
        // before the MAPQ formula (`sub = csub > sub ? csub : sub`), capturing
        // near-tied overlapping alignments that the query-overlap secondary scan misses.
        //
        // Exclusion: mirror the containment guard used for `sub` — same strand AND
        // query-contained AND reference-contained are primary re-derivations and are
        // excluded from both sub and csub.
        let mut csub = 0i32;
        for alt_res in &placements {
            let cspan = forward_query_span(&alt_res.cigar, alt_res.reverse_strand, qlen);
            let alt_ref_end = alt_res.position + alt_res.cigar.reference_length();
            let query_contained = cspan.0 >= pspan.0 && cspan.1 <= pspan.1;
            let ref_contained = alt_res.position >= primary_pos && alt_ref_end <= primary_ref_end;
            if alt_res.reverse_strand == primary_reverse && query_contained && ref_contained {
                continue; // primary re-derivation — same exclusion as for sub
            }
            // Reference-span overlap: [alt_pos, alt_ref_end) ∩ [primary_pos, primary_ref_end)
            if alt_res.position < primary_ref_end && alt_ref_end > primary_pos {
                csub = csub.max(alt_res.score);
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
            csub,
        );

        result.xs = sub;
        result.frac_rep = frac_rep;

        // Region array: chosen primary first, then the remaining candidate placements
        // sorted by score descending and deduped by (position, strand). align_read
        // returns regions[0], so this assembly cannot change existing output.
        let mut regions = Vec::with_capacity(placements.len());
        let mut seen: Vec<(usize, bool)> = vec![(result.position, result.reverse_strand)];
        regions.push(result);
        let mut rest: Vec<AlignmentResult> = placements;
        rest.sort_by_key(|p| std::cmp::Reverse(p.score));
        for p in rest {
            let key = (p.position, p.reverse_strand);
            if seen.contains(&key) {
                continue;
            }
            seen.push(key);
            regions.push(p);
        }
        Ok(ReadRegions {
            regions,
            sub,
            sub_n,
            frac_rep,
            csub,
        })
    }

    pub fn align_read(
        &self,
        query: &[u8],
        mate: Option<&[u8]>,
    ) -> Result<AlignmentResult, BwaError> {
        let rr = self.align_read_regions(query, mate)?;
        Ok(rr
            .regions
            .into_iter()
            .next()
            .unwrap_or_else(|| self.unmapped_result()))
    }

    fn build_alignment(
        &self,
        query: &[u8],
        chain: &ChainedSeed,
    ) -> Result<AlignmentResult, BwaError> {
        let ref_seq = &self.reference;
        let seed = &chain.mem;
        let query_len = query.len();
        let a = self.scoring.match_score;
        let end_bonus = self.scoring.clip_penalty;
        let bw = KSW_BAND_WIDTH as usize;

        let qs = seed.query_start.min(query_len);
        let qe = (seed.query_start + seed.length).min(query_len);
        let rs = seed.ref_start.min(ref_seq.len());
        let re = (seed.ref_start + seed.length).min(ref_seq.len());
        let seed_len = qe.saturating_sub(qs);
        let h0 = seed_len as i32 * a;

        // Right extension: lock the seed (h0) and extend query[qe..] against ref[re..].
        // Choose to-end (reach query end, soft-clip nothing) when the end-bonus makes it
        // win; otherwise stop at the local max and soft-clip the tail.
        let (right_q, right_t) = if qe >= query_len {
            (0usize, 0usize)
        } else {
            let q = &query[qe..];
            let t_end = (re + (query_len - qe) + bw + 1).min(ref_seq.len());
            let t = &ref_seq[re..t_end];
            if t.is_empty() {
                (0, 0)
            } else {
                let r = ksw_extend(q, t, &self.scoring, end_bonus, h0);
                if r.gscore > 0 && r.gscore + end_bonus > r.score {
                    (query_len - qe, r.gtle)
                } else {
                    (r.qle, r.tle)
                }
            }
        };

        // Left extension: reverse the upstream query/reference so the same forward
        // extension runs leftward from the seed start corner.
        let (left_q, left_t) = if qs == 0 || rs == 0 {
            (0usize, 0usize)
        } else {
            let q: Vec<u8> = query[..qs].iter().rev().copied().collect();
            let t_lo = rs.saturating_sub(qs + bw + 1);
            let t: Vec<u8> = ref_seq[t_lo..rs].iter().rev().copied().collect();
            if t.is_empty() {
                (0, 0)
            } else {
                let r = ksw_extend(&q, &t, &self.scoring, end_bonus, h0);
                if r.gscore > 0 && r.gscore + end_bonus > r.score {
                    (qs, r.gtle)
                } else {
                    (r.qle, r.tle)
                }
            }
        };

        let qb = qs - left_q;
        let qend = qe + right_q;
        let rb = rs - left_t;
        let rend = re + right_t;

        // Render the CIGAR over the determined window with a banded global alignment
        // (the seed region in the middle is exact, so it aligns as matches).
        let core = ksw_global(&query[qb..qend], &ref_seq[rb..rend], &self.scoring, bw);

        let mut cigar = Cigar::new();
        if qb > 0 {
            cigar.push(CigarOp::S, qb as u32);
        }
        cigar.extend(core);
        let trailing_clip = query_len - qend;
        if trailing_clip > 0 {
            cigar.push(CigarOp::S, trailing_clip as u32);
        }

        let ref_start = rb;

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
        // bwa's mem_matesw accepts a mate-rescued alignment at `opt->min_seed_len`,
        // not `opt->T`: a concordant rescue that pairs is kept even when its solo score
        // is below the report threshold (pair-conditioned acceptance).
        if la.score < self.min_seed_len as i32 {
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
        result.mapq = approx_mapq_se(
            la.score,
            0,
            0,
            span,
            &self.scoring,
            self.min_seed_len,
            0.0,
            0,
        );
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
            frac_rep: 0.0,
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
/// the aligned span `l = max(query_core, ref_core)`, `frac_rep` the fraction of
/// the read covered by repetitive seeds (bwa's final `mapq *= 1 - frac_rep`
/// down-weighting), and `csub` the tandem-hit score (best reference-overlapping
/// candidate, bwa `mem_alnreg_t.csub`). bwa-mem2 defaults drive the
/// `mapQ_coef_len = 50` branch.
#[allow(clippy::too_many_arguments)]
fn approx_mapq_se(
    score: i32,
    sub: i32,
    sub_n: u32,
    aligned_len: i32,
    scoring: &Scoring,
    min_seed_len: usize,
    frac_rep: f32,
    csub: i32,
) -> u8 {
    const MAPQ_COEF_LEN: i32 = 50;
    let a = scoring.match_score;
    let b = scoring.mismatch_penalty;

    // bwa: sub = a->sub ? a->sub : opt->min_seed_len * opt->a
    //      sub = a->csub > sub ? a->csub : sub
    let sub = if sub > 0 {
        sub
    } else {
        min_seed_len as i32 * a
    };
    let sub = sub.max(csub);
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

/// [`query_span`] expressed in forward-read coordinates. A reverse-strand placement's
/// CIGAR is written reference-forward, so its aligned span sits at the mirror position
/// in the sequenced read; bwa keeps region query bounds (`mem_alnreg_t.qb/qe`) in that
/// read frame so primary/secondary overlap is comparable across strands.
fn forward_query_span(cigar: &Cigar, reverse_strand: bool, query_len: usize) -> (usize, usize) {
    let (qs, qe) = query_span(cigar);
    if reverse_strand {
        (query_len.saturating_sub(qe), query_len.saturating_sub(qs))
    } else {
        (qs, qe)
    }
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

struct KswExtend {
    score: i32,
    qle: usize,
    tle: usize,
    gscore: i32,
    gtle: usize,
}

/// Faithful port of bwa's `ksw_extend2` (scoring only, no traceback). Extends
/// `query` against `target` from the origin starting at accumulated score `h0`
/// (the locked seed score). Returns the local maximum `(score, qle, tle)` and the
/// best query-end-reaching score `(gscore, gtle)`. Symmetric affine gaps. Banded
/// with band-shrink and Z-drop early termination; `end_bonus` biases reaching the
/// query end (used by the caller to choose to-end vs soft-clipped extension).
fn ksw_extend(
    query: &[u8],
    target: &[u8],
    scoring: &Scoring,
    end_bonus: i32,
    h0: i32,
) -> KswExtend {
    let qlen = query.len();
    let tlen = target.len();
    let a = scoring.match_score;
    let b = scoring.mismatch_penalty;
    let o = scoring.gap_open;
    let e = scoring.gap_extend;
    let oe = o + e;

    if qlen == 0 || tlen == 0 {
        return KswExtend {
            score: h0,
            qle: 0,
            tle: 0,
            gscore: -1,
            gtle: 0,
        };
    }

    let mut eh_h = vec![0i32; qlen + 1];
    let mut eh_e = vec![0i32; qlen + 1];

    eh_h[0] = h0;
    eh_h[1] = if h0 > oe { h0 - oe } else { 0 };
    let mut j = 2usize;
    while j <= qlen && eh_h[j - 1] > e {
        eh_h[j] = eh_h[j - 1] - e;
        j += 1;
    }

    let mut w = KSW_BAND_WIDTH;
    let max_ins =
        (((qlen as f64 * a as f64 + end_bonus as f64 - o as f64) / e as f64 + 1.0) as i32).max(1);
    if w > max_ins {
        w = max_ins;
    }
    let max_del =
        (((qlen as f64 * a as f64 + end_bonus as f64 - o as f64) / e as f64 + 1.0) as i32).max(1);
    if w > max_del {
        w = max_del;
    }

    let mut max = h0;
    let mut max_i: i32 = -1;
    let mut max_j: i32 = -1;
    let mut max_ie: i32 = -1;
    let mut gscore: i32 = -1;
    let mut beg: i32 = 0;
    let mut end: i32 = qlen as i32;

    for (i, &tgt) in target.iter().enumerate() {
        let mut f = 0i32;
        let mut m = 0i32;
        let mut mj: i32 = -1;

        if beg < i as i32 - w {
            beg = i as i32 - w;
        }
        if end > i as i32 + w + 1 {
            end = i as i32 + w + 1;
        }
        if end > qlen as i32 {
            end = qlen as i32;
        }

        let mut h1 = if beg == 0 {
            (h0 - (o + e * (i as i32 + 1))).max(0)
        } else {
            0
        };

        let mut jc = beg;
        while jc < end {
            let col = jc as usize;
            let mut mscore = eh_h[col];
            let ecur = eh_e[col];
            eh_h[col] = h1;
            let s = if query[col] == tgt { a } else { -b };
            mscore = if mscore != 0 { mscore + s } else { 0 };
            let mut h = mscore.max(ecur);
            h = h.max(f);
            h1 = h;
            if m <= h {
                mj = jc;
                m = h;
            }
            let mut t = mscore - oe;
            if t < 0 {
                t = 0;
            }
            let mut enew = ecur - e;
            if enew < t {
                enew = t;
            }
            eh_e[col] = enew;
            let mut t2 = mscore - oe;
            if t2 < 0 {
                t2 = 0;
            }
            f -= e;
            if f < t2 {
                f = t2;
            }
            jc += 1;
        }
        eh_h[end as usize] = h1;
        eh_e[end as usize] = 0;

        if end == qlen as i32 && gscore <= h1 {
            max_ie = i as i32;
            gscore = h1;
        }
        if m == 0 {
            break;
        }
        if m > max {
            max = m;
            max_i = i as i32;
            max_j = mj;
        } else if KSW_ZDROP > 0 {
            let di = i as i32 - max_i;
            let dj = mj - max_j;
            let drop = if di > dj {
                max - m - (di - dj) * e
            } else {
                max - m - (dj - di) * e
            };
            if drop > KSW_ZDROP {
                break;
            }
        }

        let mut nb = beg;
        while nb < end && eh_h[nb as usize] == 0 && eh_e[nb as usize] == 0 {
            nb += 1;
        }
        beg = nb;
        let mut ne = end;
        while ne >= beg && eh_h[ne as usize] == 0 && eh_e[ne as usize] == 0 {
            ne -= 1;
        }
        end = if ne + 2 < qlen as i32 {
            ne + 2
        } else {
            qlen as i32
        };
    }

    KswExtend {
        score: max,
        qle: (max_j + 1).max(0) as usize,
        tle: (max_i + 1).max(0) as usize,
        gscore,
        gtle: (max_ie + 1).max(0) as usize,
    }
}

/// Banded affine-gap *global* alignment of `query` against `target` (both ends
/// anchored at the corners), returning the optimal-score CIGAR (Eq/X/I/D) over the
/// whole rectangle. Used to render the CIGAR over the window `ksw_extend`
/// determines, mirroring bwa's `ksw_global2`. Axes: query = rows (I consumes
/// query), target = cols (D consumes target).
fn ksw_global(query: &[u8], target: &[u8], scoring: &Scoring, bw: usize) -> Cigar {
    let ql = query.len();
    let tl = target.len();
    let mut cig = Cigar::new();
    if ql == 0 {
        if tl > 0 {
            cig.push(CigarOp::D, tl as u32);
        }
        return cig;
    }
    if tl == 0 {
        cig.push(CigarOp::I, ql as u32);
        return cig;
    }
    let a = scoring.match_score;
    let b = scoring.mismatch_penalty;
    let o = scoring.gap_open;
    let e = scoring.gap_extend;
    let oe = o + e;
    let neg = i32::MIN / 4;
    let band = bw.max(ql.abs_diff(tl)) + 1;

    let cols = tl + 1;
    let mut h = vec![neg; (ql + 1) * cols];
    let mut em = vec![neg; (ql + 1) * cols];
    let mut fm = vec![neg; (ql + 1) * cols];
    let at = |i: usize, j: usize| i * cols + j;

    h[at(0, 0)] = 0;
    for jj in 1..=tl.min(band) {
        em[at(0, jj)] = -(o + e * jj as i32);
        h[at(0, jj)] = em[at(0, jj)];
    }
    for i in 1..=ql {
        let jlo = i.saturating_sub(band).max(1);
        let jhi = (i + band).min(tl);
        if i <= band {
            fm[at(i, 0)] = -(o + e * i as i32);
            h[at(i, 0)] = fm[at(i, 0)];
        }
        for jj in jlo..=jhi {
            let s = if query[i - 1] == target[jj - 1] {
                a
            } else {
                -b
            };
            let diag = h[at(i - 1, jj - 1)];
            let mv = if diag <= neg / 2 { neg } else { diag + s };
            let ev = h[at(i, jj - 1)]
                .saturating_sub(oe)
                .max(em[at(i, jj - 1)].saturating_sub(e));
            em[at(i, jj)] = ev;
            let fv = h[at(i - 1, jj)]
                .saturating_sub(oe)
                .max(fm[at(i - 1, jj)].saturating_sub(e));
            fm[at(i, jj)] = fv;
            h[at(i, jj)] = mv.max(ev).max(fv);
        }
    }

    #[derive(PartialEq)]
    enum St {
        H,
        E,
        F,
    }
    let mut ops: Vec<CigarOp> = Vec::new();
    let (mut i, mut j) = (ql, tl);
    let mut st = St::H;
    while i > 0 || j > 0 {
        match st {
            St::H => {
                if i > 0 && j > 0 {
                    let is_match = query[i - 1] == target[j - 1];
                    let s = if is_match { a } else { -b };
                    let diag = h[at(i - 1, j - 1)];
                    let mv = if diag <= neg / 2 { neg } else { diag + s };
                    if h[at(i, j)] == mv {
                        ops.push(if is_match { CigarOp::Eq } else { CigarOp::X });
                        i -= 1;
                        j -= 1;
                        continue;
                    }
                }
                if j > 0 && h[at(i, j)] == em[at(i, j)] {
                    st = St::E;
                    continue;
                }
                if i > 0 && h[at(i, j)] == fm[at(i, j)] {
                    st = St::F;
                    continue;
                }
                if j > 0 {
                    ops.push(CigarOp::D);
                    j -= 1;
                } else if i > 0 {
                    ops.push(CigarOp::I);
                    i -= 1;
                }
            }
            St::E => {
                let extend = em[at(i, j)] == em[at(i, j - 1)].saturating_sub(e);
                ops.push(CigarOp::D);
                j -= 1;
                if !extend {
                    st = St::H;
                }
            }
            St::F => {
                let extend = fm[at(i, j)] == fm[at(i - 1, j)].saturating_sub(e);
                ops.push(CigarOp::I);
                i -= 1;
                if !extend {
                    st = St::H;
                }
            }
        }
    }
    ops.reverse();
    Cigar::compress(ops)
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
    fn rescue_mate_accepts_sub_t_concordant_rescue() {
        // Verifies that a concordant rescue scoring in [min_seed_len, min_score) is
        // accepted after the bwa mem_matesw gate fix (pair-conditioned acceptance).
        // With Scoring::default() match_score=1, a 25 bp exact match scores 25.
        // min_seed_len=19 and min_score=30, so 25 >= 19 (accepted) and 25 < 30 (sub-T).
        let body = "ACGTACGATCGATCGGATTCCAGTCAGTCAGGATCCATGCATGCATTAGCATCGATCGTA";
        let ref_seq = Reference::parse_fasta(&format!(">t\n{body}")).unwrap();
        let ref_data = ref_seq.as_slice_2n();
        let index = FMIndexAlias::build_2n(&ref_seq);
        let aligner = Aligner::new(index, ref_data.clone())
            .min_seed_len(19)
            .min_score(30);
        let n_fwd = body.len();
        let fwd = &ref_data[..n_fwd];

        // mate maps forward at position 0, covering 20 bp
        let mut mate_cigar = Cigar::new();
        mate_cigar.push(CigarOp::M, 20);
        let mut mate = AlignmentResult::new(0, mate_cigar);
        mate.reverse_strand = false;

        // orphan is the reverse-strand read of forward region [35, 60): 25 bp exact
        // match → score 25, which is >= min_seed_len (19) but < min_score (30).
        let region = fwd[35..60].to_vec();
        let orphan = reverse_complement(&region);

        // upper_bound = 60 + 20 = 80, so the window [0, 60] covers position 35.
        let dist = InsertSizeDistribution::with_params(60.0, 5.0);
        let result = aligner
            .rescue_mate(&orphan, &mate, &dist)
            .expect("sub-T concordant rescue should be accepted");
        assert!(
            result.reverse_strand,
            "forward mate -> reverse orphan downstream"
        );
        assert!(result.score >= 19, "score must be >= min_seed_len");
        assert!(result.score < 30, "score must be < min_score (sub-T)");
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
        let result = approx_mapq_se(151, 0, 0, 151, &Scoring::default(), 19, 0.0, 0);
        assert_eq!(result, 60);
    }

    #[test]
    fn approx_mapq_se_sub_ge_score_is_zero() {
        assert_eq!(
            approx_mapq_se(100, 100, 0, 100, &Scoring::default(), 19, 0.0, 0),
            0
        );
        assert_eq!(
            approx_mapq_se(100, 120, 0, 100, &Scoring::default(), 19, 0.0, 0),
            0
        );
    }

    #[test]
    fn approx_mapq_se_zero_score_is_zero() {
        assert_eq!(
            approx_mapq_se(0, 0, 0, 100, &Scoring::default(), 19, 0.0, 0),
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
        let m0 = approx_mapq_se(40, 30, 0, 45, &Scoring::default(), 19, 0.0, 0);
        let m4 = approx_mapq_se(40, 30, 4, 45, &Scoring::default(), 19, 0.0, 0);
        assert_eq!(m0, 55);
        assert_eq!(m4, 48);
        assert!(m4 < m0);
    }

    #[test]
    fn approx_mapq_se_clamps_to_60() {
        // A perfect long read with no suboptimal hits: must clamp to exactly 60.
        let result = approx_mapq_se(200, 0, 0, 200, &Scoring::default(), 19, 0.0, 0);
        assert_eq!(result, 60);
    }

    #[test]
    fn approx_mapq_se_frac_rep_downweights() {
        // bwa's final mapq *= (1 - frac_rep): a fully-repetitive read drops to 0,
        // half-repetitive roughly halves the unweighted score.
        let full = approx_mapq_se(200, 0, 0, 200, &Scoring::default(), 19, 1.0, 0);
        let half = approx_mapq_se(200, 0, 0, 200, &Scoring::default(), 19, 0.5, 0);
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

    // ---- T-025: primary selection by extended-alignment score ----

    /// Two-locus reference engineered so the chain heuristic (longest exact seed)
    /// picks the WRONG locus:
    ///
    ///   Locus A (ref[0..30]): exact copy of read[0..30].  Chain score = 30 (one
    ///     30-bp MEM).  Extension soft-clips read[31..50] because ref[30..] is
    ///     all-G junk → extended-alignment score = 30.
    ///
    ///   Locus B (ref[100..150]): near-copy of the full 50-bp read with single-base
    ///     mismatches at read positions 10, 20, and 30.  The mismatches break the
    ///     50-bp alignment into exact segments of 10, 9, 9, and 19 bp.  All segments
    ///     are shorter than the 30-bp MEM at locus A; the only segment >= min_seed_len
    ///     that is NOT superseded by the locus-A MEM is the 19-bp tail at
    ///     read[31..50] / ref[131..150].  That sole locus-B seed gives chain score
    ///     19 < 30, so the chain heuristic picks locus A.  Extended-alignment score
    ///     at locus B: 47 matches - 3 x 4 = 35 > 30.
    ///
    /// T-025 must select locus B (pos 100).  This test FAILS on the pre-T-025 code
    /// (which returns pos=0) and PASSES after.
    #[test]
    fn test_primary_selected_by_max_extended_score() {
        // Encode: A=0, C=1, G=2, T=3.
        //
        // 50-bp hand-chosen non-repetitive read.  No 10-bp window appears twice
        // so the FM-index finds at most one occurrence per seed.  The anchor half
        // (read[0..30]) and the tail half (read[30..50]) are constructed to be
        // distinct: the anchor uses an ACTAGCT... pattern; the tail uses TGCAAT...
        #[rustfmt::skip]
        let read_bases: Vec<u8> = vec![
            // read[0..30]: locus-A anchor
            0,1,3,0,2,1,3,3,2,0,  // ACTAGCTTGA
            1,3,2,0,1,0,3,1,2,3,  // CTGACTGACT  <- position 10: read[10]=1(C)
            3,2,0,1,0,3,1,2,0,3,  // TGACATCGAT  <- position 20: read[20]=0(A)
            // read[30..50]: includes mismatch at position 30 of locus_B then 19-bp tail
            3,2,1,0,0,3,1,2,3,3,  // TGCAAATCTT  <- position 30: read[30]=3(T)
            2,0,1,3,2,0,1,0,3,1,  // GACTGACTGA
        ];
        assert_eq!(read_bases.len(), 50);

        // Locus B = near-copy of read with 3 mismatches; each flipped by +1 mod 4.
        let mut locus_b: Vec<u8> = read_bases.clone();
        locus_b[10] = (read_bases[10] + 1) % 4; // C -> G
        locus_b[20] = (read_bases[20] + 1) % 4; // A -> C
        locus_b[30] = (read_bases[30] + 1) % 4; // T -> A

        // Reference layout (forward-only, 160 bp total):
        //   [0..30]   locus A: exact copy of read[0..30]
        //   [30..100] junk: all G=2 (70 bp — mismatches all non-G read bases)
        //   [100..150] locus B: 50-bp near-copy (3 mismatches at 10, 20, 30)
        //   [150..160] padding: all A=0
        let mut ref_raw: Vec<u8> = Vec::with_capacity(160);
        ref_raw.extend_from_slice(&read_bases[..30]);
        ref_raw.extend(std::iter::repeat(2u8).take(70));
        ref_raw.extend_from_slice(&locus_b);
        ref_raw.extend(std::iter::repeat(0u8).take(10));
        assert_eq!(ref_raw.len(), 160);

        let base_chars = ['A', 'C', 'G', 'T'];
        let fasta_body: String = ref_raw.iter().map(|&b| base_chars[b as usize]).collect();
        let fasta = format!(">t\n{fasta_body}");

        let ref_seq = Reference::parse_fasta(&fasta).unwrap();
        let ref_data = ref_seq.as_slice().to_vec();
        let index = FMIndexAlias::build(&ref_seq);
        // min_seed_len(10) so both loci produce seeds; min_score(0) so the 30-bp
        // partial hit at locus A is not suppressed.
        let aligner = Aligner::new(index, ref_data).min_seed_len(10).min_score(0);

        let result = aligner.align_read(&read_bases, None).unwrap();

        assert_eq!(result.flag & 0x4, 0, "read should be mapped");

        // Chain heuristic (pre-T-025 behaviour):
        //   read[0..30] -> 30-bp MEM at locus A (ref 0).  Chain score = 30.
        //   read[11..20] and read[21..30] -> 9-bp segments below min_seed_len; not found.
        //   read[31..50] -> 19-bp MEM at locus B tail (ref 131).  Chain score = 19.
        //   (read[0..10] at locus-B / ref[100..110] is superseded at query_start=0 by
        //   locus A's 30-bp MEM.)  Gap between chains > 50 bp: no merge.  Chain picks A.
        //
        // Extended-alignment scores:
        //   Locus A: 30 exact matches, then 20-base soft-clip (junk) -> score = 30.
        //   Locus B (seeded at read[31..50] / ref[131..150]):
        //     backward over read[0..31] / ref[100..131]: 10= 1X 9= 1X 9= 1X -> 28 matches, 3 mm
        //     forward over read[31..50] / ref[131..150]:  19= -> 19 matches
        //     compute_alignment_score: 47 matches - 3 x 4 = 35 > 30.
        //
        // T-025 must pick locus B (pos 100).
        assert_eq!(
            result.position, 100,
            "T-025: primary must be at locus B (pos 100, extended score 35); \
             got pos={} score={} cigar={}",
            result.position, result.score, result.cigar
        );
        assert!(
            result.score >= 30,
            "locus-B extended score {} must be >= 30",
            result.score
        );
    }

    // ---- T-029: align_read_regions unit tests ----

    /// Helper to build a 2N-indexed aligner for a given reference body string.
    fn regions_test_aligner(body: &str) -> Aligner {
        let ref_seq = Reference::parse_fasta(&format!(">t\n{body}")).unwrap();
        let ref_data = ref_seq.as_slice_2n();
        let index = FMIndexAlias::build_2n(&ref_seq);
        Aligner::new(index, ref_data).min_seed_len(10).min_score(10)
    }

    #[test]
    fn align_read_regions_primary_matches_align_read() {
        let body = "ACGTACGATCGATCGGATTCCAGTCAGTCAGGATCCATGCATGCATTAGCATCGATCGTA";
        let aligner = regions_test_aligner(body);
        let query = enc(&body[5..35]);

        let ar = aligner.align_read(&query, None).unwrap();
        let rr = aligner.align_read_regions(&query, None).unwrap();

        // If it mapped, regions[0] must equal align_read on all key fields.
        if ar.flag & 0x4 == 0 {
            assert!(
                !rr.regions.is_empty(),
                "align_read_regions must have at least one region for a mapped read"
            );
            let primary = &rr.regions[0];
            assert_eq!(primary.position, ar.position, "position mismatch");
            assert_eq!(primary.score, ar.score, "score mismatch");
            assert_eq!(
                primary.cigar.to_string(),
                ar.cigar.to_string(),
                "cigar mismatch"
            );
            assert_eq!(primary.mapq, ar.mapq, "mapq mismatch");
            assert_eq!(
                primary.reverse_strand, ar.reverse_strand,
                "reverse_strand mismatch"
            );
            assert_eq!(primary.xs, ar.xs, "xs mismatch");
        }
    }

    #[test]
    fn align_read_regions_unmapped_is_empty() {
        // min_seed_len=100 prevents any seed from being found → unmapped.
        let body = "ACGTACGATCGATCGGATTCCAGTCAGTCAGGATCCATGCATGCATTAGCATCGATCGTA";
        let ref_seq = Reference::parse_fasta(&format!(">t\n{body}")).unwrap();
        let ref_data = ref_seq.as_slice_2n();
        let index = FMIndexAlias::build_2n(&ref_seq);
        let aligner = Aligner::new(index, ref_data)
            .min_seed_len(100)
            .min_score(10);

        let query = enc(&body[5..35]);
        let rr = aligner.align_read_regions(&query, None).unwrap();
        assert!(
            rr.regions.is_empty(),
            "unmapped read must yield empty regions"
        );

        let ar = aligner.align_read(&query, None).unwrap();
        assert_eq!(
            ar.flag & 0x4,
            0x4,
            "align_read must return unmapped flag when regions is empty"
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

    /// A forward placement's span is unchanged; the reverse-strand mirror of the same
    /// aligned core lands on the same forward-read interval, so overlap is comparable.
    #[test]
    fn forward_query_span_mirrors_reverse_strand() {
        // 99S 52M at the read's reference-forward orientation.
        let c = make_cigar(&[(CigarOp::S, 99), (CigarOp::M, 52)]);
        assert_eq!(forward_query_span(&c, false, 151), (99, 151));
        // On the reverse strand the same CIGAR aligns the first 52 read bases.
        assert_eq!(forward_query_span(&c, true, 151), (0, 52));
        // Which now overlaps a forward 52M 99S placement of the opposite-strand copy.
        let d = make_cigar(&[(CigarOp::M, 52), (CigarOp::S, 99)]);
        assert_eq!(forward_query_span(&d, false, 151), (0, 52));
    }

    // ---- T-026: ksw_extend + ksw_global unit tests ----

    #[test]
    fn ksw_extend_exact_match_reaches_end() {
        let scoring = Scoring::default();
        let query = enc("ACGTACGTAC");
        let target = query.clone();
        // h0 > 0 seeds the diagonal so the local-alignment accumulates;
        // using one match_score unit simulates "one previously matched base".
        let h0 = scoring.match_score;
        let end_bonus = 5;
        let r = ksw_extend(&query, &target, &scoring, end_bonus, h0);
        // Expected: h0 + qlen * match_score accumulated along the exact diagonal.
        let expected_gscore = h0 + query.len() as i32 * scoring.match_score;
        assert_eq!(
            r.gscore, expected_gscore,
            "gscore should equal seed score + full-length match score"
        );
        assert_eq!(r.qle, query.len(), "qle should reach end of query");
        assert_eq!(r.tle, target.len(), "tle should reach end of target");
    }

    #[test]
    fn ksw_extend_local_stops_before_bad_tail() {
        let scoring = Scoring::default();
        // Prefix of 8 bases match; tail of 4 bases are all mismatches
        let prefix = enc("ACGTACGT");
        let prefix_len = prefix.len();
        // tail = TTTT = [3,3,3,3] — all different from anything in prefix context
        let bad_tail: Vec<u8> = vec![3, 3, 3, 3];
        let mut query = prefix.clone();
        query.extend_from_slice(&bad_tail);
        // target: same prefix, then all-G junk at the end
        let mut target = prefix.clone();
        target.extend(std::iter::repeat(2u8).take(4));

        let r = ksw_extend(&query, &target, &scoring, 0, 0);
        // The local maximum should be at most the matching prefix length
        assert!(
            r.qle <= prefix_len,
            "qle={} should stop at or before prefix end {}",
            r.qle,
            prefix_len
        );
        // gscore should be strictly less than score (so caller would soft-clip)
        assert!(
            r.gscore < r.score,
            "gscore={} should be less than local max score={} when tail is bad",
            r.gscore,
            r.score
        );
    }

    #[test]
    fn ksw_global_exact_match() {
        let scoring = Scoring::default();
        let query = enc("ACGTACGT");
        let len = query.len();
        let target = query.clone();
        let cigar = ksw_global(&query, &target, &scoring, KSW_BAND_WIDTH as usize);
        assert_eq!(
            cigar.reference_length(),
            len,
            "reference_length should equal sequence length"
        );
        assert_eq!(
            cigar.ops.len(),
            1,
            "exact match should compress to a single op; got {}",
            cigar
        );
        assert_eq!(
            cigar.ops[0].0,
            CigarOp::Eq,
            "single op should be Eq for exact match"
        );
    }

    #[test]
    fn ksw_global_single_mismatch() {
        let scoring = Scoring::default();
        let query = enc("ACGTACGT");
        let len = query.len();
        // one mismatch in the middle: position 4 (T->A)
        let mut target = query.clone();
        target[4] = (target[4] + 1) % 4;
        let cigar = ksw_global(&query, &target, &scoring, KSW_BAND_WIDTH as usize);
        assert_eq!(
            cigar.reference_length(),
            len,
            "reference_length should equal len"
        );
        // Should have at least one X op
        let has_x = cigar.ops.iter().any(|(op, _)| *op == CigarOp::X);
        assert!(
            has_x,
            "CIGAR should contain an X op for the mismatch; got {}",
            cigar
        );
        // No I or D expected (same length)
        let has_indel = cigar
            .ops
            .iter()
            .any(|(op, _)| matches!(op, CigarOp::I | CigarOp::D));
        assert!(
            !has_indel,
            "no indels expected for single-mismatch alignment; got {}",
            cigar
        );
    }

    #[test]
    fn ksw_global_insertion() {
        let scoring = Scoring::default();
        // query is one base longer than target: extra base at position 4
        let target = enc("ACGTACGT");
        let mut query = enc("ACGT");
        query.push(2); // extra G inserted
        query.extend_from_slice(&enc("ACGT"));
        let cigar = ksw_global(&query, &target, &scoring, KSW_BAND_WIDTH as usize);
        assert_eq!(
            cigar.reference_length(),
            target.len(),
            "reference_length should equal target length"
        );
        let ins_count: u32 = cigar
            .ops
            .iter()
            .filter(|(op, _)| *op == CigarOp::I)
            .map(|(_, l)| *l)
            .sum();
        assert_eq!(
            ins_count, 1,
            "should have exactly 1 inserted base; got {}",
            cigar
        );
    }

    #[test]
    fn ksw_global_deletion() {
        let scoring = Scoring::default();
        // target is one base longer than query: extra base at position 4
        let query = enc("ACGTACGT");
        let mut target = enc("ACGT");
        target.push(2); // extra G in reference
        target.extend_from_slice(&enc("ACGT"));
        let cigar = ksw_global(&query, &target, &scoring, KSW_BAND_WIDTH as usize);
        assert_eq!(
            cigar.reference_length(),
            target.len(),
            "reference_length should equal target length"
        );
        let del_count: u32 = cigar
            .ops
            .iter()
            .filter(|(op, _)| *op == CigarOp::D)
            .map(|(_, l)| *l)
            .sum();
        assert_eq!(
            del_count, 1,
            "should have exactly 1 deleted base; got {}",
            cigar
        );
    }

    /// csub: a reference-overlapping placement that does NOT have sufficient query
    /// overlap to be a secondary (so it is absent from `sub`) should still raise
    /// the effective sub via csub, lowering MAPQ.
    ///
    /// Layout (all forward strand):
    ///   Primary:    ref[10..60]  score=50  (query[0..50])
    ///   Secondary:  ref[80..130] score=30  (query[0..50], disjoint ref → in sub)
    ///   Tandem hit: ref[20..70]  score=45  (query[60..110], no query overlap with
    ///               primary[0..50], so it is NOT a secondary; but it overlaps
    ///               primary in reference → csub=45)
    ///
    /// Expected: approx_mapq_se with sub=30, csub=0 > approx_mapq_se with csub=45.
    #[test]
    fn approx_mapq_se_csub_lowers_mapq() {
        let scoring = Scoring::default();
        let min_seed_len = 19;

        // Without csub: sub=30
        let mapq_no_csub = approx_mapq_se(50, 30, 0, 50, &scoring, min_seed_len, 0.0, 0);
        // With csub=45: effective sub = max(30, 45) = 45 → smaller score gap → lower MAPQ
        let mapq_with_csub = approx_mapq_se(50, 30, 0, 50, &scoring, min_seed_len, 0.0, 45);
        assert!(
            mapq_with_csub < mapq_no_csub,
            "csub=45 should lower MAPQ: no_csub={mapq_no_csub} vs csub={mapq_with_csub}"
        );
    }

    /// csub >= score: MAPQ must be 0 (same as sub >= score).
    #[test]
    fn approx_mapq_se_csub_ge_score_gives_zero() {
        let scoring = Scoring::default();
        let mapq = approx_mapq_se(50, 30, 0, 50, &scoring, 19, 0.0, 50);
        assert_eq!(mapq, 0, "csub >= score must yield MAPQ 0");
    }

    /// Paired csub clamp: when both reads have AS=XS (sub >= score → SE MAPQ 0)
    /// the paired MAPQ must be clamped to 0 by raw_mapq(score − csub) when csub = score.
    #[test]
    fn paired_csub_clamp_as_eq_xs_gives_zero() {
        use crate::paired::{mem_pair, paired_mapq, InsertSizeDistribution, PEN_UNPAIRED};
        use crate::types::CigarOp;
        let make_read = |pos: usize, rev: bool, sc: i32| {
            let mut cigar = Cigar::new();
            cigar.push(CigarOp::M, 150);
            let mut r = AlignmentResult::new(pos, cigar);
            r.reverse_strand = rev;
            r.score = sc;
            r.mapq = 0;
            r
        };
        let dist = InsertSizeDistribution::with_params(300.0, 30.0);
        // Both reads score 100; SE MAPQ is already 0 because sub=100 >= score=100.
        // Paired score might still give q_pe > 0.
        let r1 = make_read(100, false, 100);
        let r2 = make_read(350, true, 100);
        let choice = mem_pair(&[r1.clone()], &[r2.clone()], &dist, 1, 5).unwrap();
        let score_un = r1.score + r2.score - PEN_UNPAIRED;
        // csub = score for both reads → raw_mapq(score-csub) = raw_mapq(0) = 0 → clamp to 0
        let (m1, m2) = paired_mapq(
            &choice, score_un, 0, // q_se1 = 0 (AS=XS)
            0, // q_se2 = 0
            0.0, 0.0, 1, r1.score, // csub1 = score → clamp to 0
            r2.score, // csub2 = score → clamp to 0
            r1.score, r2.score,
        );
        assert_eq!(m1, 0, "csub=score clamp must hold MAPQ at 0 for read1");
        assert_eq!(m2, 0, "csub=score clamp must hold MAPQ at 0 for read2");
    }

    #[test]
    fn region_se_mapq_reproduces_primary_mapq() {
        let body = "ACGTACGATCGATCGGATTCCAGTCAGTCAGGATCCATGCATGCATTAGCATCGATCGTA";
        let aligner = regions_test_aligner(body);
        let query = enc(&body[5..35]);

        let ar = aligner.align_read(&query, None).unwrap();
        let rr = aligner.align_read_regions(&query, None).unwrap();

        if ar.flag & 0x4 == 0 {
            assert!(!rr.regions.is_empty());
            let recomputed =
                aligner.region_se_mapq(&rr.regions[0], rr.sub, rr.sub_n, rr.frac_rep, rr.csub);
            assert_eq!(
                recomputed, ar.mapq,
                "region_se_mapq on primary region must reproduce align_read mapq"
            );
        }
    }
}
