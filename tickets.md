# Ticket Series: bwa-rs vs bwa-mem2 Output Divergence

Tracking the investigation into why `bwa-rs` SAM output is not byte-identical to
upstream **bwa-mem2** on the bwa-mem2 performance dataset (D2 / SRR7733443
against `human_g1k_v37.fasta`).

**Status:** In progress (2026-06-01). Two blockers and one OOM resolved
(T-002, T-003, T-004 — see PR `fix-oom-and-64bit`); bwa-rs can now index the full
human genome and map both strands. Soft-clipping (T-005), a minimum-score
filter (T-017), QNAME trimming (T-006), and paired-end SAM fields (T-007) are now
also resolved. Remaining output-divergence gaps are ticketed below, ordered
roughly by impact toward byte-identical SAM.

**Verification baseline:** chr1 (249 Mbp), 300 uniquely-mapping read pairs from
SRR7733443, normalized QNAMEs, `-k 19`. POS field matches:
- After T-004: **246/300**.
- After T-005 (soft-clipping): **284/300**.

POS is unchanged by T-007. T-007 instead closed the pairing-field divergence: on
this set FLAG dropped 100% → 1.3% (8/600), RNEXT 100% → 0%, TLEN 98.7% → 1.3%.

The residual **16/300** POS mismatches break down into three independent causes
(none inside soft-clipping itself):
- **~10** — rigid seed anchoring: bwa-mem2 runs banded Smith-Waterman that can
  shift the alignment 1–10 bp relative to the seed (placing a read full-length
  where bwa-rs pins the seed and soft-clips the overhang). See **T-018**.
- **3** — low-scoring partial hits bwa-mem2 rejects. T-017 reports these unmapped
  (matching bwa-mem2's mapped/unmapped status); T-007 now also emits their
  unmapped-mate POS/FLAG convention (mapped mate's coordinate + paired bits).
- **3** — reads bwa-mem2 maps via mate rescue that bwa-rs leaves unmapped. Split
  out of T-007 to **T-019** (the 8/600 residual FLAG diffs are these rescue cases).

---

## How to use this file

- One ticket per **distinct root cause** of divergence (not per differing record).
- Order tickets by impact: how many SAM records / fields each cause accounts for.
- Status values: `OPEN` · `INVESTIGATING` · `FIXED` · `WONTFIX` · `BLOCKED`.

### Ticket template

```
### T-NNN: <short title>
- **Status:** OPEN
- **Severity:** blocker | high | medium | low
- **Affected field(s):** e.g. CIGAR, POS, MAPQ, @PG, FLAG
- **Symptom:** what differs, with a concrete example record from both outputs.
- **Suspected cause:** hypothesis (link to code: file:line).
- **Investigation:** findings as they accrue.
- **Resolution:** what fixed it / why deferred.
```

---

## Anticipated tickets (from `docs/PERFORMANCE.md`, pre-comparison)

These are hypotheses to confirm/refute once real diffs are observed:

- **@PG / header lines** differ (program ID, version, command line). Likely first diff.
- **CIGAR operators**: bwa-rs emits `=`/`X`; bwa-mem2 emits `M`. Pervasive.
- **min_seed_len default**: bwa-rs 10 vs bwa-mem2 19 — controllable via `-k 19`.
- **gap_open penalty**: bwa-rs 6 vs bwa-mem2 5 — changes scoring/alignment choice.
- **MAPQ** computation differences.
- Algorithmic differences (Occ wavelet tree, MEM finding) → seed/chain divergence.

---

## Tickets

### T-001: Index build aborts on IUPAC ambiguity codes in reference FASTA
- **Status:** OPEN (worked around for the benchmark via reference sanitization)
- **Severity:** blocker
- **Affected field(s):** N/A — fails before any SAM is produced (index build).
- **Symptom:** `bwa-mem index -r test-data/human_g1k_v37.fasta` exits with
  `Error: Parse("Invalid base: M")` after ~14s. The stock bwa-mem2 benchmark
  reference `human_g1k_v37.fasta` contains 3 IUPAC ambiguity codes
  (`1×M`, `2×R`) among its ~3.1 Gbp. bwa-rs cannot build an index from it at all.
- **Suspected cause:** `src/reference.rs:78` `encode_base()` (used by the index
  builder via `encode_bases()` at `src/reference.rs:74`) returns `Err` for any
  character outside `A/C/G/T/N`. Upstream bwa / bwa-mem2 map every non-ACGT base
  to `N` (`nst_nt4_table`) and record ambiguous runs in the `.amb`/`.ann`
  sidecars; they never error.
- **Investigation:**
  - Confirmed exactly 3 offending positions genome-wide via
    `grep -v '^>' ... | grep -oP '[^ACGTNacgtn]' | sort | uniq -c` → `1 M`, `2 R`.
  - Note the inconsistency: bwa-rs's *read* path already tolerates this —
    `src/utils.rs:18` `encode_sequence()` maps any unknown byte to `4` (N). Only
    the reference-encoding path (`reference.rs::encode_base`) is strict.
- **Resolution:** Recommended fix — make `reference.rs::encode_base` map unknown
  bases to `N` (value `0b100`) instead of erroring, matching both bwa-mem2 and
  bwa-rs's own `utils::encode_sequence`. For this benchmark run the reference was
  sanitized instead (M/R→N, identical bytes otherwise) into
  `test-data/human_g1k_v37.clean.fasta`, fed to **both** aligners so the
  comparison stays fair (bwa-mem2 maps these to N internally regardless). Released
  binary left unpatched.

### T-002: Cannot index references larger than ~2.147 Gbp (32-bit suffix array)
- **Status:** OPEN — hard blocker, no input-level workaround.
- **Severity:** blocker
- **Affected field(s):** N/A — fails during index build; **the full human-genome
  byte-identity benchmark is impossible with current bwa-rs.**
- **Symptom:** `bwa-mem index` on `human_g1k_v37.fasta` (3.1 Gbp, sanitized per
  T-001) panics after ~55s / ~8 GB RSS:
  ```
  thread 'main' panicked at libsais-rs-0.1.1/src/lib.rs:9299:39:
  input length must fit SaSint: TryFromIntError(())
  ```
- **Suspected cause:** The suffix-array / FM-index pipeline is 32-bit end to end:
  - `src/sa.rs:17,29` allocate the SA as `vec![0i32; n]` and call
    `libsais_rs::libsais` (the 32-bit entry point), which rejects any input with
    `len > i32::MAX` (2,147,483,647 ≈ 2.147 Gbp).
  - The SA is then stored as `Vec<u32>` (`sa.rs:20,31`) and serialized as 4-byte
    LE (`sa.rs:92-95`, read back at `sa.rs:98-106`; `compact.rs:483` reads a u32
    `bwt_len`). So even the on-disk format caps positions at u32 (~4.29 Gbp).
  - The SA is built over the forward sequence only (N), not forward+RC (2N) —
    `compact.rs:236`, `fm_index.rs:104-107`. So the binding limit is the i32
    libsais cap on N: **~2.147 Gbp**. The human genome's 3.1 Gbp forward sequence
    exceeds it outright. (For contrast, bwa-mem2's own log reports
    `ref seq len = 6203609478` — it indexes the 2N = 6.2 Gbp concatenation with
    64-bit indices.)
- **Investigation:** Reproduced deterministically. Confirmed via source that no
  64-bit code path exists (`libsais_rs::libsais` / `libsais_int`, never a
  `libsais64` variant).
- **Resolution:** **FIXED** (PR `fix-oom-and-64bit`). Upgraded `libsais-rs` 0.1→0.2
  and switched SA construction to the 64-bit entry point `libsais64` /
  `libsais64_omp` (`src/sa.rs`); widened the SA, `f_column`, occ counts, and the
  on-disk index format from `u32` to `u64` (`src/sa.rs`, `src/fm_index.rs`,
  `src/compact.rs`; index format VERSION 3→4). The full human genome (3.1 Gbp,
  6.2 Gbp as 2N) now indexes without panic. The `mmap_index` aux format remains
  32-bit (≤4.29 Gbp) by design; see T-012.

---

## Resolved in PR `fix-oom-and-64bit`

### T-003: `mem` exhausts >128 GB RAM and is OOM-killed during alignment
- **Status:** FIXED
- **Severity:** blocker
- **Affected field(s):** N/A — process killed before producing records.
- **Symptom:** Aligning even 2,000 read pairs to chr1 (1 thread) was SIGKILL'd
  after spiking past the 128 GB cgroup limit; only the SAM header was written.
- **Root cause:** `affine_extend_forward` (`src/alignment.rs`) allocated a *dense*
  `(query_len+1) × (ref_len+1)` DP matrix where `ref_len` was the entire remaining
  reference (up to the whole chromosome), even though the DP loop is banded — so a
  single forward extension near a chromosome start allocated tens-to-hundreds of GB.
- **Resolution:** Cap the matrix width to the band (`query_len + bandwidth`); the
  banded DP never touches columns beyond that, so results are unchanged. chr1
  alignment dropped from >128 GB (killed) to ~8.5 GB.

### T-004: No reverse-strand alignment; N-only index (does not match bwa-mem2 2N)
- **Status:** FIXED
- **Severity:** high
- **Affected field(s):** FLAG (0x10), POS, RNAME presence.
- **Symptom:** The CLI only searched the forward query, so reverse-strand reads
  were left unmapped or misplaced; the index was built over N (forward only) while
  bwa / bwa-mem2 index 2N (`ref seq len = 6203609478`).
- **Resolution:** Added `Reference::as_slice_2n()` (forward ++ reverse-complement)
  and `FMIndex::build_2n` + an `n_fwd` field (`src/reference.rs`, `src/fm_index.rs`);
  the CLI now builds the 2N index and extends in 2N space, mapping SA positions
  `≥ n_fwd` back to forward coordinates with the 0x10 flag set and the CIGAR
  reversed (`src/alignment.rs`, `src/main.rs`). **Strand now matches bwa-mem2
  297/297** on the chr1 verification set.

---

## Remaining gaps toward byte-identical SAM (OPEN)

Ordered roughly by impact. Verification set: chr1, 300 unique read pairs.

### T-005: No soft-clipping — end-to-end alignment shifts POS and CIGAR
- **Status:** FIXED (PR stacked on `fix-oom-and-64bit`, branch `t005-soft-clipping`)
- **Severity:** high (largest remaining POS/CIGAR driver)
- **Affected field(s):** POS, CIGAR.
- **Symptom:** bwa-rs forced a full-length alignment with `=`/`X` runs
  (e.g. `43=7X1=14X…`) where bwa-mem2 soft-clips poor ends (`33S66M52S`). POS then
  started further right. Accounted for the bulk of the POS mismatches (off by 1–40 bp).
- **Cause:** Seed extension (`src/alignment.rs build_alignment` /
  `affine_extend_forward`) always extended to the read ends instead of stopping at
  the max-scoring position and emitting `S` for the trimmed tail.
- **Resolution:** Added `clip_penalty` (default 5) to `Scoring`. `affine_extend_forward`
  tracks the best query-end-reaching score (`gscore`) alongside the local max and
  extends to the read end only when that costs less than `clip_penalty`; otherwise it
  stops at the max-scoring cell. `extend_seed_backward` reuses the forward DP on
  reversed slices so the 5′ end clips identically. `build_alignment` emits
  leading/trailing `S` so the CIGAR query length equals the read length.
  **POS 246 → 284/300.** Soft-clip lengths now match bwa-mem2 wherever the placement
  agrees (e.g. `4S147M` ↔ `4S8=1X19=…`, `98M53S` ↔ `80=1X17=53S`).
- **Note:** The residual `=`/`X`-vs-`M` CIGAR string diff is **T-010**; the residual
  1–10 bp POS placement diffs are **T-018** (banded SW), not soft-clipping.

### T-006: QNAME not truncated at first whitespace
- **Status:** FIXED (branch `t006-qname-trim`, stacked on `t017-min-score`)
- **Severity:** high (every record differs on QNAME)
- **Affected field(s):** QNAME.
- **Symptom:** bwa-rs emits `SRR7733443.100000 100000 length=151`; bwa-mem2 emits
  `SRR7733443.100000`. The FASTQ reader keeps the full header line as the name.
- **Cause:** `src/fastq.rs` record-name parsing kept the whole line.
- **Resolution:** Added `trim_qname` to `src/fastq.rs`: truncate the read name at
  the first ASCII whitespace, then strip a trailing `/`+single-digit read-number
  suffix (bwa `trim_readno`: `/1`/`/2`/etc.). Applied where `qname` is built in
  `next()`, so both the SAM QNAME and the `Sequence` name are trimmed. QNAME now
  matches bwa-mem2 on every record.

### T-007: Paired-end not implemented (pairing FLAG bits, RNEXT/PNEXT/TLEN)
- **Status:** FIXED (branch `t007-paired-end`, PR #52)
- **Severity:** high
- **Affected field(s):** FLAG (0x1/0x2/0x8/0x10/0x20/0x40/0x80), RNEXT, PNEXT, TLEN.
- **Symptom:** Mates were aligned independently (`align_single` per read,
  `src/main.rs`); bwa-rs emitted FLAG `144` where bwa-mem2 emits `147` (missing
  0x1 paired / 0x2 proper-pair), and RNEXT/PNEXT/TLEN were always `*`/0/0.
- **Cause:** No pairing/insert-size logic; the CLI loop never considered the mate.
- **Resolution:** Added pure pairing helpers in `src/paired.rs` — `is_proper_pair`,
  `template_length`, `fr_insert_size`/`pair_insert_size`, and `mate_fields`
  (assembling FLAG 0x1/0x2/0x8/0x10/0x20/0x40/0x80, RNEXT, PNEXT, TLEN per the bwa
  convention, plus the unmapped-mate coordinate convention where an unmapped read
  inherits its mapped mate's POS). The CLI (`src/main.rs`) now runs a two-phase
  paired loop: align all pairs and estimate the insert-size distribution, then
  emit records with proper-pair flags (output is buffered via `BufWriter`).
  Proper-pair detection classifies FR orientation by each read's 5′ end (so small
  dovetails stay FR) and gates on a mean±4σ insert window (bwa `mem_pestat`
  MAX_STDDEV). `write_paired_record` in `src/sam.rs` emits the record; added
  `AlignmentResult::is_unmapped`.
  **Measured (chr1 / 300 uniq pairs vs bwa-mem2):** FLAG 100% → 1.3% (8/600),
  RNEXT 100% → 0%, PNEXT 100% → 2.8%, TLEN 98.7% → 1.3%; proper-pair FLAG
  mismatches 18 → 0.
- **Deferred:** Mate rescue (the 8/600 residual FLAG diffs — bwa rescues an
  unmapped mate via Smith-Waterman) is split to **T-019**. Batch-wise insert-size
  estimation matching bwa exactly (this uses a single global distribution and
  buffers all alignment results) is a follow-up. RNAME stays `"ref"` (T-008) and
  SEQ/QUAL stay `*` (T-009).

### T-008: RNAME hardcoded to `"ref"`; no multi-contig coordinate mapping
- **Status:** FIXED (branch `t008-rname-contig`)
- **Severity:** high
- **Affected field(s):** RNAME, POS, RNEXT, PNEXT, TLEN (multi-contig).
- **Symptom:** Every mapped record reported RNAME `ref` instead of the contig name
  (`1`); the @SQ header already uses the real name, so output was internally
  inconsistent and invalid for multi-contig references.
- **Cause:** `AlignmentResult.position` is a global offset into the forward
  concatenation of all contigs, but both record writers hardcoded RNAME `"ref"`
  (`src/main.rs write_sam_record`, `src/sam.rs write_paired_record`) and emitted
  POS/PNEXT as the raw global offset — there was no global→contig mapping.
- **Resolution:** Added `Reference::locate(pos) -> Option<(&str, usize)>`
  (`src/reference.rs`), which walks cumulative contig lengths to map a global
  forward position to (contig name, per-contig 0-based offset). Threaded
  `&Reference` into `write_sam_record` and `write_paired_record`; both now emit
  the real contig name and per-contig 1-based POS. `write_paired_record` also
  resolves the mate from its global PNEXT, emitting RNEXT `=` only when the mate
  shares the read's contig (else the mate's contig name) and zeroing TLEN across
  contigs. Pairing logic in `src/paired.rs` stays contig-agnostic — all contig
  resolution lives in the SAM writers. On the chr1 / 300-uniq verification set
  RNAME now matches bwa-mem2 (`1`); POS/PNEXT are unchanged there (single contig,
  offset 0).

### T-009: SEQ and QUAL emitted as `*`
- **Status:** OPEN
- **Severity:** high
- **Affected field(s):** SEQ, QUAL.
- **Symptom:** Every record has SEQ `*` and QUAL `*`; bwa-mem2 emits the read
  sequence (reverse-complemented on the 0x10 strand) and qualities.
- **Suspected cause:** `src/main.rs write_sam_record` hardcodes `"*"`.
- **Resolution:** Emit the read sequence/qualities (revcomp + reversed QUAL on the
  reverse strand; hard-clipped bases excluded per CIGAR).

### T-010: CIGAR uses `=`/`X` instead of `M`
- **Status:** OPEN
- **Severity:** medium
- **Affected field(s):** CIGAR.
- **Symptom:** bwa-rs emits `151=` / `43=7X…`; bwa-mem2 collapses matches and
  mismatches to `M` (`151M`).
- **Suspected cause:** The aligner builds `Eq`/`X` ops and does not collapse them.
- **Resolution:** Collapse `=`/`X` into `M` for output (or add a CIGAR mode);
  trivial post-processing of `Cigar`.

### T-011: Optional tags missing (NM, MD, AS, XS, …)
- **Status:** OPEN
- **Severity:** medium
- **Affected field(s):** optional tag columns.
- **Symptom:** CLI records carry no tags; bwa-mem2 emits `NM:i:`, `MD:Z:`,
  `AS:i:`, `XS:i:`, `MC:Z:`, etc. (`AlignmentResult` already computes `nm` and an
  MD string, but `write_sam_record` drops them.)
- **Resolution:** Emit at least NM/MD/AS (already computed) and XS; wire the
  remaining tags as the corresponding features land.

### T-012: `index` subcommand is a no-op; `mem` rebuilds the index every run
- **Status:** OPEN
- **Severity:** medium (usability/perf, not a SAM-content diff)
- **Symptom:** `bwa-mem index -p PREFIX` builds the FM-index, prints
  `Indexed N bases`, and writes **no files**; `mem` rebuilds the whole index in
  memory on every invocation (`src/main.rs`). `FMIndex::save/load` exist and are
  now 64-bit, but the CLI doesn't call them.
- **Resolution:** Persist the index in `index` (call `save`) and load it in `mem`
  (`load`) instead of rebuilding. Note `mmap_index` is still a separate 32-bit
  format (≤4.29 Gbp); unify or document.

### T-013: CLI alignment is single-threaded
- **Status:** OPEN
- **Severity:** medium (perf)
- **Symptom:** `-t` only sizes the rayon pool used for SA construction; the read
  loop in `run_mem` (`src/main.rs`) aligns one read at a time. bwa-mem2 aligned
  200k pairs in ~9 s on 16 threads.
- **Resolution:** Batch reads and use the existing `ParallelAligner::align_batch`
  (rayon `par_iter`) instead of the sequential loop.

### T-014: Scoring/seed defaults differ from bwa-mem2
- **Status:** OPEN
- **Severity:** medium (changes which alignment wins → POS/CIGAR/MAPQ)
- **Symptom:** `gap_open` 6 vs 5; default `min_seed_len` 10 vs 19 (mitigated with
  `-k 19`). Documented in `docs/PERFORMANCE.md`.
- **Resolution:** Match bwa-mem2 defaults (gap_open 5, min_seed_len 19) or expose
  matching flags; needed for scoring-tie-break parity.

### T-015: MAPQ computation differs
- **Status:** OPEN
- **Severity:** medium
- **Affected field(s):** MAPQ.
- **Symptom:** `calculate_mapq` (`src/alignment.rs`) uses a second-best-seed ratio
  heuristic, not bwa-mem2's `mem_approx_mapq_se` formula (based on sub-optimal
  score, seed coverage, etc.). Matches on the easy MAPQ-60 unique set but will
  diverge elsewhere.
- **Resolution:** Port bwa's MAPQ formula.

### T-016: Header lines differ (@HD, @PG)
- **Status:** OPEN
- **Severity:** low
- **Affected field(s):** @HD, @PG (and @SQ ordering for multi-contig).
- **Symptom:** bwa-rs emits `@HD VN:1.0 SO:unsorted` and
  `@PG ID:bwa-rs PN:bwa-rs VN:0.1.0`; bwa-mem2 emits no @HD and
  `@PG ID:bwa PN:bwa-mem2 VN:2.2.1 CL:<command line>`. (@SQ already matches.)
- **Resolution:** Match bwa-mem2's header emission (drop @HD or match SO; emit a
  bwa-style @PG with the real command line). Note byte-identical @PG also requires
  matching the recorded command line.

### T-017: No minimum-score filter — low-scoring partial hits reported instead of unmapped
- **Status:** FIXED (branch `t017-min-score`, stacked on `t005-soft-clipping`)
- **Severity:** medium
- **Affected field(s):** FLAG (0x4), CIGAR, POS, RNAME (mapped-vs-unmapped status).
- **Symptom:** bwa-rs reported a low-quality partial alignment (e.g. a ~20 bp
  aligned core buried in soft-clips, `83S1X22=1X44S`, score ≈ 14) at a junk locus
  where bwa-mem2 reports the read unmapped. 3/300 of the chr1 verification reads.
- **Cause:** `align_read` returned `build_alignment`'s result unconditionally; bwa
  drops alignments scoring below `-T` (default 30).
- **Resolution:** Added `min_score` (default `DEFAULT_MIN_SCORE = 30`) to `Aligner`
  with a builder and a `-T` CLI flag; `align_read` returns `unmapped_result()` when
  the best alignment scores below it. The 3 reads now report unmapped, matching
  bwa-mem2's mapped/unmapped status and `*` CIGAR.
- **Note:** Their residual POS/FLAG fields still differ because bwa-mem2 assigns an
  unmapped read its mapped mate's POS and the paired FLAG bits — that is the
  unmapped-mate convention owned by **T-007**, not this ticket.

### T-018: Rigid seed anchoring prevents seed-relative shift (needs banded Smith-Waterman)
- **Status:** OPEN
- **Severity:** high (largest remaining POS driver after T-005)
- **Affected field(s):** POS, CIGAR.
- **Symptom:** bwa-mem2 reports a read full-length and shifted 1–10 bp relative to
  bwa-rs (e.g. bwa-mem2 `151M` at POS p−1 where bwa-rs emits `1S150=` at p, or
  bwa-mem2 `15S136M` where bwa-rs keeps a boundary mismatch). ~10/300 of the chr1
  verification reads; POS off by 1–10 bp.
- **Cause:** `build_alignment` treats the seed as immovable — it forces the seed
  bases as `Eq` and only extends outward from the seed boundaries
  (`extend_seed_backward` / `affine_extend_forward`). bwa-mem2 uses the seed only to
  center a DP band and runs full banded Smith-Waterman over the whole read, so the
  optimal alignment can shift within the band. A read that scores higher placed 1 bp
  left (aligning full-length) is found by bwa-mem2 but not by the seed-pinned
  extension.
- **Investigation:** Confirmed by instrumenting the extension on a 6 kb window:
  affected reads have a correct local extension, but the seed anchor is 1–10 bp off
  the placement bwa-mem2's global-in-band alignment picks. A zero-length-baseline
  tweak to the extension fixed the all-junk-tail subset but regressed the
  full-length subset (16 → 20), confirming the divergence is the rigid anchor, not
  the clip arithmetic.
- **Resolution:** Replace seed-pinned outward extension with a single banded
  Smith-Waterman over the read, centered on the seed/chain (bwa `ksw` style); related
  to scoring defaults (T-014) and MAPQ (T-015).

### T-019: Mate rescue not implemented (unmapped mate via banded Smith-Waterman)
- **Status:** OPEN (split out of T-007)
- **Severity:** medium
- **Affected field(s):** FLAG (0x4/0x8/0x2), POS, CIGAR — mapped-vs-unmapped status of a mate.
- **Symptom:** When one mate maps and the other does not, bwa-mem2 runs a banded
  Smith-Waterman around the mapped mate's expected position (`mem_matesw`) and
  often rescues the unmapped mate, marking the pair proper. bwa-rs leaves it
  unmapped. On the chr1 / 300-uniq set this is the entire 8/600 residual FLAG
  divergence after T-007 (bwa-rs `133`/`73`/`89`/`165` → bwa-mem2 `147`/`99`/`83`/`163`).
- **Suspected cause:** No `mem_matesw` equivalent; `align_read` only seeds from the
  read's own MEMs and never aligns a mate into the partner's insert-size window.
- **Resolution:** Port bwa's `mem_matesw` — align the unmapped mate within the
  insert-size window (the `InsertSizeDistribution` from T-007) around the mapped
  mate; accept if it clears the score threshold (`-T`). Related to T-018 (banded SW).
