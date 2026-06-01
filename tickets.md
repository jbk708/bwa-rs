# Ticket Series: bwa-rs vs bwa-mem2 Output Divergence

Tracking the investigation into why `bwa-rs` SAM output is not byte-identical to
upstream **bwa-mem2** on the bwa-mem2 performance dataset (D2 / SRR7733443
against `human_g1k_v37.fasta`).

**Status:** In progress (2026-06-01). Two blockers and one OOM resolved
(T-002, T-003, T-004 — see PR `fix-oom-and-64bit`); bwa-rs can now index the full
human genome and map both strands. Remaining output-divergence gaps are ticketed
below (T-005…T-016), ordered roughly by impact toward byte-identical SAM.

**Verification baseline:** chr1 (249 Mbp), 300 uniquely-mapping read pairs from
SRR7733443. After T-004, strand flag (0x10) matches bwa-mem2 **297/297**; POS
matches **246/297** (remainder driven by T-005 soft-clipping).

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
- **Status:** OPEN
- **Severity:** high (largest remaining POS/CIGAR driver)
- **Affected field(s):** POS, CIGAR.
- **Symptom:** bwa-rs forces a full-length alignment with `=`/`X` runs
  (e.g. `43=7X1=14X…`) where bwa-mem2 soft-clips poor ends (`33S66M52S`). POS then
  starts further right. Accounts for the 51/297 POS mismatches (off by 1–40 bp).
- **Suspected cause:** Seed extension (`src/alignment.rs build_alignment` /
  `affine_extend_forward`) always extends to the read ends instead of stopping at
  the max-scoring position and emitting `S` for the trimmed tail (bwa's local /
  Smith-Waterman behavior with a clip penalty).
- **Resolution:** Implement local-extension clipping: track the best-scoring cell,
  truncate the alignment there, and emit leading/trailing `S` ops.

### T-006: QNAME not truncated at first whitespace
- **Status:** OPEN
- **Severity:** high (every record differs on QNAME)
- **Affected field(s):** QNAME.
- **Symptom:** bwa-rs emits `SRR7733443.100000 100000 length=151`; bwa-mem2 emits
  `SRR7733443.100000`. The FASTQ reader keeps the full header line as the name.
- **Suspected cause:** `src/fastq.rs` record-name parsing keeps the whole line.
- **Resolution:** Trim the read name at the first whitespace (and strip a trailing
  `/1`/`/2`) as bwa does.

### T-007: Paired-end not implemented (pairing FLAG bits, RNEXT/PNEXT/TLEN)
- **Status:** OPEN
- **Severity:** high
- **Affected field(s):** FLAG (0x1/0x2/0x40/0x80, mate-reverse 0x20), RNEXT, PNEXT, TLEN.
- **Symptom:** Mates are aligned independently (`align_single` per read,
  `src/main.rs`); bwa-rs emits FLAG `144` where bwa-mem2 emits `147` (missing
  0x1 paired / 0x2 proper-pair), and RNEXT/PNEXT/TLEN are always `*`/0/0.
- **Suspected cause:** No pairing/rescue/insert-size logic; the CLI loop never
  considers the mate.
- **Resolution:** Implement proper-pair detection, mate-coordinate fields, TLEN,
  and the full pairing FLAG bits.

### T-008: RNAME hardcoded to `"ref"`; no multi-contig coordinate mapping
- **Status:** OPEN
- **Severity:** high
- **Affected field(s):** RNAME, POS (multi-contig).
- **Symptom:** Every mapped record reports RNAME `ref` instead of the contig name
  (`1`); the @SQ header already uses the real name, so output is internally
  inconsistent and invalid for multi-contig references.
- **Suspected cause:** `src/main.rs write_sam_record` hardcodes `"ref"`; the
  aligner returns a single global offset with no contig (`bns`-style) mapping.
- **Resolution:** Add a contig-offset table; map a global forward position to
  (contig, offset) and emit the contig name + per-contig POS.

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
