# Ticket Series: bwa-rs vs bwa-mem2 Output Divergence

Tracking the investigation into why `bwa-rs` SAM output is not byte-identical to
upstream **bwa-mem2** on the bwa-mem2 performance dataset (D2 / SRR7733443
against `human_g1k_v37.fasta`).

**Status:** Phase 1 (format/field/seeding, T-001 … T-024) and Phase 2 (placement
fidelity, T-025 … T-027) **complete** (2026-06-04). bwa-rs now has: full
human-genome index (64-bit SA), both strands, soft-clipping, min-score filter,
QNAME trim, paired-end fields, multi-contig coordinates, SEQ/QUAL, `M`-form CIGAR,
the optional-tag set (NM/MD/MC/AS/XS), index persistence, parallel alignment, the
`mem_approx_mapq_se` + `mem_sam_pe` MAPQ formulas, robust `mem_pestat` insert-size
estimation, mate rescue, faithful `ksw_extend` region extension, max-extended-score
primary selection, and `mem_chain_flt` chain filtering.

**Where parity stands (chr1, vs the bwa-mem2 2.2.1 binary, after T-027):**
- **Uniquely-mapping reads (300-uniq) — effectively at parity.** Per-field vs
  bwa-mem2: POS/RNAME/RNEXT/PNEXT/TLEN **0%**; FLAG **0.7%**; CIGAR/MAPQ/SEQ/QUAL
  **0.3%** (a handful of mate-rescue / soft-clip-boundary records). SAM columns
  1–10 are byte-identical to the T-026 baseline.
- **Multi-mapper / repetitive sample (sub2k) — still diverges.** POS **30.9%**,
  CIGAR **26.0%**, MAPQ **39.5%**, FLAG **33.8%**, TLEN **43.9%**, SEQ/QUAL
  **14.4%**, RNAME/RNEXT 0.6% — down from POS 40.7% / CIGAR 34.1% at the start of
  Phase 2.

**The remaining frontier — Phase 3 (T-028 … below).** The Phase-1/2 residual notes
converge on one architectural gap: **bwa-rs collapses each read to a single
placement before pairing.** `align_read` returns one max-extended-score region and
ignores its mate; the production `align_batch` path aligns every read independently;
the mate enters only *after* placement (`rescue_mate` when one end is unmapped, then
`pair_mapq` / `mate_fields`). So bwa-rs cannot (a) choose the mate-concordant locus
among equally-scoring multi-mapper placements, nor (b) surface the inexact
second-best locus that grades MAPQ down. bwa keeps the full per-read candidate-region
array (`mem_alnreg_v`) and resolves placement **jointly across the pair** in
`mem_sam_pe` / `mem_pair`. Phase 3 opens with an **investigation ticket (T-028)**
that buckets the residual `sub2k` records by root cause, then the implementation
tickets it prioritizes. Verification: `sub2k` is the Phase-3 target; the `uniq` set
stays the zero-regression gate (it is at parity, so any placement change there is a
regression).

**Verification baseline:** chr1 (249 Mbp), `-k 19 -t 16`. Two sets from SRR7733443
with normalized QNAMEs: **300 uniquely-mapping pairs** (`uniq`, the zero-regression
gate) and a **2000-pair multi-mapper sample** (`sub2k`, the Phase-2/3 target).
Compared with `bench/compare.sh` against the bwa-mem2 binary outputs
(`bench/bwamem2_uniq.sam`, `bench/bwamem2_2k.sam`).

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
- **Status:** FIXED (branch `t001-iupac-encode`)
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
- **Resolution (implemented):** `reference::encode_base` now maps any base outside
  `A/C/G/T/N` to `N` (`0b100`) instead of returning `Err`, matching
  `utils::encode_sequence`. The unsanitized `human_g1k_v37.fasta` (and any
  IUPAC-containing FASTA) now indexes without error; the encoding is byte-identical
  to the prior sanitization workaround (M/R → N).

### T-002: Cannot index references larger than ~2.147 Gbp (32-bit suffix array)
- **Status:** FIXED (PR `fix-oom-and-64bit`) — see Resolution.
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

## Phase 1: format, fields, and seeding (ALL FIXED)

Ordered roughly by impact. Verification set: chr1, 300 unique read pairs. Every
ticket in this section (T-005 … T-024) is resolved; the residuals each names are
placement-fidelity diffs now consolidated into Phase 2 (T-025 … T-027).

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
- **Status:** FIXED (branch `t009-t010-seq-qual-cigar`, bundled with T-010)
- **Severity:** high
- **Affected field(s):** SEQ, QUAL.
- **Symptom:** Every record had SEQ `*` and QUAL `*`; bwa-mem2 emits the read
  sequence (reverse-complemented on the 0x10 strand) and qualities.
- **Cause:** `write_sam_record` (`src/main.rs`) and `write_paired_record`
  (`src/sam.rs`) hardcoded `"*"`, and the buffered paired loop dropped the read
  bases/qualities before the emit pass.
- **Resolution:** Added `oriented_seq_qual(bases, qual, reverse)` in `src/sam.rs`,
  which decodes the 2-bit read bases to SEQ and orients to the reference strand
  (reverse-complement SEQ + reversed QUAL when `reverse_strand`); empty qual → `*`.
  Threaded the read bases + raw quality string through both writers; the paired
  loop now buffers a `ReadAln { qname, bases, qual, result }` so the emit pass has
  the sequence. **Measured (chr1 / 300-uniq vs bwa-mem2): SEQ 100% → 0.5% (3/600),
  QUAL 100% → 0.5%.** The residual 3 are T-019 mate-rescue reads bwa-mem2 places on
  the opposite strand (the A/B sequences are reverse complements) — an alignment
  difference, not a SEQ bug.

### T-010: CIGAR uses `=`/`X` instead of `M`
- **Status:** FIXED (branch `t009-t010-seq-qual-cigar`, bundled with T-009)
- **Severity:** medium
- **Affected field(s):** CIGAR.
- **Symptom:** bwa-rs emitted `151=` / `43=7X…`; bwa-mem2 collapses matches and
  mismatches to `M` (`151M`).
- **Cause:** The aligner builds `Eq`/`X` ops and the writers rendered them via
  `Display` (which keeps `=`/`X`).
- **Resolution:** Added `Cigar::to_sam_string()` (`src/types.rs`) — renders `=`,
  `X`, and `M` all as `M` with adjacent runs merged, other operators unchanged —
  and switched both SAM record writers to it. `Display for Cigar` is unchanged
  (`=`/`X`), so internal callers and tests that rely on operator-level CIGAR are
  untouched. **Measured (chr1 / 300-uniq vs bwa-mem2): CIGAR 99.3% → 5.3% (32/600),
  `=`/`X` CIGARs 594 → 0.** The residual 32 are 1–3 bp boundary shifts owned by
  T-018 (banded SW), not CIGAR formatting.

### T-011: Optional tags missing (NM, MD, AS, XS, …)
- **Status:** FIXED (branch `t011-optional-tags`) — XS deferred to T-015.
- **Severity:** medium
- **Affected field(s):** optional tag columns.
- **Symptom:** CLI records carry no tags; bwa-mem2 emits `NM:i:`, `MD:Z:`,
  `AS:i:`, `XS:i:`, `MC:Z:`, etc. (`AlignmentResult` already computes `nm` and an
  MD string, but `write_sam_record` drops them.)
- **Resolution:** Both SAM record writers now append optional tags on mapped
  records, in bwa-mem2's column order: `NM:i` · `MD:Z` · `MC:Z` (paired only) ·
  `AS:i`. `NM`/`AS`/`MD` reuse the already-computed `AlignmentResult.nm`/`.score`/
  `.md_tag` (previously dropped at write time). `MC:Z` (mate CIGAR, M-form) is
  carried on a new `MateFields.mc: Option<String>` populated in
  `paired::mate_fields` from the mate's `Cigar::to_sam_string()` (present only when
  the mate is mapped). `write_sam_record` (`src/main.rs`, single-end) emits
  NM/MD/AS (no mate → no MC); `write_paired_record` (`src/sam.rs`) emits all four.
  Unmapped records emit no tags. **XS is intentionally not emitted:** bwa-rs has no
  true second-best *alignment* score (only the MEM-score MAPQ heuristic), so a
  byte-matching XS depends on the region-scoring/MAPQ port — deferred to **T-015**.
  New tag-order/absence tests in `src/main.rs`, `src/sam.rs`, `src/paired.rs`; all
  tests pass, `cargo clippy` (lib+bin) clean.
- **Measured (chr1 / 300-uniq vs bwa-mem2):** NM matches 592/594, MC 586/588, AS
  505/594 (the AS gap is scoring-defaults divergence → **T-014**), all residuals
  tracing to the known T-019 mate-rescue reads. The comparison also exposed that
  the already-computed MD string was malformed and reverse-strand-incorrect
  (186/594 differing) — a pre-existing computation defect, not a tag-emission bug —
  fixed in follow-up **T-020** (MD differing 186 → 2).

### T-012: `index` subcommand is a no-op; `mem` rebuilds the index every run
- **Status:** FIXED (branch `t012-index-persist`)
- **Severity:** medium (usability/perf, not a SAM-content diff)
- **Symptom:** `bwa-mem index -p PREFIX` built the FM-index, printed
  `Indexed N bases`, and wrote **no files**; `mem` rebuilt the whole 2N index in
  memory on every invocation (`src/main.rs`). `FMIndex::save/load` exist and are
  64-bit, but the CLI didn't call them.
- **Resolution:** `index` now persists the 2N index via `FMIndex::save` to
  `<prefix>.bwarsidx` (single file; `-p` prefix defaults to the reference path,
  matching bwa-mem2's reference-keyed index convention). `mem` gained a matching
  optional `-p` and now **loads** the index instead of rebuilding: if
  `<prefix>.bwarsidx` exists it calls `FMIndex::load` and verifies
  `index.n_fwd == reference.total_len()` (a present-but-mismatched/stale index
  errors rather than silently rebuilding); a corrupt index propagates the load
  error; an **absent** index falls back to the in-memory `build_2n` with a stderr
  note (non-breaking). `mem` still reads the FASTA for the `Reference` (header
  `@SQ`, `locate`, and the `as_slice_2n` ref data) — only the expensive suffix-array
  build is skipped. A pure `index_path_for(prefix)` helper (unit-tested) derives the
  `.bwarsidx` path. Per the ticket's "unify or document": `src/mmap_index.rs` now
  documents the two on-disk formats (64-bit `.bwarsidx` used by the CLI vs the
  experimental 32-bit zero-copy mmap format, ≤4.29 Gbp, unused by the CLI) — left
  separate, not unified.
- **Save/load performance fix:** `FMIndex::save`/`load` wrote/read the suffix array
  through an **unbuffered** `File` (`SuffixArray::write_to`/`read_from` do one 8-byte
  syscall per entry — ~498 M syscalls for chr1's 2N SA), making a load *slower* than
  an in-memory rebuild. Wrapped both in `BufWriter`/`BufReader` (`src/fm_index.rs`).
  On chr1 (249 Mbp / 498 Mbp 2N, 4.49 GB index): index build+save **856 s → 132 s**,
  `mem` load **409 s → 44 s** (vs the prior ~106 s in-memory rebuild per run — load
  now wins and amortizes across runs).
- **Measured (chr1 / 300-uniq vs bwa-mem2, `bench/compare_t012.md`):** the load-path
  SAM is **byte-identical** to the T-013 in-memory-build baseline
  (`bench/bwars_2n_uniq_t013.sam`) except the `@PG CL:` output-filename field, and the
  per-field divergence vs bwa-mem2 is unchanged — POS 0, RNAME 0, RNEXT 0, PNEXT 0;
  FLAG 4, CIGAR 4, TLEN 4, MAPQ 2, SEQ 2, QUAL 2. **Zero regressions.** CLI smoke tests
  confirm all four paths (save, load, stale `n_fwd` guard → error, missing-index
  fallback build) and load==build byte-identity on a tiny reference. All unit tests
  pass; `cargo clippy -- -D warnings` and `cargo fmt --check` clean.

### T-013: CLI alignment is single-threaded
- **Status:** FIXED (branch `t013-parallel-align`)
- **Severity:** medium (perf)
- **Symptom:** `-t` only sized the rayon pool used for SA construction; the read
  loop in `run_mem` (`src/main.rs`) aligned one read at a time via `align_single`.
  bwa-mem2 aligned 200k pairs in ~9 s on 16 threads.
- **Resolution:** `run_mem` now drives the existing `ParallelAligner::align_batch`
  (rayon `par_iter`) instead of the per-read loop. The paired path reads the whole
  file into owned `RawRead`s, flattens to one ordered query list (read1 then read2
  per pair), aligns the entire batch in parallel, then walks the results back in
  order to rebuild the two-phase `buf` and accumulate the `InsertSizeDistribution`
  (the whole file must be buffered here — proper-pair/`pair_mapq` need the complete
  distribution, so it cannot be chunked). The single-end path aligns in
  `BATCH_SIZE = 65536` chunks to bound memory. `par_iter().collect()` preserves
  input order, so emit order — and thus output — is unchanged. A small `align_all`
  helper unwraps `Vec<Result<…>>` → `Result<Vec<…>>` for both paths. `align_read`
  is pure per-read (no mate context passed, matching the prior `align_single`), so
  results are identical; `align_batch_paired`/`align_batch_with_mates` were *not*
  used because they pass mate context to `align_read` and would change output.
- **Measured (chr1 / 300-uniq vs bwa-mem2, `bench/compare_t013.md`):** output is
  **byte-identical to the T-022 baseline** (only the `@PG CL:` output-filename
  field differs) — per-field divergence unchanged: POS 0, RNAME 0, RNEXT 0, PNEXT 0;
  FLAG 4, CIGAR 4, TLEN 4, MAPQ 2, SEQ 2, QUAL 2. **Zero regressions.** All 324 unit
  tests pass; `cargo clippy -- -D warnings` and `cargo fmt --check` clean.
- **Note:** Wall time on this small set is dominated by the in-memory index rebuild
  (T-012), not alignment; large-set throughput measurement remains gated on the
  repetitive-read seeding blow-up (**T-023**).

### T-014: Scoring/seed defaults differ from bwa-mem2
- **Status:** FIXED (branch `t016-headers-t014-defaults`)
- **Severity:** medium (changes which alignment wins → POS/CIGAR/MAPQ)
- **Symptom:** Believed `gap_open` 6 vs 5 and default `min_seed_len` 10 vs 19.
- **Investigation:** Verified against the actual benchmark binary
  (`bench/tools/bwa-mem2-2.2.1_x64-linux/bwa-mem2 mem`, run with default options):
  its defaults are `-A 1 -B 4 -O [6,6] -E [1,1] -L [5,5] -k 19 -T 30`. The
  `gap_open 6 vs 5` claim in this file and `docs/PERFORMANCE.md` was WRONG —
  bwa-mem2 uses `gap_open = 6`, which bwa-rs's `Scoring::default()` already matched.
  The library `DEFAULT_MIN_SEED_LEN` was also already 19. The only real divergence
  was the CLI `-k` flag, whose default was still 10 (so the binary needed `-k 19`
  to match).
- **Resolution:** Changed the CLI `-k` default 10 → 19 (`src/main.rs`) so the
  binary's default seeding matches bwa-mem2 without an explicit flag. Corrected the
  false `gap_open 6 vs 5` row and the `min_seed_len` row in `docs/PERFORMANCE.md`.
  All other scoring defaults (match 1, mismatch 4, gap_open 6, gap_extend 1,
  clip 5, min_score 30) already matched bwa-mem2 and are unchanged.

### T-015: MAPQ computation differs
- **Status:** FIXED (branch `t015-mapq`)
- **Severity:** medium
- **Affected field(s):** MAPQ.
- **Symptom:** `calculate_mapq` (`src/alignment.rs`) used a second-best-seed ratio
  heuristic, not bwa-mem2's `mem_approx_mapq_se` formula (based on sub-optimal
  score, seed coverage, etc.). Because it only ever saw the single best chain its
  effective output was MAPQ 60 for every mapped read; the rescue path used a
  separate `rescue_mapq` linear placeholder. Matched on the easy MAPQ-60 unique set
  but diverged on rescued mates and multi-mappers.
- **Resolution:** Ported bwa's `mem_approx_mapq_se` as a pure `approx_mapq_se`
  (`src/alignment.rs`) — the bwa-mem2 default `mapQ_coef_len = 50` branch, with
  `csub`/`frac_rep` not modeled (treated as 0) — plus an `aligned_span` helper
  (`l = max(query_core, ref_core)`). `align_read` now derives the suboptimal
  alignment score `sub` (and `sub_n`) from the best-scoring secondary chain whose
  placement does not overlap the primary (building up to `MAX_SUB_CHAINS = 32`
  alternatives), then sets MAPQ via the formula; primary selection is unchanged, so
  POS/CIGAR/FLAG are byte-identical. The dead `calculate_mapq` and the `rescue_mapq`
  placeholder were removed; `rescue_mate` now uses the same formula. New unit tests
  cover every formula branch (`approx_mapq_se_*`, `aligned_span_*`).
  **Measured (chr1 / 300-uniq vs bwa-mem2, `bench/compare_t015{,_summary}.md`):** no
  placement regressions (field divergence identical to T-019). The two mate-rescue
  reads moved from the flat placeholder 12 to formula values tracking bwa-mem2
  (`100124` 12→27 vs 30; `100170` 12→33 vs 59); the residual gap is the differing
  CIGAR core (T-014/T-018), not the formula. On reads with **identical placement** in
  the multi-mapper sub2k set (1376 reads) MAPQ exact-match vs bwa-mem2 went 72.5% →
  72.7% and mean |MAPQ−bwa| 14.65 → 13.09.
- **Residual (other root cause — follow-up):** 370/375 of the remaining
  identically-placed MAPQ diffs are bwa-rs scoring *higher* (e.g. 60 vs 0) on
  heavily soft-clipped reads whose aligned fragment is repetitive: bwa-mem2 finds the
  secondary copy and grades MAPQ down, but bwa-rs's seeding surfaces no non-overlapping
  secondary chain, so `sub = 0` → baseline → 60. The bottleneck is now
  **secondary-hit detection in seeding** (the same machinery gating the deferred XS
  tag, T-011), not the MAPQ formula. Paired-end MAPQ recalculation (`mem_sam_pe`) and
  the repetitive-read seeding time/memory blow-up (pre-existing; T-013) are also
  follow-ups.

### T-016: Header lines differ (@HD, @PG)
- **Status:** FIXED (branch `t016-headers-t014-defaults`)
- **Severity:** low
- **Affected field(s):** @HD, @PG (and @SQ ordering for multi-contig).
- **Symptom:** bwa-rs emitted `@HD VN:1.0 SO:unsorted` and
  `@PG ID:bwa-rs PN:bwa-rs VN:0.1.0`; bwa-mem2 emits no @HD and a
  `@PG ID:bwa-mem2 PN:bwa-mem2 VN:2.2.1 CL:<command line>` line. (@SQ already matches.)
- **Resolution:** `write_header` (`src/main.rs`) now drops the `@HD` line entirely
  (bwa-mem2 emits none) and emits a bwa-style `@PG` line carrying the crate version
  (`env!("CARGO_PKG_VERSION")`) and a `CL:` field with the actual command line
  (`std::env::args()`), replacing the hardcoded `VN:0.1.0` and absent CL. Header
  construction was factored into a pure `header_text(reference, cmdline)` helper with
  a unit test (`header_omits_hd_and_emits_pg_with_cmdline`).
- **Intentional residual divergence:** the `@PG` `ID`/`PN` stay `bwa-rs` rather than
  impersonating `bwa-mem2`, so the `@PG` line is *not* byte-identical by design — an
  honest program identity is preferred over byte-matching the header. Everything else
  in the header (no `@HD`, `@SQ`, bwa-style `@PG` shape with `VN`/`CL`) now matches.

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
- **Status:** FIXED (branch `t018-banded-sw`)
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
- **Resolution:** Stopped pinning the whole seed as forced `Eq`. `build_alignment`
  now anchors only at the seed's *start* corner `(query_start, ref_start)` and folds
  the seed region into the forward affine DP (`affine_extend_forward` over
  `query[seed.query_start..]` from `seed.ref_start`), mirroring bwa's `ksw_extend`
  model — the seed interior is now free to absorb a mismatch/gap so the banded DP can
  shift the placement within the band instead of locking it to the seed. The forced
  middle `Eq(seed.length)` push is removed; `ref_start` derives from the backward
  extension's leftmost consumed base. Backward extension is unchanged. All 238 lib
  tests pass plus two new ones (`test_exact_match_regression`,
  `test_seed_interior_free_full_length_alignment`).
  **Measured (chr1 / 300-uniq vs bwa-mem2, `bench/compare_t018.md`): POS 17 → 2
  (0.3%), CIGAR 32 → 4 (0.7%), PNEXT 17 → 2.** No regressions; FLAG (8), TLEN (8),
  MAPQ (2), SEQ/QUAL (3) are unchanged. The 2 residual POS diffs (reads
  `100170`/`100124`) are **T-019** mate-rescue reads (bwa-rs leaves a mate unmapped,
  `CIGAR *` / MAPQ 0), not seed anchoring. Of the 4 residual CIGAR diffs, 1 is that
  same T-019 case and ~3 are ≤2 bp soft-clip-boundary shifts attributable to scoring
  tie-breaks (**T-014**, `gap_open` 6 vs 5) and the not-yet-implemented joint
  optimization across the seed corner — tracked as follow-ups along with MAPQ
  (**T-015**).

### T-019: Mate rescue not implemented (unmapped mate via banded Smith-Waterman)
- **Status:** FIXED (branch `t019-mate-rescue`)
- **Severity:** medium
- **Affected field(s):** FLAG (0x4/0x8/0x2), POS, CIGAR — mapped-vs-unmapped status of a mate.
- **Symptom:** When one mate maps and the other does not, bwa-mem2 runs a banded
  Smith-Waterman around the mapped mate's expected position (`mem_matesw`) and
  often rescues the unmapped mate, marking the pair proper. bwa-rs left it
  unmapped. On the chr1 / 300-uniq set this was the entire 8/600 residual FLAG
  divergence after T-007.
- **Cause:** No `mem_matesw` equivalent; `align_read` only seeds from the read's own
  MEMs and never aligned a mate into the partner's insert-size window.
- **Resolution:** Added `local_align` (affine-gap **local** Smith-Waterman over a
  reference window with free reference start/end and soft-clipped query ends) and
  `Aligner::rescue_mate` (`src/alignment.rs`). When exactly one mate maps,
  `rescue_mate` derives the FR insert-size window in forward coordinates from the
  `InsertSizeDistribution` (T-007) around the mapped mate — downstream + reverse
  strand if the mate is forward, upstream + forward strand if the mate is reverse —
  orients the orphan onto the opposite strand, runs `local_align`, and returns a
  forward-coordinate alignment gated by `min_score` (`-T`). `ParallelAligner`
  delegates; the paired emit loop in `src/main.rs` invokes rescue before computing
  proper-pair/mate fields. The dead, semantically-wrong `rescue_orphan` stub (it set
  the 0x8 mate-unmapped bit on the orphan) was removed from `src/paired.rs`.
  **Measured (chr1 / 300-uniq vs bwa-mem2, `bench/compare_t019.md`), T-020 → T-019:
  POS 2 → 0, PNEXT 2 → 0, FLAG 8 → 4, TLEN 8 → 4, SEQ 3 → 2, QUAL 3 → 2.** The two
  rescuable mates (`100124`, `100170`) are now placed at the exact bwa-mem2
  coordinates with matching proper-pair FLAGs (`147`/`163`); no regressions (only
  those two reads changed behaviour, unmapped → correctly placed).
- **Residuals (other tickets, not rescue defects):** the rescued reads' MAPQ
  (`12` vs bwa `30`/`59`) is **T-015** (MAPQ formula); their CIGAR-core exactness
  (`44S37M1I2M1D58M9S` vs `44S98M9S`, `19S75M57S` vs `19S88M44S` — `local_align`
  gap/clip tie-breaks vs bwa's `ksw_extend` end-bonus model) traces to T-014/T-015;
  the 4 remaining FLAG diffs (`100198`/`100277`) are reads bwa-mem2 *also* leaves
  unmapped — an unmapped-mate strand-flag (0x10/0x20) convention follow-up; the one
  1 bp CIGAR boundary shift (`100142`) is **T-014**/**T-018**.

### T-020: MD:Z string is malformed and computed on the wrong strand
- **Status:** FIXED (branch `t020-md-correctness`, stacked on `t011-optional-tags`)
- **Severity:** high (byte-affecting on ~31% of records once T-011 emits the tag)
- **Affected field(s):** MD.
- **Symptom:** Surfaced by the T-011 chr1 comparison: MD differed from bwa-mem2 on
  **186/594** records. Two independent defects in `AlignmentResult::mdz_string`
  (`src/types.rs`):
  1. **Missing `0` separators.** The SAM MD grammar
     (`[0-9]+(([A-Z]|\^[A-Z]+)[0-9]+)*`) requires a number before *every* mismatch
     base and at the string's end, even when zero. bwa-rs only emitted the count
     when `> 0`, so adjacent mismatches produced `…4AA18` instead of `…4A0A18`, a
     leading mismatch produced `C119…` instead of `0C119…`, and a mismatch after a
     deletion lost its separator.
  2. **Wrong strand.** MD was computed at the seed-extension call site *before* the
     2N→forward remap, so a reverse-strand read recorded reverse-complement-half
     reference bases in mirrored order — e.g. bwa-rs `132A18` where bwa-mem2 emits
     `18T132` (`A` = complement of `T`, position mirrored).
- **Cause:** `mdz_string` flushed `match_count` only when `> 0`; and
  `build_alignment` (`src/alignment.rs`) called `mdz_string(query, ref_seq)` in the
  2N indexed space rather than against the forward reference.
- **Resolution:** (1) `mdz_string` now always flushes the running match count before
  each mismatch/deletion and once at the end, yielding spec-compliant `0` separators
  and trailing number. (2) `build_alignment` computes MD *after* the strand remap;
  for a reverse-strand read it walks the forward-oriented read
  (`reverse_complement(query)`) along the forward CIGAR against the forward
  reference. NM/AS are mismatch/indel *counts* (strand-invariant) and are unchanged.
  Two existing unit tests that encoded the malformed format were corrected
  (`2T2^C` → `2T2^C0`, `TGC` → `0T0G0C0`); added `test_mdz_leading_mismatch`,
  `test_mdz_adjacent_mismatches`, and `test_mdz_reverse_strand_uses_forward_reference`.
  **Measured (chr1 / 300-uniq vs bwa-mem2): MD differing 186 → 2 (0.3%).** The 2
  residuals (`100142`, `100189`) are not MD defects — they are 1 bp soft-clip-boundary
  shifts (`139M12S` vs `138M13S`) owned by **T-014**/**T-018**; the MD correctly
  reflects bwa-rs's own CIGAR. Core fields and NM/MC/AS unchanged; no regression.

### T-021: Suboptimal-hit detection in seeding (drives MAPQ `sub` and XS)
- **Status:** FIXED (machinery + `frac_rep` + XS; bulk MAPQ residual reassigned to T-024) — branch `t021-seed-secondary`
- **Severity:** medium
- **Affected field(s):** MAPQ, XS (optional tag).
- **Symptom:** Surfaced by T-015. On reads where the primary placement matches
  bwa-mem2, MAPQ still diverges on 375/1376 identically-placed sub2k reads — and
  **370 of those are bwa-rs scoring *higher*** (e.g. MAPQ 60 vs bwa 0 on `95M56S`,
  `121S30M`, `45M106S`). These are heavily soft-clipped reads whose aligned fragment
  is repetitive: bwa-mem2 finds the secondary copy and grades MAPQ down, but bwa-rs's
  seeding surfaces no non-overlapping secondary chain, so the T-015 `sub` term is 0 →
  baseline → MAPQ 60. The same missing second-best *alignment* score is why XS is not
  emitted (deferred from **T-011**).
- **Suspected cause:** `find_mems` + `filter_mems` + `chain_seeds`
  (`src/seed.rs`, `src/chaining.rs`) produce fewer/different secondary chains than
  bwa's seeding (re-seeding of long MEMs, seed dropping, chain de-overlap), so a
  repetitive fragment that bwa maps to ≥2 loci yields a single bwa-rs chain. T-015's
  `align_read` then computes `sub` over an incomplete candidate set.
- **Resolution:** Ported four pieces of bwa's secondary-region machinery:
  (1) **re-seeding** — a rare long SMEM (occ ≤ `SPLIT_WIDTH = 10`, len ≥
  `SPLIT_FACTOR×min_seed_len`) is split at its midpoint and re-seeded for a *more
  frequent* sub-seed (`longest_match_min_occ`, bwa `min_intv = parent_occ + 1`) so a
  repeated core surfaces the repeat's loci (`src/mem_finder.rs`);
  (2) an **occurrence cap** `DEFAULT_MAX_OCC = 500` (`FMIndex::find_all_capped`) —
  over-frequent SMEMs are not enumerated for chaining;
  (3) **faithful `sub`/`sub_n`** — `align_read` marks secondaries `mem_mark_primary_se`-
  style (a distinct-locus candidate whose query span overlaps the primary by more than
  `MASK_LEVEL = 0.5` contributes its score to `sub`), replacing the ref-space
  non-overlap heuristic (`src/alignment.rs`);
  (4) **`frac_rep`** — `mapq *= (1 − frac_rep)`, the merged query fraction covered by
  over-frequent (occ > max_occ) SMEMs (`repetitive_fraction`, threaded into
  `approx_mapq_se`).
  `XS:i` is now emitted on every mapped record in both writers (after `AS:i`),
  unconditionally as bwa-mem2 does — mechanically closing the **T-011** XS deferral.
- **Measured (chr1 / 300-uniq, `bench/compare_t021{,_summary}.md`):** **no regressions** —
  field divergence byte-identical to T-019/T-020 (POS 0, FLAG 4, MAPQ 2, CIGAR 4, …);
  594 unique reads stay MAPQ 60; `frac_rep` down-weights 0 of them; XS now present.
  On repetitive reads the machinery fires (XS non-zero; MAPQ drops to 0/22 where a
  secondary's score nears the primary).
- **Residuals (other root cause — follow-up):** the bulk of T-015's over-scoring reads
  are **not** closed. Their secondary locus is *inexact* (bwa DP-extends a short shared
  seed across mismatches/indels, e.g. read `SRR7733443.26`'s 72 bp core maps elsewhere
  as `96M55S`), which bwa-rs's exact-MEM seeding never surfaces → `sub = 0` → MAPQ 60;
  and they are only *moderately* repetitive (occ ≤ 500) so `frac_rep = 0` doesn't apply.
  Faithful approximate-secondary seeding is split to **T-024**. Full multi-mapper
  measurement remains **OOM-blocked (T-023)**. XS *values* are not byte-identical
  (partial secondary set) — same dependency on T-024.

### T-022: Paired-end MAPQ recalculation (`mem_sam_pe`) not implemented
- **Status:** FIXED (branch `t022-pe-mapq`) — exact rescued-read values gated on T-024/T-018.
- **Severity:** low
- **Affected field(s):** MAPQ (paired / rescued reads).
- **Symptom:** After T-015 the two chr1 mate-rescue reads report MAPQ `27`/`33` where
  bwa-mem2 reports `30`/`59`. bwa derives the single-end MAPQ with
  `mem_approx_mapq_se` (now ported) and then *recomputes* a paired MAPQ in
  `mem_sam_pe` from the pairing score `o`/`subo` (combining both mates and the
  insert-size likelihood). bwa-rs applied only the SE formula per mate.
- **Cause:** No PE pairing-score MAPQ pass; `src/main.rs` emitted each mate's
  SE MAPQ unchanged. Part of the residual is also the differing CIGAR core
  (T-014/T-018), which lowers the SE score feeding the formula.
- **Resolution:** Ported `mem_sam_pe`'s paired-MAPQ recalculation as pure helpers in
  `src/paired.rs`: `raw_mapq` (bwa `raw_mapq`), a from-scratch `erfc` (Abramowitz &
  Stegun 7.1.26, no new deps), `insert_bonus` (bwa `mem_pair`'s insert-size
  log-likelihood term, `.721·ln(2·erfc(|ns|/√2))·a`), and `pair_mapq` which assembles
  the pairing score `o = score1 + score2 + insert_bonus`, `subo = max(0, score1 +
  score2 − PEN_UNPAIRED)` (PEN_UNPAIRED = 17), `q_pe = raw_mapq(o−subo, a)` clamped and
  down-weighted by `(1 − ½(frac_rep1+frac_rep2))`, then caps each mate's SE MAPQ via
  `q_se > q_pe ? q_se : min(q_pe, q_se+40)`. `frac_rep` is now persisted on
  `AlignmentResult` (mirrors `xs`; set in `align_read`, 0 for rescued/unmapped); added
  `match_score()` accessors on `Aligner`/`ParallelAligner`. The paired emit loop in
  `src/main.rs` recomputes both mates' MAPQ via `pair_mapq` when the pair is proper.
- **Measured (chr1 / 300-uniq vs bwa-mem2, `bench/compare_t022{,_summary}.md`):** **no
  regressions** — per-field divergence byte-identical to T-021 (POS 0, FLAG 4, MAPQ 2,
  CIGAR 4, …); the 592 identically-placed reads stay **100% MAPQ-exact**. The two
  mate-rescue reads pick up pairing evidence: `100170` 33 → **60** (bwa 59, off-by-1),
  `100124` 27 → **60** (bwa 30, overshoots).
- **Residual (other root cause — follow-up):** `100124` overshoots because bwa's `q_pe`
  is low there (`o−subo ≈ 5`) — bwa finds a *competing* pairing at the rescued read's
  repetitive secondary locus, raising `subo`. bwa-rs keeps a single placement per end
  (no `mem_alnreg_v` candidate array), so `subo = score_un` and `q_pe ≈ 60` for any
  near-mean concordant pair. Reconstructing the competing pairing needs the
  candidate-region pairing score — **T-024** (which also feeds XS). The rescued-read
  CIGAR core (`44S37M1I2M1D58M9S` vs bwa `44S98M9S`) is **T-018/T-014** and lowers the
  score feeding the formula.

### T-023: Repetitive-read seeding time/memory blow-up
- **Status:** FIXED (branch `t023-rescue-window-cap`)
- **Severity:** medium
- **Affected field(s):** N/A — throughput/RSS, not SAM content.
- **Symptom:** A single full `sub2k` (2000-pair) alignment against chr1 peaks at
  ~68 GB RSS and crawls on repetitive reads; the T-015 sub2k measurement could only
  complete a partial set (two of eight parallel shards were OOM-killed). Per-read wall
  time on repeat-heavy reads is seconds, not microseconds.
- **Suspected cause:** Pre-existing — `find_mems` (`src/seed.rs`) materialises the full
  SA occurrence list for high-frequency k-mers (a repeat seed can occur at thousands
  of loci), with no occurrence cap / `max_occ` cutoff like bwa's. T-015's per-read
  secondary-chain builds (`MAX_SUB_CHAINS = 32`) add to per-read cost but are not the
  OOM root.
- **Resolution:** Cap seed occurrences (bwa `-c` / `max_occ`) and avoid collecting all
  positions for ultra-frequent seeds; combine with batched/parallel alignment
  (**T-013**). Independent of the byte-identity goal but blocks large-set verification.
- **Update (T-021):** An occurrence cap `DEFAULT_MAX_OCC = 500`
  (`FMIndex::find_all_capped`, used by `find_supermaximal_mems`) is now in place —
  over-frequent SMEMs are no longer enumerated. This is a *partial* mitigation: a
  positional `sub2k` slice still OOM-killed mid-run during T-021 verification (the
  reseed/secondary-chain builds and per-read DP on repeat-heavy reads remain
  uncapped), so large-set measurement is still blocked. Batched/parallel alignment
  (**T-013**) and a per-read work cap remain to do.
- **Investigation (T-023 fix):** The suspected seeding cause was *refuted* as the
  residual blocker. With T-021's `max_occ` cap and T-013's batched alignment, a full
  `sub2k` run reaches the **emit phase**: the master binary writes 913 of 4014 records
  before being OOM-killed. Records are written only in the single-threaded step-2 emit
  loop (`src/main.rs`), so primary alignment — seeding + extension — had already
  completed for all 2000 pairs. The kill is in **mate rescue**: `rescue_mate` →
  `local_align` allocates a dense `3 × (orphan_len+1) × (window+1)` i32 DP matrix, and
  the window width is `dist.upper_bound() = mean + 4·std_dev`. `InsertSizeDistribution`
  estimated mean/std with a **naive running accumulator over all FR-concordant inserts,
  with no outlier filtering**, so a few concordant-by-orientation but mismapped pairs
  (mates megabases apart) inflate `std_dev` to millions → `upper_bound` → tens of
  millions → window approaches the whole chromosome (≥55 M bp) → a single rescue
  allocates >100 GB → OOM. This is the T-003 dense-matrix bug class in the rescue path.
  bwa-mem2 avoids it via robust `mem_pestat` estimation (its log on this data:
  `# candidate unique pairs (FF,FR,RF,RR): (0,859,0,1)`; `boundaries for computing mean
  and std.dev: (1, 689)`; `mean and std.dev: (342.42, 112.45)`).
- **Resolution:** Ported bwa's `mem_pestat` robust insert-size estimation
  (`InsertSizeDistribution::from_insert_sizes`, `src/paired.rs`), replacing the naive
  streaming `add`: sort the FR inserts, take p25/p75 (bwa's `(q·n + .499)` index), and
  **filter outliers to `[p25 − 2·IQR, p75 + 2·IQR]` before computing mean/std**;
  proper-pair bounds are `[p25 − 3·IQR, p75 + 3·IQR]` widened to include `mean ± 4σ`
  (`OUTLIER_BOUND = 2.0`, `MAPPING_BOUND = 3.0`, `MAX_STDDEV = 4.0`, verified against two
  bwa-mem2 log blocks). `src/main.rs` collects the FR inserts into a vec (the pairs are
  already buffered) and builds the distribution once before the emit loop. No hard cap is
  added — the bounded rescue window falls out of correct estimation, exactly as in
  bwa-mem2; the rescue geometry in `rescue_mate` is unchanged. `with_params` is retained
  (test-only, documented) as a parameter-supplied constructor using the simpler
  `mean ± 4σ` bound.
- **Measured (chr1, `-k 19 -t 16`, 128 GB allocation; `bench/compare_t023_summary.md`):**
  **`sub2k` master OOM-killed (913/4014 records, >128 GB) → T-023 exit 0
  (4000/4000 records, peak RSS 15.21 GB, wall 42 s).** The **300-uniq verification set
  is byte-identical to the T-022 baseline** (excl. `@PG CL:`) — zero regression, since
  unique pairs carry no outliers so robust ≈ naive. sub2k distribution is robust, not
  the degenerate fallback: median |TLEN| 325 (≈ bwa-mem2 mean 342.42), 2702 proper-pair
  flags. All unit tests pass; `cargo clippy -- -D warnings` and `cargo fmt --check` clean.
- **Residual (other tickets):** the `sub2k` per-field divergence vs bwa-mem2 (POS 40.7 %,
  mapped 3594 vs 3987) is the pre-existing placement/seeding gap on a random multi-mapper
  sample — primary POS/CIGAR are independent of the insert distribution. Closing it is
  **T-024** (approximate secondary seeding) / **T-018** (banded SW); T-023 only unblocks
  the measurement.

### T-024: Approximate (inexact) secondary seeding to surface bwa's second-best loci
- **Status:** FIXED (branch `t024-approx-secondary-seeding`)
- **Severity:** medium
- **Affected field(s):** MAPQ, XS (optional tag).
- **Symptom:** Split from T-021. After T-021's re-seeding + `mem_mark_primary_se`-style
  `sub`/`sub_n` + `frac_rep`, the bulk of T-015's over-scoring reads still report MAPQ
  60 where bwa-mem2 reports 0–52 (heavily soft-clipped, moderately repetitive reads).
  Their second-best locus is *inexact* — bwa places it by DP-extending a short shared
  seed across mismatches/indels (e.g. `SRR7733443.26`'s 72 bp core maps at a second
  locus as `96M55S` / `65M1I35M1D26M24S`, not an exact 72 bp MEM). bwa-rs's exact-MEM
  seeding never surfaces that locus, so `sub = 0` → MAPQ 60 and XS values diverge.
- **Suspected cause:** Two gaps in the T-021 `sub`/XS set. (1) **Tandem-repeat copies**
  (e.g. the telomeric `CCCTAA` hexamer) collapse into one chain — `chain_seeds` merges
  seeds within 50 bp and `filter_mems` drops ref-overlapping MEMs — so a read whose core
  maps to N nearby copies yields one placement, `sub = 0`. (2) **Inexact secondary loci**
  are never seeded: rounds 1 (SMEM) + 2 (midpoint reseed) miss a short shared exact core.
- **Resolution:** Ported bwa's **third seeding round** (`bwt_seed_strategy1`,
  `MAX_MEM_INTV = 20`) as `collect_short_seeds` (`src/mem_finder.rs`): scan the read
  left-to-right, collect at each position the shortest exact seed (len ≥ min_seed) whose
  SA interval ≤ 20, advance past it. It is kept **out** of the primary seeding path
  (`find_supermaximal_mems`) so primary placement is unchanged. `align_read`
  (`src/alignment.rs`) now derives `sub`/`sub_n` over the **unfiltered** MEMs (keeping
  tandem copies) **plus** the third-round seeds: each candidate is extended via
  `build_alignment`, duplicates of the primary dropped by bwa's containment rule
  (contained in query AND reference, same strand), a candidate at a distinct reference
  locus overlapping the primary in query space counted as a secondary. Candidates are
  deduped by (query_start, ref_start) before the DP and by (position, strand) after, and
  capped at `MAX_SUB_CANDIDATES = 64`. A guard excludes any candidate scoring *above*
  the chosen primary (bwa's primary is the max-scoring region, so `sub ≤ score`; a
  higher-scoring candidate is a primary-placement gap, **T-018**, not a secondary). XS
  is emitted from the resulting `sub`.
- **Measured (chr1, `-k 19 -t 16`, `bench/compare_t024{,_summary}.md`):** **zero
  regression** — the 300-uniq set is byte-identical to the clean-master baseline on SAM
  columns 1–10, and sub2k placement (cols 1–4,6) is frozen. MAPQ vs bwa-mem2 over the
  sub2k both-mapped records: exact-match **62.3 % → 64.0 %**, mean |MAPQ−bwa| **18.40 →
  16.25**, over-scoring reads (bwa-rs > bwa) **1329 → 1252**, big over-scorers (>+20)
  **1238 → 1128**; identically-placed mean |MAPQ−bwa| 11.90 → 11.14. Throughput/RSS
  unchanged (47 s, 15.2 GB). All unit tests pass; `cargo clippy -- -D warnings` and
  `cargo fmt --check` clean.
- **Residual (other root causes — follow-up):** under-scoring rises (26 → 41) — reads
  where bwa-rs's *primary* is mis-placed (lower AS than bwa-mem2), now exposed by better
  secondary detection. Byte-matching XS *values* and closing the remaining MAPQ gap need
  max-alignment-score primary selection (**T-025**) and faithful banded (`ksw`)
  region-extension scoring at each candidate locus (**T-026**); both change placement.
  Overlaps **T-022** (faithful pairing score).

---

## Phase 2: placement fidelity (T-025, T-026, T-027 FIXED)

The remaining sub2k divergence (POS 40.7 %, CIGAR 34.1 %, with FLAG/TLEN/MAPQ
cascading) is one root cause — final alignment placement on repetitive,
heavily-soft-clipped, multi-mapping reads — broken into the three pieces of bwa's
`mem_align1_core` that bwa-rs still approximates. **T-025** is the smallest and
unblocks the others; **T-026** is the deep core; **T-027** refines the candidate
set. **T-025 + T-026 cut sub2k POS 40.7 % → 31.1 % and CIGAR 34.1 % → 25.9 %** with
zero uniq regression; **T-027** (chain filtering) then corrected the XS
over-counting the faithful placement exposed (sub2k POS 31.1% → 30.9%, CIGAR
≈flat). The residual multi-mapper placement gap is carried forward to **Phase 3**.
Verification: the `sub2k` set (Phase-2 target); the `uniq` set stays the
zero-regression gate (it is already at parity, so any placement change there is a
regression).

### T-025: Select the primary by extended-alignment score, not chain score
- **Status:** FIXED (branch `t025-max-score-primary`)
- **Severity:** high (largest remaining POS/CIGAR driver)
- **Affected field(s):** POS, CIGAR, MAPQ, FLAG, TLEN (all cascade from primary placement).
- **Symptom:** bwa-rs picks the primary from the highest-scoring *chain*
  (`chains.max_by(|c| c.score)`, `src/alignment.rs`), where `chain.score` is a
  seed-length/gap heuristic from `chain_seeds`. bwa picks the region with the highest
  *extended* alignment score. The two are not comparable, so bwa-rs sometimes anchors a
  lower-scoring placement: e.g. `SRR7733443.249` is emitted `101S50M` (AS 50) while a
  candidate MEM DP-extends to AS 57 at another locus — the candidate *is* the alignment
  bwa reports. T-024 had to add a `sub <= score` guard (skip candidates scoring above the
  chosen primary) purely to stop this mis-selection from forcing MAPQ to 0.
- **Suspected cause:** `align_read` chooses `best_chain` before extension and never
  compares the *extended* scores of competing chains/candidates; `build_alignment` is run
  once on the winner. `chain_seeds` → `finalize_chain` also collapses each chain to one
  representative seed, discarding alternative anchors within a chain.
- **Resolution (proposed):** Extend the top-N chains/candidates (the T-024 candidate set
  already does this for `sub`), then select the primary as the highest *extended-score*
  region, with bwa's tie-breaks (lower coordinate, then strand). Remove the T-024
  `sub <= score` guard once the primary is guaranteed maximal. This will move POS/CIGAR on
  the sub2k set, so it must be measured against `uniq` for zero regression and against
  `sub2k` for net POS/CIGAR improvement.
- **Resolution (implemented):** `align_read` (`src/alignment.rs`) now keeps the chained
  SMEM as a `baseline` and extends the existing T-024 candidate set (unfiltered MEMs +
  `collect_short_seeds`, deduped by `(query_start, ref_start)`, capped at
  `MAX_SUB_CANDIDATES = 64`) into `placements`. The primary is the highest *extended*-score
  placement; the baseline wins exact score ties so uniquely-mapping reads only move when a
  candidate scores **strictly** higher (no coordinate/strand tie-break needed). The
  `if alt_res.score > result.score { continue; }` guard is removed — the primary is now
  maximal, so a distinct-locus candidate overlapping the primary in query space is a
  genuine secondary. `build_alignment` reads only `chain.mem`, so baseline and candidates
  extend on equal footing. New test `test_primary_selected_by_max_extended_score`.
- **Measured (chr1, `-k 19 -t 16`, `bench/compare_t025{,_summary}.md`):** **uniq placement
  improved with no regression** — POS 0, CIGAR 4 → **2**, TLEN 4 → **0**; the two changed
  reads (`100142`, `100189`) move *to* bwa-mem2's exact CIGAR/TLEN (`138M13S`/`146M5S`,
  the higher-scoring region is bwa's). **sub2k net placement improvement:** POS 40.7% →
  **38.0%** (−108), CIGAR 34.1% → **30.0%** (−167), RNAME/RNEXT 8.8% → **2.0%** (−272),
  FLAG/TLEN/SEQ/QUAL all down; identically-placed reads 2283 → **2405** (+122). MAPQ ticks
  up (sub2k 38.9% → 39.9%; uniq 2 → 4; identically-placed exact 76.7% → 73.3%) — the
  newly-correctly-placed reads are heavily-soft-clipped repeats bwa grades down via deeper
  region scoring / chain filtering; bwa-rs over-scores them (`bwa-rs lower` 5 → 20).
  Throughput/RSS unchanged (sub2k 37 s, 15.39 GB). All tests pass; clippy/fmt clean.
- **Residual (other root causes):** byte-matching the remaining POS/CIGAR and the MAPQ/XS
  *values* needs faithful `ksw_extend`/`mem_chain2aln` region extension (**T-026**) and
  `mem_chain_flt` chain pruning (**T-027**). T-025 makes the primary maximal over the
  current candidate set; T-026 makes each candidate's extended score/boundaries faithful.
- **Dependencies:** Builds on the T-024 candidate machinery. Feeds **T-026**.

### T-026: Faithful `mem_chain2aln` / `ksw_extend` region extension (end-bonus banded SW)
- **Status:** FIXED (branch `t026-ksw-extend`)
- **Severity:** high (the deep core — remaining 1–10 bp POS shifts and soft-clip boundaries)
- **Affected field(s):** POS, CIGAR, AS, NM, MD, MAPQ.
- **Symptom:** Even where bwa-rs and bwa agree on the locus, the extended alignment
  differs: bwa-rs scores lower and clips/shifts differently (e.g. `SRR7733443.623`
  bwa-rs AS 90 / `104M1D2M2D42M` vs bwa AS 124; the T-018/T-014 residual 1–10 bp POS
  shifts and `±` soft-clip-boundary CIGARs). bwa-rs's `affine_extend_forward` /
  `build_alignment` anchor at the seed start corner and extend outward; bwa's
  `ksw_extend` runs banded SW from each seed with a query-end *bonus* (`opt->w`,
  `zdrop`, end-bonus) and reconciles left+right extensions, which picks different
  boundaries and resolves scoring ties differently.
- **Suspected cause:** `src/alignment.rs::affine_extend_forward` / `extend_seed_backward`
  are a simplified banded affine DP without bwa's end-bonus, `zdrop`, or the
  `mem_chain2aln` per-chain multi-region reconciliation; `Scoring` tie-breaks
  (`gap_open`/`clip_penalty`) don't match `ksw` exactly on equal-score paths.
- **Resolution (proposed):** Port `ksw_extend` (banded SW with end-bonus + `zdrop`) and
  `mem_chain2aln` (extend each chain into one-or-more `mem_alnreg_t` regions), replacing
  the single-seed `build_alignment`. Reuse the T-024 candidate set as the region inputs.
  This is the largest single change and subsumes the T-018 (banded SW) and T-014 (scoring
  tie-break) residuals. Measure POS/CIGAR/AS on `sub2k`; guard `uniq` for zero regression.
- **Resolution (implemented):** Replaced the seed-start-anchored banded affine DP
  (`affine_extend_forward` / `extend_seed_backward`) on the `build_alignment` path with
  bwa's region-extension model (`src/alignment.rs`):
  (1) **`ksw_extend`** — a faithful port of bwa's `ksw_extend2` (scoring only): extends
  from the locked seed score `h0 = seed_len·match`, returning the local max
  `(score, qle, tle)` and the best query-end-reaching score `(gscore, gtle)`. Banded
  (`KSW_BAND_WIDTH = 100`) with bwa's per-row band shrink and `KSW_ZDROP = 100` X-drop.
  (2) **`ksw_global`** — banded affine global alignment with affine traceback, rendering
  the Eq/X/I/D CIGAR over the window `ksw_extend` determines (bwa's `ksw_global2`).
  (3) **`build_alignment`** now locks the exact seed `[qs,qe]×[rs,re]`, right-extends
  `query[qe..]`/`ref[re..]` and left-extends the reversed upstream slices (each with `h0`
  and `end_bonus = clip_penalty = 5`), chooses to-end vs local per side via
  `gscore + end_bonus > score`, assembles the window `[qb,qend]×[rb,rend]`, and runs
  `ksw_global` for the core CIGAR with soft-clipped flanks. The nm/score/2N→forward-remap/MD
  tail and `align_read`'s T-024/T-025 candidate-set + max-extended-score primary selection
  are unchanged. The old extension fns stay (`pub`, unit-tested) but off the hot path.
- **Measured (chr1, `-k 19 -t 16`, `bench/compare_t026{,_summary}.md`):** **zero
  regression on uniq, large sub2k improvement.** 300-uniq: POS 0, CIGAR 2, TLEN 0 (all
  unchanged); MAPQ 4 → **2** — the only two changed reads (`100142`/`100189`) move *to*
  bwa-mem2's MAPQ 60, POS/CIGAR/FLAG/TLEN byte-identical to T-025. sub2k: POS 38.0% →
  **31.1%** (−276), CIGAR 30.0% → **25.9%** (−163), RNAME/RNEXT 2.0% → **0.6%**, PNEXT
  −276, FLAG 36.0% → 33.8%, TLEN 48.4% → 44.2%, SEQ/QUAL 17.3% → 14.5%; identically-placed
  reads **2405 → 2578 (+173)**. Throughput/RSS unchanged (35 s, 15.72 GB). 286 unit tests
  pass (6 new `ksw_extend_*` / `ksw_global_*`); `cargo clippy -- -D warnings` and
  `cargo fmt --check` clean.
- **Residual (other root cause — follow-up):** identically-placed MAPQ exact-match dips
  73.3% → 69.5% with bwa-rs scoring **higher** on 719 (27.9%) and lower on only 67 (2.6%) —
  the newly-correctly-placed reads are heavily soft-clipped repeats bwa grades **down** via
  chain filtering (`mem_chain_flt`) and deeper region scoring. Closing it is **T-027**.
- **Dependencies:** **T-025** (max-score primary selection over the extended regions).

### T-027: Chain filtering (`mem_chain_flt`) to prune weak/contained chains
- **Status:** FIXED (branch `t027-chain-filtering`)
- **Severity:** medium
- **Affected field(s):** MAPQ, XS (and indirectly POS via the candidate set).
- **Symptom:** bwa-rs has no `mem_chain_flt`, so weak chains that bwa drops survive into
  the candidate set. This surfaced in T-024 as XS over-counting on unique reads
  (`SRR7733443.100004` bwa-rs XS 15 where bwa emits XS 0) — bwa discards the short,
  largely-contained chain before scoring it.
- **Suspected cause:** `chain_seeds` (`src/chaining.rs`) emits every greedy chain;
  there is no drop by chain weight / containment (bwa's `chain_drop_ratio = 0.50`,
  `min_chain_weight`, overlap-vs-better-chain test).
- **Resolution (implemented):** Ported bwa's chain-weight and chain-filter as pure
  functions in `src/chaining.rs`: `mem_chain_weight` (min of query-span and ref-span
  seed-union coverage, capped at `(1<<30)-1`); `build_candidate_chains` (groups
  collinear seeds into `Chain`s via bwa's `test_and_merge` diagonal-band test,
  `CHAIN_BAND_WIDTH = 100`, `MAX_CHAIN_GAP = 10000` — cross-strand/cross-locus seeds
  fall out because their reference offset breaks the band, contained seeds are
  absorbed); and `mem_chain_flt` (sort chains by weight desc, greedily keep, drop a
  chain that overlaps a kept higher-weight chain in query space by `>= CHAIN_MASK_LEVEL
  = 0.50` and is much weaker — `w < best_w * CHAIN_DROP_RATIO (0.50)` **and**
  `best_w - w >= 2 * min_seed_len`). `align_read` (`src/alignment.rs`) now runs
  `mem_chain_flt(build_candidate_chains(cand_mems))` on the unfiltered MEM +
  third-round (`collect_short_seeds`) candidate set and extends only the survivors.
  Because the baseline (chained-SMEM) placement is still pushed separately and the
  primary is the max-extended-score region, pruning narrows only the secondary
  (`sub`/XS) set — it never moves the primary placement, so uniq stays frozen.
- **Measured (chr1, `-k 19 -t 16`, `bench/compare_t027{,_summary}.md`):** **zero uniq
  regression** — SAM cols 1–10 byte-identical to the T-026 baseline (POS 0, CIGAR 2,
  TLEN 0, FLAG 4, MAPQ 2, SEQ/QUAL 2; 0 records differ on cols 1–10). The XS
  over-counting this ticket targets is closed: uniq records with non-zero `XS:i` drop
  **94 → 13**, and the documented read `SRR7733443.100004` moves `XS:i:19 → XS:i:0`
  (bwa-mem2 emits 0). sub2k placement is net-neutral-to-better (POS 1243 → 1234,
  PNEXT 1243 → 1234, TLEN 1768 → 1754, SEQ/QUAL 579 → 578; CIGAR 1036 → 1040, FLAG
  flat); identically-placed reads 2578 → 2582; MAPQ exact-match flat (69.5% → 69.3%),
  bwa-rs-lower 67 → 60. 347 tests pass (7 new chaining tests); `cargo clippy
  -- -D warnings` and `cargo fmt --check` clean.
- **Residual (other root cause — follow-up):** the bulk of the sub2k MAPQ over-scoring
  is the *inverse* of this ticket — bwa finds an *inexact* secondary locus (DP-extends a
  short shared seed across mismatches/indels) that bwa-rs's seeding never surfaces, so
  `sub = 0` → MAPQ 60. `mem_chain_flt` only removes candidates, so closing it needs
  richer approximate-secondary seeding (the **T-024** class), not chain filtering.
- **Dependencies:** Composes with **T-025**/**T-026**; landed independently as a
  candidate-set refinement.

---

## Phase 3: multi-mapper placement & pairing (T-028 FIXED; T-029 … T-032 OPEN)

After Phase 2 the `uniq` set is at parity but `sub2k` still diverges (POS 30.9%,
CIGAR 26.0%, MAPQ 39.5%, FLAG 33.8%, TLEN 43.9%). The Phase-1/2 residual notes
(T-021, T-022, T-024, T-026, T-027) all point at the same architectural gap:
**bwa-rs collapses each read to a single placement before pairing.** `align_read`
(`src/alignment.rs`) builds a transient candidate-region set (unfiltered MEMs +
`collect_short_seeds`, chain-filtered) only to derive `sub`/XS, then returns **one**
max-extended-score region; the `_mate` argument is ignored and the production
`align_batch` path aligns every read independently. The mate enters only *after*
placement — `rescue_mate` (one end unmapped) and `pair_mapq`/`mate_fields`. bwa
instead keeps the whole per-read region array (`mem_alnreg_v`) and, in
`mem_sam_pe`/`mem_pair`, **chooses each mate's locus jointly** to maximize the pair
score (both-mate score + insert-size log-likelihood vs an unpaired penalty). On a
multi-mapper whose read maps equally well to several loci, that joint choice is what
puts the read next to its mate — e.g. `SRR7733443.732#128` lands at POS 249,240,229
in bwa-rs but at 10,027 (beside its mate) in bwa-mem2.

Phase 3 is structured **investigate → implement**: **T-028** buckets the residual
records by root cause and quantifies each bucket, fixing the priority of the
implementation tickets (**T-029 … T-032**), which are written as hypotheses gated on
its findings. The `uniq` set remains the zero-regression gate throughout.

### T-028: Characterize & bucket the residual `sub2k` divergence (INVESTIGATION)
- **Status:** FIXED (investigation complete, 2026-06-04) — `bench/bucket_t028.pl`,
  `bench/compare_t028_buckets.md`
- **Severity:** high (gates and prioritizes all of Phase 3)
- **Affected field(s):** POS, CIGAR, FLAG, TLEN, MAPQ, RNAME (diagnosis only — no code change).
- **Goal:** Turn the aggregate `sub2k` percentages into a per-record root-cause
  breakdown so the implementation tickets are ordered by real impact, not hypothesis.
- **Method:**
  1. Start from the T-027 comparison (`bench/bwars_2n_2k_t027.sam` vs
     `bench/bwamem2_2k.sam`, `bench/compare.sh`). For every record that differs on
     POS/CIGAR/FLAG/MAPQ, classify it into one bucket:
     - **B1 mate-disambiguated placement** — bwa-rs places the read at a different
       (often distant or opposite-strand) locus of equal/near-equal score, while
       bwa-mem2's locus is the one concordant with the already-agreed mate. Detect:
       bwa-mem2 POS is within the insert window of the mate and bwa-rs POS is not, or
       strand flips (`111S40M` vs `40M111S`). Seed examples: `732#128`, `654#64`,
       `921#128`, `987#64`.
     - **B2 unmapped-vs-mapped** — bwa-rs emits `*` (RNAME/POS unset) where bwa-mem2
       maps the read. Seed examples: `182#64`, `182#128`, `247#64`.
     - **B3 MAPQ-only over-scoring** — POS+CIGAR agree but bwa-rs MAPQ > bwa-mem2
       (missing inexact secondary locus → `sub = 0`). Seed examples: `732`, `730`,
       `91` (MAPQ 10/17/40 vs 0). This is the T-024/T-027 residual.
     - **B4 soft-clip-boundary / small POS shift** — same locus, CIGAR boundary or
       1–10 bp POS shift (the T-026 residual tail).
     - **B5 other** — anything unclassified (record and inspect).
  2. Report, per bucket: record count, share of the POS / CIGAR / FLAG / MAPQ
     divergence it accounts for, and 2–3 worked examples with both tools' fields.
  3. Cross-check B1 by re-running bwa-rs with pairing forced off in bwa-mem2 (or by
     inspecting whether bwa-mem2's chosen locus is the mate-concordant one) to confirm
     the bucket is genuinely pairing-driven and not a scoring gap.
- **Deliverable:** `bench/compare_t028_buckets.md` with the table above, and this
  ticket updated with the bucket sizes + the resulting priority order for T-029…T-032
  (e.g. if B1 dominates, T-029+T-030 lead; if B3 dominates, T-032 leads).
- **Dependencies:** None (read-only analysis on existing SAMs). Gates T-029 … T-032.
- **Findings (`bench/compare_t028_buckets.md`, T-027 sub2k, 4000 read records):**
  1765 records are identical on the compared columns; **2235 diverge.** Bucketed
  (each record assigned to one bucket by priority B2 → B1 → B4 → B3 → B5):

  | Bucket | Meaning | Count | % of divergent | MAPQ A>B | MAPQ A<B |
  |---|---|---|---|---|---|
  | **B1** | different-locus placement | **936** | **41.9%** | 334 | 83 |
  | B3 | MAPQ-only (placement identical) | 793 | 35.5% | 733 | 60 |
  | B4 | same-locus ≤100 bp shift / soft-clip boundary | 395 | 17.7% | 305 | 44 |
  | B2 | unmapped-vs-mapped | 74 | 3.3% | 0 | 23 |
  | B5 | FLAG/PNEXT-only at identical placement | 37 | 1.7% | 0 | 0 |

  **The mate-disambiguation hypothesis is confirmed.** Of the 936 B1 records,
  **827 (88.4%) are proper-pair-flag flips — bwa-mem2 marks the pair proper (0x2)
  and bwa-rs does not** — i.e. bwa places the read at the mate-concordant locus and
  bwa-rs at a discordant equal-scoring copy. RNAME never differs (single contig);
  537 B1 records flip strand, 399 are same-strand >100 bp away. B5 (37) is the
  downstream cascade: same POS/CIGAR/MAPQ as bwa-mem2 but a different proper-pair
  FLAG / PNEXT because the *mate* is a B1 record. B2 (74) is the same root cause at
  the extreme — bwa pairs/rescues the read to a locus bwa-rs leaves unmapped.
  → **B1 + B5 + B2 ≈ 47% of the divergence is the single mate-aware-placement gap.**

  **B3 (35.5%, 733 over-scoring) is the inexact-secondary MAPQ residual** — identical
  placement, bwa-rs MAPQ far higher (e.g. `100#64` 80M71S MAPQ 54 vs 5) because bwa's
  inexact second-best locus gives `sub > 0` while bwa-rs's seeding yields `sub = 0`.
  Pure MAPQ; no placement move. **B4 (17.7%)** is a mix: the pair shifts together by
  tens of bp (e.g. `1#64` POS 10043 vs 9999, same CIGAR shape) with MAPQ also
  over-scored — part will resolve with mate-aware placement, the rest is the T-026
  extension-fidelity tail. MAPQ over-scoring totals 334 (B1) + 733 (B3) + 305 (B4)
  = **1372 records ≈ the 39.5% aggregate MAPQ divergence**, so fixing placement
  (T-030) also reclaims the 334 + ~305 placement-driven over-scores, leaving B3 as the
  standalone MAPQ work.
- **Resulting priority for Phase 3 (evidence-based):**
  1. **T-029 → T-030 (mate-aware placement)** — *highest*. Directly closes B1 (41.9%),
     cascades into B5 (1.7%) and most of B2 (3.3%), and reclaims ~640 placement-driven
     MAPQ over-scores. The 88.4% proper-flip rate makes this the dominant lever.
  2. **T-032 (inexact-secondary MAPQ)** — *second*. B3 (35.5%, 733 over-scores); pure
     MAPQ, independent of placement, the largest standalone MAPQ contributor.
  3. **T-031 (unmapped recovery)** — B2 (3.3%); likely folds into T-030's pairing/rescue.
  4. **B4 residual** — re-measure after T-030; the remainder is T-026-class extension
     fidelity, not a new ticket yet.

### T-029: Expose the per-read candidate-region array (`mem_alnreg_v`)
- **Status:** OPEN
- **Severity:** high (structural prerequisite for T-030/T-031/T-032)
- **Affected field(s):** none directly — enables the tickets that move POS/FLAG/MAPQ.
- **Symptom:** `align_read` already computes the extended candidate regions (it needs
  them for `sub`/XS) but discards all but the single chosen primary. The paired driver
  in `src/main.rs` therefore has nothing to choose among — it can only rescue or
  recompute MAPQ on a fixed placement.
- **Suspected cause:** `Aligner::align_read` → returns one `AlignmentResult`;
  `ParallelAligner::align_batch` and the `align_all` flat-list path (`src/main.rs`)
  are built around one-region-per-read.
- **Resolution (proposed):** Return the deduped, score-sorted candidate-region list
  (bwa's `mem_alnreg_v`: position, strand, score, query span, CIGAR, `sub`/`sub_n`,
  `frac_rep`) from a new entry point (keep the single-result API as a thin wrapper
  that takes region 0, so the single-end path and all current tests stay byte-stable).
  Plumb the array through `align_batch`. **Zero output change on its own** — verified
  by the `uniq` gate and `sub2k` byte-identity to T-027; this ticket only changes the
  data available downstream.
- **Dependencies:** Feeds **T-030**, **T-031**, **T-032**. Independent of T-028's
  outcome (the array is needed by any of the fixes), so it can start in parallel.

### T-030: Paired-end mate-aware placement selection (`mem_sam_pe` / `mem_pair`)
- **Status:** OPEN
- **Severity:** high (**confirmed largest driver** by T-028: B1 41.9% of divergence,
  88.4% proper-pair-flag flips, plus the B5/B2 cascade)
- **Affected field(s):** POS, RNAME, FLAG, TLEN, PNEXT, RNEXT, SEQ/QUAL orientation, MAPQ.
- **Symptom:** Bucket **B1** — on multi-mappers bwa-rs picks a per-read max-score
  locus that need not be the one pairing concordantly with the mate; bwa-mem2 picks
  the pairing-maximizing combination. Cascades into FLAG (proper-pair bit), TLEN,
  PNEXT/RNEXT, and SEQ/QUAL orientation.
- **Suspected cause:** No `mem_pair` pass. `src/main.rs` pairs two already-fixed
  single placements; `align_read` ignores `_mate`.
- **Resolution (proposed):** Port `mem_sam_pe`'s region pairing over the T-029 arrays:
  for each candidate region of mate 1 × mate 2, score the pair as `score1 + score2 +
  insert_bonus(d)` (reuse T-022's `insert_bonus`/`pair_mapq`), pick the max-scoring
  concordant combination above the unpaired alternative (`PEN_UNPAIRED`), and emit
  those two placements; fall back to the current independent best when no pair beats
  unpaired. The T-022 paired-MAPQ recompute then runs on the *selected* pair (with the
  real competing-pair `subo`, closing the T-022 `100124` overshoot residual). Guard:
  `uniq` stays byte-identical (unique pairs have one region each, so the selection is a
  no-op there); measure `sub2k` POS/FLAG/TLEN.
- **Dependencies:** **T-029** (region arrays). Subsumes part of **T-022**'s residual.

### T-031: Unmapped-vs-mapped recovery on multi-mappers
- **Status:** OPEN
- **Severity:** medium (T-028: B2 = 3.3% of divergence; likely folds into T-030's pairing/rescue)
- **Affected field(s):** RNAME, POS, FLAG, RNEXT, PNEXT (records bwa-rs leaves `*`).
- **Symptom:** Bucket **B2** — bwa-rs emits `*` where bwa-mem2 places the read
  (`182#64`, `182#128`, `247#64`). Either the read's only acceptable locus scores
  below `-T` on its own but is rescued/accepted via the mate, or seeding surfaces no
  region at all for it.
- **Suspected cause:** `result.score < self.min_score → unmapped` is applied per read
  with no pair context; `rescue_mate` only fires when the *other* mate already mapped,
  and mate rescue runs after independent placement rather than as part of pairing.
- **Resolution (proposed):** Determine via T-028 whether these are (a) pair-rescue
  cases (mate maps, orphan should be `mem_matesw`-rescued — check the existing
  `rescue_mate` trigger/threshold) or (b) sub-`T` placements that bwa keeps because the
  *pair* clears threshold. Fix accordingly: extend rescue to fire from the T-030
  selected pair, and/or apply bwa's pair-conditioned acceptance. Guard `uniq` (no
  unmapped records there) for zero regression.
- **Dependencies:** **T-029**; likely **T-030** (pairing selects the locus to rescue toward).

### T-032: Residual MAPQ over-scoring from inexact secondary loci
- **Status:** OPEN
- **Severity:** high (**confirmed second driver** by T-028: B3 = 35.5% of divergence,
  733 MAPQ over-scores; the largest standalone MAPQ contributor)
- **Affected field(s):** MAPQ, XS.
- **Symptom:** Bucket **B3** — POS+CIGAR agree but bwa-rs MAPQ is too high because the
  second-best locus is *inexact* (bwa DP-extends a short shared seed across
  mismatches/indels; bwa-rs's exact-MEM + third-round seeding never surfaces it), so
  `sub = 0` → MAPQ 60. Documented across T-021/T-024/T-027 as the inverse of chain
  filtering — only richer secondary discovery can close it.
- **Suspected cause:** The candidate set fed to `sub`/XS still misses inexact
  second-best alignments; `collect_short_seeds` (third round) catches some but not the
  cross-mismatch/indel loci bwa reaches by extension.
- **Resolution (proposed):** Re-measure after T-029…T-031 (mate-aware placement may
  already supply the competing locus via the pair). If a gap remains, widen secondary
  discovery toward bwa's region set — e.g. extend more sub-optimal chains / seed the
  primary's flanks against repeat copies — feeding `sub`. Guard `uniq`; target the
  identically-placed-read MAPQ exact-match rate on `sub2k`.
- **Dependencies:** **T-029**; re-scoped by **T-028** and the T-030/T-031 outcome
  (may shrink or close before standalone work is needed).
