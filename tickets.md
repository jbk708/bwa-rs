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
- **Status:** OPEN
- **Severity:** medium
- **Affected field(s):** MAPQ, XS (optional tag).
- **Symptom:** Split from T-021. After T-021's re-seeding + `mem_mark_primary_se`-style
  `sub`/`sub_n` + `frac_rep`, the bulk of T-015's over-scoring reads still report MAPQ
  60 where bwa-mem2 reports 0–52 (heavily soft-clipped, moderately repetitive reads).
  Their second-best locus is *inexact* — bwa places it by DP-extending a short shared
  seed across mismatches/indels (e.g. `SRR7733443.26`'s 72 bp core maps at a second
  locus as `96M55S` / `65M1I35M1D26M24S`, not an exact 72 bp MEM). bwa-rs's exact-MEM
  seeding never surfaces that locus, so `sub = 0` → MAPQ 60 and XS values diverge.
- **Suspected cause:** T-021 surfaces secondaries only where the repeat is an *exact*
  MEM at ≥2 loci. bwa's candidate set additionally includes chains whose seed is short
  and shared but whose full alignment is inexact (found via `ksw`/banded extension at
  each chain), producing a richer `sub`/XS set.
- **Resolution:** Build alignment candidates from *all* chains (not just exact-repeat
  reseeds) and let banded extension score each at its locus, so an inexact second-best
  contributes to `sub`/XS; tighten the `mem_mark_primary_se` clustering to bwa's. This
  unblocks the remaining T-015 MAPQ residual and XS *value* byte-matching (T-011).
  Overlaps with **T-022** (faithful pairing score) and is gated on **T-023** for
  large-set measurement.
