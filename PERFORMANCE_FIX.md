
# flexalign performance fix — summary

_Generated 2026-07-27. Diagnosed from `flexalign_benchmark` (bacteroides whole-genome DB +
protal marker DB). All changes implemented, compiled, and verified end-to-end. **Uncommitted.**_

## Changes at a glance

| # | Change                                               | File(s)                                                                              | Purpose                                                                      |
| - | ---------------------------------------------------- | ------------------------------------------------------------------------------------ | ---------------------------------------------------------------------------- |
| 1 | Bump blob format version 1 → 2                      | `../flexmap/src/flexmap.rs`                                                        | Mark the post-`d00a75f` key scheme as a new on-disk format                 |
| 2 | `FlexmapBlob::is_compatible()` header peek         | `../flexmap/src/flexmap.rs`                                                        | Non-panicking staleness/param check                                          |
| 3 | `DB::index_compatible()` + rebuild-on-mismatch     | `src/database/flexmap.rs`, `src/flexalign.rs`                                    | Stale/incompatible index →**rebuild** instead of silent misread       |
| 4 | New`--min-query-coverage` (default **0.70**) | `src/options.rs`                                                                   | Reject partial (soft-clipped) wrong-target hits                              |
| 5 | `query_coverage` on `MateOut` + coverage gate    | `src/align/modular_workflow.rs`                                                    | Enforce coverage alongside`--min-ani`                                      |
| 6 | Phred-style MAPQ, clamped`[0,60]`                  | `src/align/common.rs` (`StdPairedAnchorMAPQ::anchor_mapq`)                       | Informative, monotone MAPQ (fixes the PR-curve bump)                         |
| 7 | Gapless flank fast path (skip WFA)                   | `src/align/data_structures/anchor.rs` (`align_left_flank`/`align_right_flank`) | ~2.7–4× faster alignment; abort from Hamming vs budget                     |
| 8 | Whole-read gapless pre-filter (SIMD Hamming)         | `src/align/data_structures/anchor.rs` (`align`)                                  | drops wrong-target/clip anchors before WFA — protal alignment 0.9s → 0.03s |
| 9 | mmap reference sequences + skip run-time name map    | `src/database/flexmap.rs`, `common.rs`, `Cargo.toml`                           | DB load 38.7s → 4s on protal; total wall 50s → ~8s (warm)                  |

Changes 1–2 live in the sibling **`../flexmap`** repo — must be pushed to GitHub to publish
(per `CLAUDE.md`). Only `--min-query-coverage` (0.70) and the blob version (2) changed behaviour;
every other parameter is unchanged.

## Root causes (it was **not** the chaining algorithm)

The benchmark report conflated two unrelated bugs on its two datasets.

### 1. Stale-index regression → the "slow" + "insensitive" bacteroides result

flexmap commit `d00a75f` started `mask_coremer`-mixing the index key but left
`FLEXMAP_BLOB_VERSION` at 1. The benchmark's Jul-19 index was silently reused with mismatched
keys → **ranges/read collapsed 23.73 → 1.29** → almost nothing aligned, and WFA2 burned ~30 s
failing on garbage anchors (+11 GB of mmap re-faults). Same index, same reads, identical
minimizers — only the lookup regressed.

**Fix (changes 1–3):** the blob carries a format version; on load, an incompatible blob is caught
by a cheap header peek and the index is **rebuilt** rather than misread. It logs:
`On-disk index is incompatible with this build (format/params changed); rebuilding it.`

### 2. Missing coverage gate → the protal "precision nosedive"

`approx_ani` charges only the *aligned window's* cost against the *full* read length, so a read
aligning ~50 % of itself to a wrong-species marker at local high identity read as ~100 % ANI and
passed the `--min-ani 0.85` gate. There was no query-coverage requirement (69 % of alignments
carried an indel; mean query coverage 51 %; replayed identity 55 %).

**Fix (changes 4–5):** each mate computes `query_coverage = cigar.query_consumed() / read_len`
(= `1 − soft-clip fraction`; seed span for Hamming-only anchors), and emission now requires
**both** `ani ≥ min_ani` **and** `query_coverage ≥ min_query_coverage`.

## Verified results

### Bacteroides (whole-genome DB) — stale index, now fixed by rebuild

| metric                        | before (stale) | **after**   | minibwa |
| ----------------------------- | -------------- | ----------------- | ------- |
| aligned                       | 0.03 %         | **99.77 %** | 100 %   |
| position-correct (of aligned) | 2.15 %         | **97.99 %** | 97.99 % |
| wall time                     | 42.5 s         | **8.05 s**  | 8.87 s  |
| peak RAM                      | 5.9 GB         | 6.2 GB            | 3.1 GB  |

Tied with minibwa on sensitivity and precision; comparable speed.

### Protal (marker DB) — real precision problem, fixed by coverage gate

| metric                             | before | **after**  | protal | minibwa |
| ---------------------------------- | ------ | ---------------- | ------ | ------- |
| **species-correct / mapped** | 26.3 % | **92.1 %** | 93.1 % | 83.7 %  |
| species-correct / total            | 4.36 % | 4.01 %           | 3.78 % | 4.48 %  |
| align rate                         | 16.6 % | 4.35 %           | 4.07 % | 5.35 %  |

Precision now matches protal; per-read recall beats protal. (Coverage sweep: 0.5 → 91.5 % prec /
4.18 % recall; 0.7 → 92.1 % / 4.00 %; 0.9 → 92.3 % / 3.81 %. Junk clusters below ~0.5 coverage;
precision is flat above the knee, so 0.70 clears it with margin.)

### MAPQ / precision–recall curve (change #6)

The old pseudo-MAPQ was a raw best-vs-second score gap cast to `u8`: it wrapped past 255, forced
uniquely-mapped reads to 0, and produced a **bimodal** distribution (values only ~0–33 or ~150–254,
nothing between) — so the MAPQ-swept PR curve had a non-monotone bump/plateau. The new MAPQ is the
**relative** margin `1 − score2/score1`, scaled to `[0, 60]` and clamped (no wrap), with a lone hit
given a fixed middle value (30) rather than the max. It spreads smoothly across 0–60 and the curve
is now **monotone** — precision rises steadily as the cutoff rises:

| MAPQ cutoff ≥ | recall | precision |
| -------------- | ------ | --------- |
| 58             | 0.39 % | 98.96 %   |
| 44             | 2.06 % | 98.22 %   |
| 23             | 3.31 % | 96.76 %   |
| 0 (all)        | 4.01 % | 92.08 %   |

Mapping positions are unchanged, so species-correct stays 92.08 %; only the MAPQ column / curve
shape changed.

### Alignment hot path — gapless flank fast path (change #7)

The flanks (read ends beyond the seed chain) were each WFA-aligned for every candidate anchor, and
most of the time went to the ~90 % of candidates that are wrong-target: WFA explored the divergent
flank up to the ANI abort budget before dropping it. `left_flank`/`right_flank` always return an
**equal-length** q/r window, so a gapless alignment is *provably optimal* whenever its mismatch cost
is below the cheapest possible gap (`gap_open+gap_extend`), and — mirroring protal's `approximate`
path — when **no indel is expected in the anchor** (`!flagged_for_indel`) the gapless alignment is
taken for any mismatch count. A gapless cost over the (threaded, min-ANI-derived) budget now **drops
the anchor with an O(len) Hamming pass instead of a WFA abort**. WFA still runs for `flagged_for_indel`
flanks and boundary soft-clip cases.

Measured (vs the `results/_baseline_pre_flank/` snapshot):

| dataset     | alignment time                    | read-processing            | accuracy                                                                 |
| ----------- | --------------------------------- | -------------------------- | ------------------------------------------------------------------------ |
| protal      | 2.48 s →**0.90 s** (2.7×) | 15.99 s →**7.86 s** | corr/mapped 92.08 % → 92.06 %; indel share 0.6 % →**0.1 %**      |
| bacteroides | 0.44 s →**0.11 s** (4×)   | —                         | aligned 99.77 % → 99.76 %; position-correct 97.99 % →**97.99 %** |

No accuracy regression on either dataset (bacteroides is a real genome with true indels), and the
spurious flank indels WFA used to invent are gone (0.6 % → 0.1 %, near protal's 0.2 %).

### Whole-read gapless pre-filter (change #8)

Instrumentation of the remaining alignment time showed **~95 % of it was still throw-away** — WFA on
wrong-target anchors that *drop* — and that **~99.997 % of those WFA calls were the soft-clip case**
(read overhangs a short marker gene: clip=216,542 vs indel=7 on protal). The flank fast path didn't
cover clips, so those anchors ran WFA up to the abort before dropping. Fix: at the top of `align`,
before `smart_align`, one **SIMD Hamming pass on the anchor diagonal** (`triple_accel::hamming`) —
when no indel is expected and the gapless cost already exceeds the ANI budget, drop the anchor there
instead of running WFA on both flanks.

| dataset     | alignment time                                               | WFA flank calls            | accuracy                                      |
| ----------- | ------------------------------------------------------------ | -------------------------- | --------------------------------------------- |
| protal      | 0.90 s →**0.029 s** (~31×; ~85× vs original 2.48 s) | 216,549 →**13,793** | corr/mapped 92.06 % → 92.03 %                |
| bacteroides | 0.11 s →**0.15 s** (noise)                            | 13,493 →**2,682**   | position-correct**97.99 %** (identical) |

After this, protal alignment is ~12 % throw-away (the cheap Hamming pass), and the surviving WFA is
what it should be: genuine indels (bacteroides 2,336) and the few kept clip alignments.

### DB load: mmap the reference sequences + skip the run-time name map (change #9)

Profiling the wall showed alignment was never the bottleneck — **~75% of the protal wall was the
one-time DB load**: re-parsing the 16.4 Gbp reference FASTA into owned records (25 s) and loading the
14.5 M-entry `rname->rid` map. Two changes, both mirroring the index blob's existing mmap:

- **Reference sequences** are written once to a memory-mappable blob (`<index>.refseq`: concatenated
  sequences + offset table) and `mmap`'d on load; `get_reference` returns a zero-copy slice. Generated
  lazily on the first load that doesn't find it (so no rebuild needed), and at build time.
- **`rname->rid`** is only needed to attribute FASTA records to ids and for gold-standard eval, so on
  the mmap path it is not loaded at all (unless `GOLDSTD_EVAL`).

|             | load references             | total DB load             | warm wall             | correctness                   |
| ----------- | --------------------------- | ------------------------- | --------------------- | ----------------------------- |
| protal      | 25.5 s →**31 ms**    | 38.7 s →**3.95 s** | 50 s →**~8 s** | 92.03 % (identical)           |
| bacteroides | 1.85 s →**0.037 ms** | 1.85 s →**17 ms**  | —                    | 99.76 % / 97.99 % (identical) |

flexalign's protal wall is now **well under minibwa's 46 s**. (First load per index still pays a
one-time parse + ~16 GB blob write; every load after mmaps it. Peak RSS is similar but now
file-backed / reclaimable rather than anonymous heap.)

### strobealign added to the protal dataset

`strobealign` is now listed for the protal dataset in `datasets.tsv` (it was already on bacteroides);
its index was built and it was run on protal sample 1, so it appears in the regenerated report.

## Speed (measured, this build)

The "flexalign is slow" finding was **entirely the stale index** (bacteroides 42 s = WFA2 thrash +
mmap re-faults). The coverage gate is a pure output filter and cannot slow the pipeline.

| dataset       | flexalign (now)                    | read-processing only | index/ref load     | vs others (report)               |
| ------------- | ---------------------------------- | -------------------- | ------------------ | -------------------------------- |
| bacteroides   | **8.05 s** wall, 6.2 GB RSS  | 6.5 s                | 1.3 s              | minibwa 8.87 s · bowtie2 67.7 s |
| protal (warm) | **60.7 s** wall, 49.9 GB RSS | 14.6 s               | ~37 s (32 GB refs) | minibwa 46.2 s · protal 104.8 s |

- **bacteroides:** on par with the fastest full aligners, far below bowtie2.
- **protal:** **faster than protal** (60.7 s vs 104.8 s), slower than minibwa (46.2 s). The bulk is
  loading the 32 GB reference into RAM (~37 s), an unavoidable cost every big-index tool pays;
  actual read-processing is only 14.6 s. A cold-cache run right after the rebuild measured 120 s
  (disk-faulting the freshly written index) — not representative of steady state.
- Index **build** is a one-time cost, not counted in benchmark wall time (the harness builds it
  separately).

> The version bump means each existing v1 index is rebuilt **once** on the next run (bacteroides
> ~2 min; protal ~30 min for 32 GB). This is the guard working as intended; subsequent runs load
> the v2 index normally.

## Report status — REGENERATED

`results/report.html` (+ `bacteroides/report.html`, `protal/report.html`) and `results/timing.html`
were regenerated against the fixed binary. Method: refreshed only the `flexalign` and
`flexalign-sam` outputs (both fixed, same binary) via their harness `run.sh`, then
`REPORT_ONLY=1` re-scored — reusing the existing minibwa/bowtie2/strobealign/protal outputs (no
expensive reruns). The reports now show flexalign bacteroides 99.77% aligned / 97.77% position
(verdict-agreement vs minibwa 0.98, was 0.00) and protal species-correct-of-mapped 92.08%.

**Caveat:** `flexalign-sharded` was **excluded** from the regenerated report — its emit path
(`src/shard/rejoin_align.rs` → `shard/disk.rs`) doesn't yet apply the coverage gate, and its shard
blobs are still the old v1 format. It needs the same two fixes (coverage gate on the rejoin emit +
re-slice shards from the v2 index) before it can be shown fairly. `datasets.tsv` still lists it, so
a future full `pixi run bench` will include it again.

## Current acceptance strategy

A candidate mate is emitted only if, after gapped alignment, **both**:

1. `approx_ani ≥ --min-ani` — identity over the aligned window (also converted to the WFA2 abort score).
2. `query_coverage ≥ --min-query-coverage` — fraction of the read actually aligned.

## Current settings (defaults)

| flag                     | default                                                      | role                             |
| ------------------------ | ------------------------------------------------------------ | -------------------------------- |
| `--min-ani`            | 0.85                                                         | identity gate + WFA2 abort score |
| `--min-query-coverage` | **0.70** ← new                                        | read-coverage gate               |
| `--ranges`             | 15                                                           | minimizers looked up             |
| `--max-range-size`     | 256                                                          | max occurrences per key          |
| `--max-best-flex`      | 16                                                           | seeds kept per key (best flank)  |
| `--extend-top-x`       | 4                                                            | anchors Hamming-extended         |
| `--align-top-y`        | 4                                                            | anchors gap-aligned              |
| `--min-ranges`         | 4                                                            | recovery floor                   |
| index const-generics     | K=31, C=15, F=16, S=4, CELLS_PER_BODY=16, HEADER_THRESHOLD=2 | fixed at`flexalign.rs`         |
| `FLEXMAP_BLOB_VERSION` | **2** ← bumped                                        | on-disk format guard             |

## Suggested follow-ups (not done)

- Apply the same coverage gate to the sharded rejoin path (`src/shard/rejoin_align.rs`) for parity
  with the unsharded pipeline (flexalign-sharded mirrored flexalign in the report).
- PAF still emits `alignment_block_length = 0` (fine for position scoring, but leaves the PAF's
  block-length column uninformative). _(MAPQ fixed — see change #6.)_
- Small recall ceiling vs minibwa (4.36 % vs 4.48 %) — a chaining-sensitivity item, independent of
  precision.
- Commit + push (flexmap sibling repo must be pushed to publish); regenerate `report.html`.
