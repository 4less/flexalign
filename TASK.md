# TASK — push flexalign to best-in-class sensitivity, precision and speed

_Self-assigned brief. Written 2026-08-04 on branch `Optimization` (identical to `shards`:
`git rev-list --left-right --count shards...Optimization` = `0 0`)._

## Scope directive (2026-08-04)

**`flexalign-batch` is the contender under optimisation — work on it solely.** The other
flexalign variants (`flexalign`, `-paf`, `-sharded*`) are not separate tuning targets: they must
**resolve anchors with the exact same code**. So every accuracy change lands in the shared
pipeline (`src/align/modular_workflow.rs` + `src/align/process/*`), never in a variant-specific
path, and the variants are expected to move together. `--extend-top-z` has exactly one call site,
which is the shape to preserve. A variant diverging in accuracy (see workstream D) is a bug in
that variant, not a tuning opportunity.

Benchmark runs therefore use `TOOLS=flexalign-batch`; other contenders are reference points read
from existing results, not re-run.

## The task

Improve `flexalign` so that, on the `flexalign_benchmark` suite, it is simultaneously

1. **the most sensitive** contender (highest per-read recall of correctly-placed reads),
2. **the most precise** contender (highest correct-of-mapped), and
3. **the fastest** contender (lowest wall time, at competitive RSS),

**without regressing any metric on any dataset.** Every change must be justified by a
before/after benchmark run, not by reasoning alone.

Hard constraints:

- **Do not modify `../flexmap`.** Any change that would require it — index key scheme, syncmer
  selection, blob layout, build-side canonicalisation — is written up in **`.claude/FLEXMAP.md`**
  as a proposal (problem, exact call sites, proposed change, expected effect, invalidation cost)
  and left unimplemented. Note that build/query symmetry means most index-side ideas (e.g. the
  whole of `TODO.md`'s syncmer rework) fall under this rule: the query side lives here, the build
  side lives in flexmap, and they must move together. So: `TODO.md` items 1–3 → `.claude/FLEXMAP.md`.
- `kmerrs` / `bioreader` are likewise sibling repos — treat a change there as the same kind of
  cross-repo commitment and document it before doing it.
- **stdout is the alignment stream.** No `println!`/color/log on stdout, ever
  (`modular_workflow.rs` has historically leaked debug prints into PAF/SAM).
- Nightly toolchain; `cargo +nightly build --release` is the only meaningful build for timing.

## Where I stand (baseline to beat)

From `flexalign_benchmark/results/*/summary.tsv` as of this writing. `recall_pct` =
correctly-placed reads / all reads; `correct_pct` = correct-of-mapped (precision).

> **Provenance warning.** These numbers were read from `results/*/summary.tsv` at 22:06 on
> 2026-08-04, *before* any run of this task's own. A benchmark run **overwrites the row of whatever
> contender it re-runs**, so the live file stops matching this table as soon as work begins — the
> protal `flexalign-batch` row was overwritten at 22:11 by the first (cold, invalid) run of the
> tie-extension change. The table here is the baseline of record; snapshots live in
> `results/_snapshot_*`. Snapshot `results/<ds>/` **before** any run that touches a baseline row.

### bacteroides (whole-genome DB, 1 sample)

| tool | align% | precision% | **recall%** | wall s | RSS MB |
| --- | --- | --- | --- | --- | --- |
| flexalign | 99.76 | 98.28 | **98.04** | **4.83** | 6003 |
| flexalign-paf | 99.76 | 98.28 | 98.04 | 4.70 | 5984 |
| flexalign-sharded | 99.70 | 98.30 | 98.01 | 12.54 | 3254 |
| minibwa | 100 | 98.28 | **98.28** | 6.88 | 3094 |
| bowtie2 | 100 | 98.25 | 98.25 | 60.47 | 1211 |
| strobealign | 100 | 98.01 | 98.01 | 5.72 | 4795 |

→ **Fastest already; recall is 0.24 pp short of minibwa, and all of that gap is the 0.24 % of
reads flexalign does not align at all.** Precision is at parity. RSS is 2× minibwa.

### protal (marker DB, 10 samples, species-level scoring)

| tool | align% | precision% | **recall%** | wall s | RSS MB |
| --- | --- | --- | --- | --- | --- |
| flexalign | 4.82 | 87.28 | 4.204 | 106.1 | 47709 |
| flexalign-batch | 4.85 | 86.93 | 4.214 | **38.2** | 48392 |
| flexalign-sharded-s2 | 4.76 | 87.79 | 4.182 | 554.1 | 17399 |
| flexalign-sharded-s2-batch | 4.76 | 87.79 | 4.182 | 101.9 | 17926 |
| minibwa | 5.43 | 80.82 | **4.391** | 212.8 | 36551 |
| protal | 4.21 | **89.15** | 3.756 | — | — |
| strobealign | 81.01 | 5.60 | 4.533 | 1135.4 | 76453 |

→ **This is the dataset that decides the task.** flexalign is the fastest by 5.6× but sits
_between_ the two references on both axes: recall 4.204 vs minibwa 4.391 (−4.3 % relative),
precision 87.28 vs protal 89.15 (−1.9 pp). Winning means **recall ≥ 4.40 and precision ≥ 89.2 at
the same time**, at ≤ 38 s.

⚠ strobealign's 4.53 recall / 5.6 % precision is a **scoring artifact**, not a target: it emits
~141 k bare ~25-base strobemer anchors that gated tools never emit. Do not chase its recall number
and do not "fix" precision by matching its emission policy. Any comparison against it needs a
common minimum-aligned-bases floor.

### markersim (marker DB, reference-level truth, 20 M pairs)

| tool | align% | precision% | **recall%** | wall s | RSS MB |
| --- | --- | --- | --- | --- | --- |
| flexalign | 99.12 | **98.92** | 98.04 | 98.2 | 48580 |
| flexalign-batch (`--mate-rescue`) | 99.58 | 98.75 | 98.34 | **82.6** | 48581 |
| minibwa | 100 | 98.76 | **98.76** | 251.4 | 38423 |
| strobealign | 100 | 98.76 | 98.76 | 242.2 | 80291 |

→ Fastest by 3×, best precision, **recall 0.42–0.72 pp short**. `--mate-rescue` buys +0.30 pp
recall for −0.17 pp precision — the right direction, not yet the whole gap.

## Diagnosis I already trust (do not re-derive)

- **The recall gap is a half-pair loss.** 86 % of the mates flexalign misses have their *partner*
  emitted; the funnel counters were pair-level and blind to it. `--mate-rescue` (lazy, score-frac
  selected) closes part of it. See `memory/flexalign-recall-halfpair.md`.
- **Acceptance is `end_to_end` + a minimum of 35 *aligned bases*** (`--min-query-coverage 35`);
  the old 0.70 coverage *fraction* is gone. This is what took protal precision 26 % → 92 % under
  the old scoring. See `memory/flexalign-precision-protal.md`.
- Perf history and the nine landed changes are in `PERFORMANCE_FIX.md`. Don't redo them.

## Workstreams

Ordered by expected payoff. Each one ends with a benchmark run and a line in the results log.

### A. Sensitivity — close the half-pair gap (biggest lever)

- Make `--mate-rescue` good enough to be the **default**, then default it. It currently costs
  precision (98.92 → 98.75 on markersim); find out whether the loss is rescue *placements* or
  rescue *acceptance* (`--mate-rescue-score-frac`, `--mate-rescue-margin`,
  `--min-query-coverage-mate`) and remove it rather than trading it away.
- `STRATEGY.md` (a): `--extend-top-z` truncates at a fixed count. **Extend past the cutoff while
  anchors are score-tied**, stopping at the first strictly worse anchor. A tie broken by list
  order is a coin-flip that costs both recall and precision.
- `STRATEGY.md` (b): rescue a mate whose expected position — from the observed insert-size
  distribution — still lands on the reference, even when the index reported nothing there.
  (`c92c8a8` started this; check what it does *not* yet cover.)
- Attack the 0.24 % (bacteroides) / 0.88 % (markersim) of reads that produce **no alignment at
  all**. Take a sample, replay them, and classify: no seeds → seed extraction; seeds but no anchor
  → chaining; anchor but gated → acceptance. Fix the dominant class only.

### B. Precision — reach protal's 89.2 % on the marker DB without giving back recall

- The gate is a single global `min-ani` + aligned-bases floor. Look at whether a **score-margin /
  ambiguity** criterion (best vs second-best across *references*, which MAPQ already computes)
  discriminates better than a fixed identity floor on a marker DB where near-identical markers
  from sibling species are the whole problem.
- Verify precision is not being spent on **multi-mapping reads emitted at an arbitrary one of
  several equally good markers**. If so, the fix is emission policy, not the gate.
- Re-check that the MAPQ→precision curve is still monotone after any change (it was fixed in
  `PERFORMANCE_FIX.md` change #6 and is the cheapest regression canary).

### C. Speed and memory

- **`flexalign-paf` on protal takes 446 s wall / 273 s analysis vs `flexalign` (SAM) at 106 s /
  46.9 s on identical inputs.** Same pipeline, different output stage — that is a bug or a
  pathological output path, and it is the single largest unexplained cost in the suite. Find it.
- **Sharding is 5–10× slower in analysis** (protal s2: 498 s vs 46.9 s unsharded) for a 2.7×
  RSS win. Batched form recovers most of it (88 s). Understand where sharded analysis time goes;
  sharding should cost bandwidth, not 10×.
- RSS: 47.7 GB vs minibwa's 36.6 GB on protal, 6.0 GB vs 3.1 GB on bacteroides. The mmap work
  (`memory/flexalign-db-mmap.md`) already fixed load time; peak resident set is the remaining item.
- Keep the batch path's advantage: 38 s vs 106 s on protal is index-load amortisation and it is
  real. Consider making batch the documented default for multi-sample runs.

### D. Sharded-path parity

`PERFORMANCE_FIX.md` flags that the sharded rejoin emit (`src/shard/rejoin_align.rs`) has not
always carried the same acceptance gate as the unsharded pipeline. Confirm parity holds today —
sharded and unsharded must produce the *same* accuracy numbers, differing only in time and RSS.
Current protal numbers differ (4.182 vs 4.204 recall; 87.79 vs 87.28 precision), which means they
are **not** identical pipelines. Explain that difference or remove it.

## How to measure (and what "no regression" means)

```bash
# build (the harness does this too, and checks the flags it needs exist)
cargo +nightly build --release

cd ../flexalign_benchmark
pixi run check                       # prerequisites, scorer, dataset paths
pixi run bench bacteroides           # fast loop: 1 sample, ~5 s/tool
pixi run bench markersim             # recall-sensitive, reference-level truth
pixi run bench protal                # the decisive one, 10 samples
TOOLS=flexalign,flexalign-batch pixi run bench protal   # re-run only my tool
pixi run results protal && pixi run report protal       # re-score + render, no re-run
```

- **Build once, then `FLEXALIGN_NO_BUILD=1`.** Every tool's `build.sh` re-enters
  `scripts/build_flexalign.sh`, so a plain `pixi run bench` pays a cargo invocation per contender —
  measured at ~25 min of rebuilding to run a 38 s benchmark, and a failed first build additionally
  triggers a full libwfa2 (WFA2-lib C) recompile. The loop is only viable as:

  ```bash
  cargo +nightly build --release                  # once, in ../flexalign
  cd ../flexalign_benchmark
  FLEXALIGN_NO_BUILD=1 TOOLS=flexalign-batch pixi run bench markersim
  ```

  `FLEXALIGN_NO_BUILD=1` verifies the existing binary and compiles nothing. Never edit + rebuild
  while a bench is running: `bin/flexalign` is a symlink into `target/release`, so a mid-run
  rebuild swaps the binary under the measurement.
- Don't pipe a backgrounded bench through `tail` — output buffers until exit and the run is opaque
  while it matters. Redirect to a log and read that.
- `benchmark.conf` pins `FLEXALIGN_HOME=../flexalign` and `FLEXALIGN_REF=shards`. My work is on
  `Optimization`, which is currently the same commit. **Do not let `FLEXALIGN_PULL=1` reset the
  checkout out from under uncommitted work**, and if the branches diverge, point `FLEXALIGN_REF`
  at the branch actually under test.
- Timing hygiene, learned the hard way: a "hot" run whose working set exceeds RAM silently lies.
  **mincore-verify residency** before quoting a hot number, and quote cold and hot separately
  (`pixi run bench-cold`). See `memory/benchmark-cold-hot-cache.md`.
- A results row records the *path* of a tool script, not its *contents* — check the artifacts
  (the actual `.sam`/`.paf`/`.time` files), not the config, when comparing across runs. See
  `memory/benchmark-settings-consistency.md`.

### Per-iteration protocol (2026-08-04 directive)

Every iteration runs all three, in this order — cheapest signal first:

| fixture | what it answers | cost |
| --- | --- | --- |
| `hardset-protal`, `hardset-markersim` — **NOT BUILT YET**: the extractor exists, the datasets.tsv rows and `data/hardset/` do not, so this line is a plan, not a runnable step | **did the change rescue what it targeted?** Only reads where flexalign-batch failed *and* a peer succeeded — built by `scripts/make_hardset.py` from benchpro's `--per-read` table | seconds |
| `realgut-1` | **speed**, real reads, 1 sample, no truth, no scoring | ~1 min |
| `markersim` | **regression**, full accuracy at reference-level truth | ~2 min + scoring |

The hardset is a **directional instrument, not a score**. Its recall starts near 0 by construction,
so every point gained is a genuinely rescued read — but a change can rescue hard reads and wreck
easy ones in the same edit. Never quote hardset recall as flexalign's recall, and never ship on
hardset evidence alone: the full-dataset regression check is what decides.

**Regression gate — a change ships only if, versus the baseline tables above:**

| dataset (contender) | recall% | precision% | wall s | RSS |
| --- | --- | --- | --- | --- |
| bacteroides (`flexalign`) | ≥ 98.0430 | ≥ 98.2762 | ≤ 4.83 (+5 % tol) | ≤ 6003 MB |
| protal (`flexalign-batch`) | ≥ 4.21380 | ≥ 86.9284 | ≤ 38.16 s (+5 % tol) | ≤ 48392 MB |
| protal (`flexalign`) | ≥ 4.20351 | ≥ 87.2785 | ≤ 106.11 s (+5 % tol) | ≤ 47709 MB |
| markersim (`flexalign-batch`) | ≥ 98.3410 | ≥ 98.7525 | ≤ 82.58 s (+5 % tol) | ≤ 48581 MB |

Thresholds are quoted at the artifact's own precision on purpose: rounding the protal-batch gate to
`≥ 4.214 / ≥ 86.93` lets a run at 4.2135 / 86.929 pass while being a real regression.

**bacteroides is gated on `flexalign`, not `-batch`** — `datasets.tsv` does not list `flexalign-batch`
among that dataset's contenders, so there is no batch row to compare. Since every variant resolves
anchors through the same code (scope directive), the unsharded row is the valid witness there.

plus: `cargo +nightly test` green, `cargo +nightly clippy --release` no new warnings, no output on
stdout other than alignments, and PAF and SAM outputs agreeing on placement for the same run.

A deliberate trade (recall up, precision down) is allowed **only** if the product improves on both
reference tools at once and is recorded explicitly in the results log with the numbers.

## Target (what "done" looks like)

| dataset | recall% | precision% | wall s |
| --- | --- | --- | --- |
| bacteroides | ≥ 98.28 (match minibwa) | ≥ 98.28 | ≤ 4.83 |
| protal | ≥ 4.40 (beat minibwa) | ≥ 89.2 (beat protal) | ≤ 38 (batch) |
| markersim | ≥ 98.76 (match minibwa) | ≥ 98.76 | ≤ 82.6 |

## Working rules

- **One change per benchmark run.** Bundled changes have already produced one wrong published
  result in this project's history.
- Keep a running log at the bottom of this file: change → dataset → before/after → keep/revert.
- Suspect the *index and the harness* before the algorithm. Both of the last two "algorithm bugs"
  in this repo were a stale index and a scorer artifact (`memory/flexalign-stale-index-regression.md`,
  `memory/benchmark-strobealign-unfiltered.md`).
- flexmap ideas go to `.claude/FLEXMAP.md`. Nothing in `../flexmap` gets edited.

## Results log

| date | change | dataset | before → after | verdict |
| --- | --- | --- | --- | --- |
| 08-04 | tie-extend the `--extend-top-z` cut (STRATEGY a) | protal, flexalign-batch | precision 86.9284 → 86.9673 (+0.039 pp); recall 4.21380 → 4.21573 (+0.002 pp); **+1 alignment of 95,589** | accuracy: negligible |
| 08-04 | same, timing only — paired A/B on one warm binary via `--extend-tie-cap` | protal, 10 samples, DB load 5.7–5.8 s (warm) | cap=1 (off): 37.84 / 37.69 s · cap=2: 39.10 / 39.27 s · cap=64: 39.16 / 38.95 s | **costs ~1.4 s (+3.7 % wall)** |

**Reading of the first two rows.** The change buys essentially nothing on protal (+1 alignment) and
costs 3.7 % wall. It is not justified by this dataset. Whether it survives depends on `markersim`,
where truth is resolved at REFERENCE level and a tie-break between two near-identical markers can
actually change a verdict — on protal's species-level scoring, tie-breaking between markers of the
same species is invisible by construction, so protal cannot see this change's upside even if it exists.

**Two false alarms this produced, both worth remembering.** A first protal run reported
38.16 s → 88.53 s, and a second → 70.5 s: both were **cold-cache artifacts**. The tell is
`db_load 5.6 s → 30.1 s` with "45.6 GB actually read from disk" in the `.time` file — a change to
anchor selection cannot touch DB load, so any run whose load moved is not comparable. The same
binary also produced analysis times of 58.42 s and 42.91 s on two consecutive runs, i.e. run-to-run
noise larger than the effect. The instrument that worked: **one binary, two flag settings
(`--extend-tie-cap 1` vs `2`), alternated, with DB load verified warm** — `cap 1` reproduces the
pre-change code path exactly, so no rebuild and no cache difference sits between the arms.
