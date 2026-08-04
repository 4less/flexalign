# Night review 10 — tie-extension cap, `--extend-tie-cap`, make_hardset rewrite, TASK.md baselines

_night-reviewer, 2026-08-04. Verdict: **STOP**. Saved verbatim; disposition at the bottom._

## Central question: are the baselines stale? — YES, and the cause is bigger than I thought

`--extend-tie-cap 1` IS a bit-exact no-op in release (`tie_cap == screen_z`, loop guard
`screen_z < screen_z` false; `debug_assert!` compiled out since `[profile.release]` does not set
`debug-assertions`; funnel changes are behind the `GOLDSTD_EVAL` compile-time const). So the
bacteroides delta cannot come from the change.

Eliminated: rebuilt index (unchanged since 07-27), different reads/truth (07-19/07-21), changed tool
script (no diff), `--mate-rescue` (`Mate rescue attempted 0`), thread nondeterminism (both caps emit
exactly 1 987 115), benchpro drift (same binary, same session).

**The real cause is commit vintage PLUS thread count.** The baseline artifact
`results/bacteroides/single/flexalign-paf.time` is dated 2026-07-28 09:39 and its command line ends
`... 8 scripts/tools/flexalign/params` — **`-t 8`**. Tonight's run is `-t 16`. The binaries differ by
inspection (old prints `Memory-mapping index took`, new prints the `Database load` block and a
mate-level funnel). Five commits sit between, and `c92c8a8` is dispositive:

```
-    #[arg(long = "min-query-coverage", default_value_t = 0.70)]
+    #[arg(long = "min-query-coverage", default_value_t = 35)]
```

That gate change alone explains 1 986 814 → 1 987 115 and the ~0.026 pp precision drop that goes
with emitting 301 more marginal records.

## CRITICAL

**C1. `options.rs:84,90` state a measured bacteroides regression that my own artifacts refute — and
it is the entire justification for `default_value_t = 1`.**

```
cap=2 (23:07): 99.7778  98.2499  98.0316   1987115 records
cap=1 (23:11): 99.7778  98.2498  98.0315   1987115 records
```

The tie run on bacteroides is **+0.0001 pp / +0.0001 pp / zero record delta** — a wash, marginally
positive. The documented `-0.0114 / -0.0263 / "301 more records placed worse"` is cap=1 minus the
**July-28 stale baseline**, i.e. the delta of commits `14aefc7…c92c8a8`, attributed in shipped
`--help` text to a flag that demonstrably did not cause it. The paralog/repeat mechanism paragraph
is an explanation invented for a non-existent effect.

**C2. The cap costs accuracy on protal and buys no measurable time.**

| protal, flexalign-batch | precision | recall | wall (2 reps) |
| --- | --- | --- | --- |
| cap=1 | 86.9284* | 4.21380* | 37.84 / 37.69 |
| cap=2 | 86.9557 | 4.21517 | 39.10 / 39.27 |
| cap=64 (≈uncapped) | **86.9673**† | **4.21573**† | 39.16 / 38.95 |

`cap=64` is not slower than `cap=2` yet is 0.0116 pp more precise. "Uncapped, one such read can cost
more than a thousand ordinary ones" is unsupported: cap=64 adds **4** successful alignments over
cap=2 (165 218 → 165 222), not thousands.

**C3. TASK.md's baselines and gate are stale for every flexalign row, and the gate now fails HEAD
against itself.** bacteroides rows are 07-28 artifacts at `-t 8`; protal unsharded/sharded and
markersim `flexalign` are 08-02, pre-`c92c8a8`. The gate demands bacteroides
`≥ 98.0430 / ≥ 98.2762 / ≤ 4.83 s / ≤ 6003 MB`; HEAD with the change disabled measures
`98.0315 / 98.2498 / 5.75 s / 6514 MB` — a four-way "regression" against a July-28, 8-thread binary.

**C4. The results log records a code version that no longer exists** (the `86.9673` arm is the
UNCAPPED run; shipped capped code gives `86.9557 / 4.21517`) **and omits three runs that were made**,
including the cap2-vs-cap1 bacteroides pair that refutes C1.

**C5. `results/{protal,bacteroides}/summary.tsv` are internally inconsistent** — one 16-thread 08-04
row beside seven 8-thread 07-28 rows. Read as a leaderboard they give wrong answers.

## WARNING

**W1.** `modular_workflow.rs:1060` "The window still ends at the first worse anchor" is false when
the cap binds — it then ends mid-tie-run, restoring the arrival-order coin flip at 2Z instead of Z.

**W2.** Funnel counters now measure a window the pipeline does not use (`screen_z` counted,
`extend_top_z` extended), so with cap>1 they report "truth lost at the screen" for reads whose truth
was screened. Keep both: `funnel_true_in_screen` on `screen_z` AND `funnel_true_in_extended`.

**W3.** The markersim delta in `options.rs:82` is an across-binary comparison (08-03 vs 08-04, with
`c92c8a8` in between), not the one-binary A/B documented as the working instrument.

**W4.** Disk at 99 %, 23 GB free, holding 29 GB of scratch SAMs; a truncated SAM scores silently.
Orphaned `until ! pgrep …` polling shells from killed sessions are still running.

**W5.** `make_hardset.py` failure paths leave the bad fixture on disk: the desync/size guards fire
*after* both FASTQs and the truth file are written to their final names.

**W6.** `_Writer` is not a context manager while `_Reader` is, so an exception in the read loop skips
`out.close()`, leaving an unflushed pigz child and a partial file.

**W7.** `cargo clippy --release` was never run this session although TASK.md requires it; and tests
ran only `--release`, where the new `debug_assert!` is compiled out.

## NOTE

**N1.** The review-5 fixes in `make_hardset.py` are all correct (per-mate keying, desync guards,
pigz status on both sides, peer-less detection, documented `single`→`1`).
**N2.** `_Reader` on a plain uncompressed file works; review-5 N5 resolved.
**N3.** `_Reader.__exit__` closes the pipe before `wait()`; latent SIGPIPE misreport on early break.
**N4.** `extract_truth` checks `nt != 0`, not `nt == 2 * len(ids)`.
**N5.** `--max-pairs` truncates lexicographically — a biased prefix, not a sample.
**N6.** The table's `dataset` column is still ignored (latent).
**N7.** `pair_anchor_score` is the third copy of the same expression; promote to
`AnchorScore::score_paired`.
**N8.** The `debug_assert!` scans the whole slice after the loop; cheaper over `anchors[..tie_cap]`
before it.
**N9.** `--extend-tie-cap 0` silently means 1.
**N10.** The memory file omits the night's biggest lesson: a `summary.tsv` row can be a month old and
run at a different `-t`, so the baseline must be re-derived on current HEAD before any A/B.

---

## Disposition

| id | action |
| --- | --- |
| C1 | **Accepted in full.** The false bacteroides deltas and the invented mechanism are deleted from `options.rs`; replaced with the real paired cap2-vs-cap1 numbers. The default flip to 1 rested on this and is now **unjustified** — held at 1 pending the two paired accuracy runs (markersim scoring in flight, protal cap=1 to run), because flipping a default twice on unsound evidence is how this happened. |
| C2 | Accepted. The "thousand ordinary reads" claim is deleted from both `modular_workflow.rs` and `options.rs`; the cap is described as an unmeasured safety belt. Whether to drop it entirely is decided with C1's numbers. |
| C3 | Accepted. Every flexalign baseline row in TASK.md is marked INVALID with its artifact date and thread count; the gate is marked unusable until re-derived on HEAD. |
| C4 | Accepted — results log corrected and the missing runs added. |
| C5 | Accepted — recorded in TASK.md as a property of those files. |
| W1 | Fixed — the invariant sentence now states the cap case. |
| W2 | Fixed — added `funnel_true_in_extended` alongside the screen counter. |
| W3 | Accepted — the markersim line is marked as corroboration until the in-flight paired scoring lands. |
| W4 | Fixed — scratch SAMs deleted after scoring; orphaned pollers killed. |
| W5/W6 | Fixed — write to `.tmp` and `os.replace` only after every check passes; `_Writer` is a context manager. |
| W7 | Fixed — clippy run; debug-profile test run added. |
| N1–N10 | N7/N8/N9/N10 actioned; N3/N4/N5/N6 recorded as latent. |
