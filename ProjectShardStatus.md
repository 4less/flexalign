# ProjectShard — implementation status

Tracks what is built against the design in [ProjectShard.md](ProjectShard.md). Section numbers
(§) refer to that document. Branch: `shards`. Everything here is **additive** — see the constraint
note at the bottom.

Last updated: 2026-07-17.

---

## Summary

The **format** and the **cut** are built and unit-tested: a shard pass could produce evidence and
the rejoin could rebuild seeds from it, and the walk that turns ranges into seeds now exists once
and is shared by both the sharded and unsharded paths. What is missing is the **plumbing that runs
them** — the per-pass drivers, the physical index slice, and the read-id support in bioreader that
lets evidence re-pair with reads. Nothing is wired into a runnable subcommand yet.

37 shard unit tests pass; the library builds warning-free. One unrelated pre-existing test failure
(`distance::indel_detection`) is documented below.

---

## Built

| Module | What it is | § | Origin |
|---|---|---|---|
| [`src/shard/record.rs`](src/shard/record.rs) | Wire format: `GroupWriter`/`GroupReader`, `RangeRecord`, `GroupRecord`, LEB128 varints, monotone cursor | §7 | pre-existing (bug fixed, see below) |
| [`src/shard/chunkfile.rs`](src/shard/chunkfile.rs) | Self-indexed blob container; atomic-bump concurrent writes; `read_at` positioned reads | §8 | pre-existing |
| [`src/shard/db.rs`](src/shard/db.rs) | `ShardedDB<D>` — key-range view over a built index; `ShardRange::split_even` | §2 | new |
| [`src/shard/manifest.rs`](src/shard/manifest.rs) | Run manifest + `verify_compatible` (the hurdle-6 assertion) | §6, §8 | new |
| [`src/shard/emit.rs`](src/shard/emit.rs) | `emit_ranges`: `Range` → `RangeRecord`; applies flex filter, drops discarded | §4, §7 | new |
| [`src/shard/replay.rs`](src/shard/replay.rs) | `ReplaySeedExtractor`: `RangeRecord` → `[Seed]`, the rejoin's seed rebuild | §9 | new |

### The shared walk (§13.4, hurdle 11)

The single most important structural change. The `retrieve_seeds` walk was factored out of
`StdSeedExtractor` into [`seed_extractor.rs`](src/align/process/seed_extractor.rs):

- `walk_ranges<K,C,F,R: SeedRange>(...)` — the budget/cutoff loop, written once.
- `trait SeedRange` — a range that can emit its seeds and report its verdict
  (`Headered` / `Headerless` / `Discarded`).
- `impl SeedRange for Range` (unsharded, live map) and `impl SeedRange for RangeRecord`
  (rejoin, off the wire) — the only thing that differs between the two modes.

This is **behaviour-preserving** for the unsharded pipeline: same counters, same `break`
placement, same header re-scan. The design flags a duplicated walk as *the main correctness risk*,
because with exactness relaxed (§4) nothing would catch the two copies drifting.

### Bug fixed in `record.rs`

`min_dist` is a hamming distance over `F=16` flanking bases, so it spans `0..=16` — 17 values,
needing **5 bits**. The constants encoded 4 bits (`DIST_SHIFT=4, DIST_MAX=0x0F`) and *rejected*
`min_dist == 16`, so an all-mismatch flank would have errored out mid-pass. `min_dist` lands in
`Seed::mismatch` and feeds anchor scoring, so 16 must not alias to 0. Fixed to
`DIST_SHIFT=3, DIST_MAX=0x1F`; the test that had frozen the bug (asserting 16 is rejected) was
replaced with round-trip coverage over `0..=16` plus a flags-collision test.

---

## Not built yet

Ordered roughly by dependency. Items 1–3 are the design's "do this regardless" measurement steps;
4 is the upstream dependency everything runtime-facing needs.

1. **Measurement (§11, §13.1–3).** Histogram `r` (ranges/mate), `positions.len()`, and
   `count_at_min_dist` on a real sample; measure what fraction of reads trigger the recovery path
   (decides §4 vs §4a); check whether the current aligner is already gzip-bound (§10a). This is the
   go/no-go and it needs no sharding code. Runnable on `main` today.

2. **Fix the recovery-path seed duplication on `main` (§3.4).** `retrieve_seeds` pushes without
   clearing, so the `--min-ranges` recovery walk re-emits every seed the first walk produced. It
   pollutes the baseline the sharded run will be measured against. *(Note: this is a `main` fix,
   outside the additive-shard scope — track separately.)*

3. **Flex-threshold sweep (§13.3).** "No recovery, `--max-best-flex ∈ {16,24,32}`" vs unsharded
   via the existing `GOLDSTD_EVAL` harness. Settles §4 before any pass driver exists.

4. **bioreader read-id plumbing (§6).** ~30-line upstream change: a batch counter on
   `FastqPairedByteReader`, surfaced through `PairedFastqReader`, plus a
   `read_fastq_paired_end_state_par` variant whose closure receives `(rec1, rec2, read_id, state)`.
   `read_id = (batch_id, idx_in_batch)`. **Everything below depends on this** — it is how evidence
   re-pairs with reads. Prove determinism with the 1-vs-32-thread identical-map test (§6, hurdle 7).

5. **Physical index slice + `flexalign shard` subcommand (§2, §13.6).** Slice an existing
   `.flex.index` into N shard files + manifest via the control-head prefix-sum binary search (for
   byte balance, vs the key-count `split_even` we have now). Verify
   `∀k: sharded.get_vrange(k) == unsharded.get_vrange(k)`.

6. **Shard pass driver (§13.7).** Run `StdKmerExtractor` + `StdRangeExtractor` against a
   `ShardedDB`, feed each mate's ranges through `emit_ranges`, group per `(read, mate)`, and write
   via `GroupWriter` → `ChunkFileWriter`. One chunk per bioreader batch. Report only per-pass IO
   stats, not the pipeline `Stats` (hurdle 9).

7. **Rejoin driver (§13.8).** Drive the existing parallel FASTQ reader; per batch, `read_at` the N
   shard blobs, walk reads in `read_id` order with N monotone `GroupReader` cursors, feed
   `ReplaySeedExtractor`, then continue into today's `ModularPE` from `AnchorExtractor` onward
   unchanged. Build the **in-RAM `EvidenceSink` path first** (§10) — no intermediate files, easier
   to debug.

8. **Acceptance (§13.9).** `GOLDSTD_EVAL` sensitivity/precision of sharded vs unsharded on a
   simulated sample, plus a seed-level diff at the cut line confirming the only differences are the
   ones §4 predicts. With exactness relaxed, this comparison *is* the spec.

9. **Scale-out (§13.10).** Disk sink, parallel rejoin, multi-sample shards-outer.

### Deferred by design, not forgotten

- **§4a deferred-read fixup pass** — only if step-1 measurement shows the recovery path is
  load-bearing. The format already tolerates it (`n_cells = 0` header-only records round-trip).
- **§7 optional wins** — per-range cell delta-varints (~30%), per-batch zstd-1. Both need
  measurement first.
- **§10a syncmer partition pass 0** — the structurally right fix for the N× inflate cost, but a
  second record format; build only if inflate + syncmer extraction actually dominate.

---

## Test status

- **37 shard unit tests pass**; `cargo build --lib` is warning-free for the new code.
- `distance::indel_detection::tests::test_offset_similarities` **fails**, but identically on a
  clean checkout with none of the shard work present (verified in a throwaway worktree). It is
  pre-existing WIP breakage in a module unrelated to sharding, not a regression from this work.

## Constraint

The shard feature is **additive**: it must not change the behaviour of existing code. The only
edits to pre-existing files are `pub mod shard;` in `src/lib.rs`, the additive `serde_json`
dependency, and the behaviour-preserving `walk_ranges` refactor (approved specifically to avoid
duplicating the walk). Note the `shards` branch *also* carries an unrelated savefile→bincode
migration in `src/database/*` and `src/misc.rs` (commit `3bdf36c`) — that is a **different
workstream** sharing the working tree, not part of ProjectShard, and it *does* alter existing
save/load behaviour. Consider committing the two separately.
