# TODO — syncmer selection (s=4) rework

Diagnosis and fixes for the core-mer syncmer selector. All items change **which cores get
indexed**, so they are global result changes: they must land in **both** `flexalign` and the
`flexmap` crate, they invalidate every existing index (rebuild required), and they should **not**
ride along on the `shards` branch (that work is meant to be additive / behaviour-preserving).

## Context you need before touching anything

- Params: `K=31`, `C=15`, `F=16`, `S=4`, `L = C−S+1 = 12` — [src/flexalign.rs:36-40](src/flexalign.rs#L36-L40).
  `L` is derived, so changing `S` updates `L` automatically.
- The selector: `ClosedSyncmer::is_minimizer` in the `kmerr` crate,
  `~/.cargo/git/checkouts/kmerr-2a352a43fd59c860/ed923ec/src/syncmer/closed_syncmer.rs:40-48`.
  (Source repo: https://github.com/4less/kmerr.git — fix it there, not in the cargo checkout.)
- **Build/query symmetry is load-bearing.** Query side selects in
  [src/align/process/kmer_extractor.rs:41](src/align/process/kmer_extractor.rs#L41); build side runs the
  *identical* `ClosedSyncmer<C,S,L>` in `flexmap` `src/build.rs` (canonicalization at lines 100-104 /
  162-166 / 258-261 / 360-364). Any change to the selection rule MUST be made once, in the shared
  `kmerr` crate, so both sides move together. After changing, confirm build and query still call the
  same code path.

## Important correction to the mental model

The deployed selector is **NOT a closed syncmer**. Line 47 reads:

```rust
min == &self.smers[0]// && hash < 33554432// || min == self.smers.last().unwrap()
```

The `|| min == self.smers.last()` is commented out. What ships is an **open syncmer with offset
t=0**: selected iff the minimum s-mer value sits at position 0. Nominal density is therefore
`1/(C−s+1) = 1/12`, not `2/12`. Keep this straight when quoting density targets below.

---

## Item 1 — Hash the s-mer before comparing (fixes A-rich compositional bias)

**Problem.** `load()` packs raw 2-bit s-mers (closed_syncmer.rs:30) and they are compared by raw
`<`. Since `A = 0b00`, `AAAA = 0x00` is the global lexicographic minimum, so any core whose first
4-mer is AAAA is unconditionally selected, with a systematic bias toward A-rich cores.

**Fix.** Apply a bijection on the s-mer value before the min/compare (a fixed permutation / integer
hash over the `2*S`-bit space; must be deterministic and identical build-side and query-side). A
bijection removes AAAA's privilege without changing the value set.

**Where.** `kmerr` `closed_syncmer.rs`: hash inside `load()` (or in `is_minimizer` before the min),
so every comparison is over hashed values.

**Note.** This does NOT reduce tie probability (a bijection maps equal inputs to equal outputs), and
does NOT touch the homopolymer case (Item 3). Necessary but not sufficient.

## Item 2 — Raise s (fixes the tie rate)

**Problem.** s=4 → only 256 distinct s-mer values across 12 overlapping windows. Ties are frequent
(~20% of cores have a tie somewhere), which skews selection density off nominal.

**Fix.** Set `S = 5` or `S = 6` in [src/flexalign.rs:39](src/flexalign.rs#L39). `L` updates
automatically (line 40).
- `S=5` → 11 windows, 1024 values, P(tie) ≈ 5%, open density 1/11.
- `S=6` → 10 windows, 4096 values, P(tie) ≈ 1%, open density 1/10. Comfortably clean.

**Note.** Density figures are the **open** form (1/L). Only use 2/L if you also re-enable the
`last()` branch (Item 4). The inline comment on flexalign.rs:39 records an empirical sweep
(`//6 0.31 //5 0.29 //4 0.289 …`) — re-measure that metric after Items 1 & 3, since the old sweep
was done with the bias and ties present.

## Item 3 — Decide tie-breaking explicitly (fixes homopolymer / low-complexity over-indexing)

**Problem — this is the big one, and Items 1 & 2 do NOT fix it.** The check `min == smers[0]` is a
**value** equality, not an argmin. In a homopolymer / low-complexity tract all 12 s-mers are equal,
so `smers[0]` trivially holds the min and the core is selected — density → 1.0, in exactly the
regions that already produce the largest value blocks (poly-A, alpha satellite, simple repeats get
indexed at up to ~12× the intended rate). Raising s and hashing don't help: equal s-mers stay equal.

**Fix.** Make the tie-break a deliberate decision instead of the current implicit "position 0 counts
if it ties the min" (the most permissive option in low-complexity). Options:
- **Strict / last-occurrence tie-break** so a run does not select at every position, or
- an explicit **low-complexity mask** (skip cores below an entropy/complexity threshold).

**Where.** `kmerr` `closed_syncmer.rs:45-47` for the tie-break; a mask could live at the call site
in [kmer_extractor.rs:41](src/align/process/kmer_extractor.rs#L41) (and the matching build-side loop
in `flexmap` `src/build.rs`).

**Priority.** If only one item ships, ship this one — it targets the largest-block problem directly.

## Item 4 — Closed vs open: make it an explicit choice

**Problem.** The struct is named `ClosedSyncmer` but the `|| min == self.smers.last()` branch is
commented out (closed_syncmer.rs:47), so it runs as an open syncmer at t=0. This is silent — decide
it on purpose.

**Fix.** Either (a) re-enable the closed branch (density → 2/L, more/denser seeds), or (b) keep it
open but consider moving the offset to the middle `t=(C−s)/2` for better mutation robustness.
- Note (b) does **not** lower density vs current (already open at 1/L); it only relocates the biased
  position.
- `t=(C−s)/2` is integer only when `C−s` is even → pairs cleanly with **s=5** (t=5, center of 11).
  For s=4 (L=12) there is no single center.
- Conservation benefit of the middle offset is weak here: a core SNP kills the seed regardless of
  where t sits.

## Item 5 — Micro-optimization (do last, depends on Item 4)

**Current.** `is_minimizer` uses `.iter().min()` — a value reduction over 12 lanes (not argmin;
`index_min()` exists but is not on the hot path).

**Fix.** Once the selection rule is fixed:
- Open t=0 rule: test `s[0] <= min(s[1..])` with early-exit the moment any interior lane beats
  `s[0]`.
- Closed rule: `min(s[0], s[L-1]) <= min(s[1..L-1])`.

Both avoid a full scan in the common reject case. Pick the form matching Item 4.

---

## Suggested order & validation

1. Land Items 1 + 3 in `kmerr` (bias + homopolymer). 2. Set `S` (Item 2). 3. Decide Item 4.
4. Micro-opt (Item 5) last.

**Validate after each change:**
- Build side and query side compile against the same `kmerr` selector (grep both for
  `ClosedSyncmer<`).
- **Rebuild all indexes** — old databases are invalid.
- Measure selection **density in low-complexity tracts** (poly-A / satellite) and overall
  sensitivity before/after. Re-do the flexalign.rs:39 empirical sweep under the new selector.
- Do this on a dedicated branch, off `shards`.

---

# Part B — lookup semantics & database layout

Separate body of work from Part A (that one changes *which* cores are selected; this one changes
*how* matches are retrieved and *how* the index is stored). Same guardrails apply: result-affecting
changes (B1–B2) invalidate indexes and belong off `shards`; the storage items (B3–B7) are
memory/perf and should be **gated on measurement (B-GATE)**.

Part A items 1–3 of the reviewer's list are already covered above:
- "raise s + hash + re-measure density" = Items 1 + 2.
- "low-complexity rejection on the core" = the mask option in Item 3.
- "explicit tie-break + `min(ends)<=min(interior)`" = Items 3 + 5.

## B1 — Return all entries within a distance threshold, not the argmin

**Problem / opportunity.** The flex filter today returns only cells that tie at the *exact* minimum
flank distance: `min_flank_dist` then `cells_at_dist(min_dist)` —
[src/align/process/seed_extractor.rs:84-114](src/align/process/seed_extractor.rs#L84-L114). Exact-min
can miss a true hit whose flank carries a mismatch, costing sensitivity.

**Fix.** Generalize exact-min to "within a distance threshold `d`": return every cell whose flank
distance `<= min_dist + d` (or `<= d` absolute — decide which). This is an evolution of the existing
filter, not new machinery. Interacts with `max_best_flex` (the count budget in
`emit_seeds`/`retrieve_seeds`) — widening the distance widens the candidate set, so re-tune that
budget together.

**Caveat.** This is a sensitivity/specificity tradeoff and a **result change** — only meaningful with
B2's calibration. Measure false-seed rate before/after.

## B2 — Make the threshold a function of block size

**Problem.** A fixed distance threshold admits far more chance matches in a large block than a small
one (more flanks = more that fall within `d` by luck).

**Fix.** Scale the B1 threshold by block size, calibrated so the **expected number of chance hits
stays below ~0.1** per lookup. Smaller (tighter) threshold for large blocks, looser for rare/small
blocks. Depends on B1 and on the block-size distribution from B-GATE.

## B-GATE — Instrument the real reference BEFORE implementing anything below

**Do this first.** Measure, on the actual reference:
- the empirical distribution of **block sizes**, and
- **distinct-flex-mers-per-block**.

Everything below (B3–B7), and the calibration in B2, depends on these numbers. In particular it
tells you whether the singleton-majority assumption (B4/B6) and the sparsity assumption (B3) actually
hold.

## B3 — Replace the direct-indexed key table with rank/select + dense offsets

**Problem.** `FMKeys` is direct-indexed over the whole core space: `kmer_to_index` →
`data[4^C ...]` ([flexmap keys.rs:53-57, 97-98](file)). At C=15 that's ~2^30 `KCell`s (~2 GB)
regardless of how many cores are actually occupied.

**Fix.** Replace the dense table with a **rank/select bitvector over occupied cores + a dense offset
array**, so storage scales with occupancy, not with `4^C`.

**GATE.** *Value is contingent on B-GATE occupancy.* If the reference nearly saturates 2^30 distinct
cores, the dense table is already near-optimal and this saves little; if occupancy is sparse it is a
large win. **Measure occupancy (B-GATE) before building this** — do not treat it as a free win.

## B4 — Sort + run-length-encode headers for large blocks; parallel arrays for singletons

**Fix.** For blocks above a size threshold, **sort and RLE the header** (`HeaderSeq` flank codes,
[flexmap values.rs:19,108-110](file)). Keep the singleton majority in separate parallel arrays so the
common case stays flat and cheap. Depends on the block-size skew from B-GATE.

## B5 — Cap run lengths but keep a true-count flag for MAPQ

**Problem.** RLE / capping in B4 loses repeat multiplicity, which MAPQ needs to down-weight repeats.

**Fix.** Cap run lengths for storage but store a **true-count flag / field** so MAPQ can still see the
real multiplicity of a repeat. Pairs with B4.

## B6 — Drop inline flex for singleton blocks; verify from the 2-bit reference

**Fix.** For singleton blocks (one position), don't store the flex-mer inline — **recompute/verify
the flanks from the 2-bit-packed reference** at that position. Saves memory across the singleton
majority. Requires the 2-bit reference to be retained with cheap random access. Depends on B-GATE
confirming singletons dominate.

## B7 — Batched prefetching across the key → header → body chain

**Problem.** Lookup chases pointers: key table → header → body (`VRange.positions`), a
cache-hostile pattern over the large tables above.

**Fix.** Process lookups in batches with **software prefetch** across the key→header→body chain to
hide memory latency. Structurally independent, but best applied after the layout (B3/B4/B6) is
settled so you prefetch the final structure.

## Part B order

1. **B-GATE first** — nothing below it is worth building unmeasured.
2. B1 + B2 (lookup semantics; result change; re-measure sensitivity + false-seed rate).
3. B3 (only if B-GATE shows sparse occupancy).
4. B4 + B5 (large-block compression + MAPQ multiplicity), then B6 (singleton flex).
5. B7 (prefetch) last, over the settled layout.
