# ProjectShard — key-sharded flexmap with deferred rejoin

Status: design draft. Nothing implemented. Numbers marked `~` are estimates derived from the
current constants; formulas are given so they can be re-derived from a real index.

Scope: shard the flexmap **by key** (core-mer), run N prefilter passes emitting per-read
evidence, then rejoin the evidence and run the reference-based stages once. Reference memory is
out of scope.

**Exactness budget (decided):** sharded output does *not* need to be bit-identical to unsharded.
Same-or-better sensitivity, re-tuned, is acceptable. §4 spends that budget in exactly one place
and quantifies what it buys; everything else stays exact because staying exact is free there.

**Not negotiable:** read-id determinism (§6). That is a correctness property — evidence must
re-pair with the right read — not a heuristic. Do not let the relaxed sensitivity budget leak
into it.

The hard parts are §4, §7 (record format) and §9 (rejoin). Everything before that is groundwork.

---

## 1. Where the pipeline splits

`ModularPE::run` ([src/align/modular_workflow.rs](src/align/modular_workflow.rs)) is already
decomposed into traits ([src/align/common.rs](src/align/common.rs)):

```
KmerExtractor    read  -> [(qpos, Kmer<31>)]              no DB access
RangeExtractor   kmer  -> VRange (db.get_vrange(cmer))    map only     <-- shard pass
SeedExtractor    VRange -> [Seed]  sorted by (rval,rpos)  map only     <-- shard pass
=========================== cut =====================================
AnchorExtractor  [Seed] -> [AnchorPair]                   no DB access
AnchorExtender   hamming vs reference                     references only
Align (WFA2)     top-y anchors                            references only
Output           PAF/SAM
```

`modular_workflow.rs:563` already marks this line in a comment. The cut is real and it is free:
a shard pass reuses `StdKmerExtractor` + `StdRangeExtractor` unchanged and replaces
`StdSeedExtractor` with an emitter; the rejoin replaces `StdSeedExtractor` with a replayer and
leaves everything below the cut untouched.

Parameters (`flexalign.rs`, `options.rs`):
`K=31, C=15 (the key), F=16, S=4, CELLS_PER_BODY=16, HEADER_THRESHOLD=2`;
`--ranges 15`, `--max-range-size 256`, `--max-best-flex 16`, `--min-ranges 4`.

## 2. Slicing the map is a post-processing step

`FMKeys` is **not** a hash table — it is direct-addressed over all `2^30` core-mers, in key order:

```
kmer_to_ctrl_block_index(k) = (k >> 4) * 20      // 4 head cells + 16 body cells per 16 keys
kmer_to_index(k)            = (k >> 4) * 20 + 4 + (k & 15)
```

- Key table is **fixed at ~2.68 GB** (`2^30 + 2^28 + 4` cells × 2 B); it does not grow with the DB.
- Values are ordered by key; each key's headers and positions are adjacent (`FMValues::get_range`).
- The control head is an **absolute** u64 offset into values; the per-key `KCell` is a
  **block-relative** u16.

So a shard is a key range `[lo, hi)` (multiples of 16) = a contiguous slice of *both* arrays:

```
keys   : [(lo>>4)*20, (hi>>4)*20 + 4)   // +4 = trailing ctrl head, read as the end sentinel
values : [ctrl_head(lo), ctrl_head(hi))
fixup  : subtract ctrl_head(lo) from each ctrl head in the slice; KCells untouched
```

The ctrl-head array **is already a prefix sum over value offsets**, so choosing N byte-balanced
boundaries is a binary search, not a scan. Cost per block: `20/16 * 2 B per key + 8 B per value cell`.

**`build.rs` does not change.** A `flexalign shard --index X --splits N` subcommand reads an
existing `.flex.index` and writes N shard files plus a manifest. `ShardedDB` implementing
`FlexalignDatabase` is ~20 lines:
`get_vrange(cmer) = if (lo..hi).contains(&cmer) { inner.get_vrange(cmer) } else { None }`.

Key-sharding is **bit-exact per key**: a key lives entirely in one shard, so `positions.len()`,
`HEADER_THRESHOLD` and the flank headers all see exactly what they see today. Every semantic
question in this design is therefore about *cross-key* decisions only.

## 3. The contract the rejoin must reproduce

The rejoin hands `StdPairedAnchorExtractor::generate` what `StdSeedExtractor::generate` hands it
today: **`&[Seed]` per mate, sorted by `(rval, rpos)`** via `glidesort`. Everything below the cut
runs unmodified.

Four details in [src/align/process/seed_extractor.rs](src/align/process/seed_extractor.rs) that a
naive merge gets wrong:

1. **Ranges are ordered by `positions.len()` ascending** (rarest first), in
   `StdRangeExtractor::generate`. For a headered block that is `size - header_size`, not the raw
   block size. This ordering is **global** across a read's ranges — the one real cross-key
   dependency, and the reason §7 emits at range granularity.

2. **Headerless ranges do not consume the `--ranges` budget.** In `retrieve_seeds`,
   `retrieved_matches += 1` happens **only in the `Some(headers)` arm**. The `None` arm
   (≤ `HEADER_THRESHOLD` positions) pushes seeds and increments nothing. The `break` never fires
   on them, but they *are* skipped if they sort after the range that trips it.

3. **`sort_unstable_by_key` is not stable**, so its tie-break depends on the pre-sort order.
   Reproducing that exactly would force the rejoin to re-sort by qpos first to reconstruct kmer
   order. Since exactness is not required: **sort by `(n_positions, qpos)`** — one sort, an
   explicit tie-break, and deterministic, which the current code arguably is not across
   refactors. This is a strict improvement and costs nothing.

4. **The `--min-ranges` recovery path appends to a non-cleared buffer.** `generate` calls
   `self.seeds.clear()` once, then `retrieve_seeds(..., max_best_flex, ...)`, then — if
   `retrieved < min_ranges && discarded > 0` — `retrieve_seeds(..., 128, ...)` again.
   `retrieve_seeds` pushes without clearing, so the second walk **re-emits every seed the first
   walk already emitted**. Recovery reads currently feed duplicate seeds to the anchor extractor.
   This looks like a bug. It matters here mainly because the unsharded run is the **baseline you
   will measure the sharded run against** — fix it on main first, or you are tuning against a
   polluted control.

## 4. Spending the exactness budget: drop the recovery path

This is the whole relaxation, and it is worth being precise about why it is the right place to
spend it.

The flex filter (`count_at_min_dist <= max_best_flex`) is key-local in its *test* but not in its
*threshold*: the recovery path retroactively raises 16 → 128 based on a **global** condition
(`retrieved < min_ranges && discarded > 0`). Staying exact therefore forces shards to emit the
min-dist set capped at **128** for every range, plus `count`, plus zero-cell records for
discarded ranges so the rejoin can compute the tally.

**That is the expensive part of the format, and it is expensive in the worst place.** A range with
`16 < count <= 128` carries up to 128 cells ≈ 1 KB — against a ~21 B typical range. These are
exactly the repetitive k-mers, and their payload is emitted on every read that touches them, to be
discarded unread except on the small fraction of reads that trigger recovery.

Dropping it:

- Shards apply `count <= max_best_flex` (16) themselves and **emit nothing** for ranges that fail.
- `count_at_min` leaves the wire. Zero-cell discarded records leave the wire. The `discarded`
  tally disappears from the rejoin.
- Plausibly **20–50 % off the intermediate** — needs measurement (§11), but it is the only field
  with a 50× tail.

The sacrifice: reads whose ranges are dominated by flex-discarded repeats lose their fallback and
get less evidence. Compensate by **re-tuning `--max-best-flex`** (16 → 24/32), which buys back
sensitivity at a bounded, uniform cost instead of a rare 1 KB spike.

**You can measure this directly.** The `GOLDSTD_EVAL` harness already exists
(`stats.gold_std_evaluation`, `Stats::plot_mapq`, the `.fp` read output in
[process_fastq.rs](src/align/process_fastq.rs)), so the sensitivity delta of "no recovery, flex=X"
is a parameter sweep on a simulated sample against unsharded — **runnable on main today, before
any sharding code exists**.

Do this first (5 lines of instrumentation): **what fraction of reads actually trigger recovery?**
If it is ~0, delete the path from sharded mode and stop reading this section. If it carries real
sensitivity, the escape hatch is §4a rather than reinstating the 128-cap emission.

### 4a. Escape hatch: deferred-read fixup pass

If recovery turns out to matter, it can be made exact for near-zero steady-state cost, because
**the rejoin has the read in hand at the moment it discovers the need** (it is streaming the
FASTQ; see §9):

- Shards emit a 4-byte header-only record for flex-discarded ranges (qpos, n_positions, count,
  `n_cells=0`) — cheap, no payload.
- The rejoin evaluates the global recovery condition. Reads that trigger it get their sequence
  copied into a side buffer and are **deferred** rather than finished.
- After the rejoin, one fixup pass loads each shard once and resolves only the deferred reads at
  flex=128. At 6 GB/shard and ~2 GB/s that is ~30 s of shard loading total, over a read set that
  is presumably ~1 % of the sample.

The cost is a whole extra phase and a second code path, so do not build it on spec — build it if
and only if the §4 measurement says the recovery path is load-bearing.

## 5. What is decided where

| Decision | Code | Scope | Where |
|---|---|---|---|
| `positions.len()` (sort key) | `StdRangeExtractor` | key-local | shard (emit it) |
| header present? | `FMValues::get_range` (`size > HEADER_THRESHOLD`) | key-local | shard |
| `min_dist` over flank headers | `retrieve_seeds` | key-local | shard (emit it) |
| flex filter `count <= max_best_flex` | `retrieve_seeds` | key-local once the threshold is fixed | **shard** (§4) |
| range ordering by `positions.len()` | `StdRangeExtractor` | **global** | rejoin |
| `break` at `--ranges` | `retrieve_seeds` | **global** | rejoin |
| `--min-ranges` recovery | `generate` | **global** | **dropped** (§4) |
| `glidesort` by `(rval, rpos)` | `generate` | **global** | rejoin |

With recovery gone, the only surviving cross-key dependency is the rank-order cutoff — and it is
free: the rejoin sorts and walks ~20 records per mate.

Note the cutoff **barely binds** at N=10. A 150 bp read yields ~23 syncmers → `r ≈ 20` ranges, of
which the headered ones are capped at 15; each shard sees ~2. So a shard emitting everything it
has (post-flex-filter) is close to what the global rank would keep anyway. The rejoin still
applies the cutoff properly — it costs nothing — but do not expect shard-side pruning to pay.

### Distributed top-k (a free cap, not a win at N=10)

A range in the global 15-rarest must also be in its own shard's 15-rarest, so each shard may cap
itself at 15 *headered* ranges per read (headerless ranges are always emitted; they are the rarest
and cost ≤2 cells). At ~2 ranges/shard/read it never binds — but it bounds the pathological case
for free.

## 6. Read identity — the batch is the unit of everything

Evidence needs a read id, and the rejoin must re-pair it with the read. In
[bioreader](https://github.com/4less/bioreader) `parallel/fastq.rs`, worker threads pull batches
from a `Mutex<FastqPairedByteReader>` via `fill_buf`. **Batch boundaries depend only on
`buffer_size` (2^24) and the byte stream — not on thread count, not on scheduling** (also true
under `GzDecoder`; batching happens on the decompressed stream). Therefore:

- **`read_id = (batch_id: u32, idx_in_batch: u32)` packed into a u64 is stable across passes** and
  across thread counts, and sorts in file order.
- A thread runs a batch to completion before taking the next, so **the records it writes for that
  batch are contiguous and in `read_id` order**.

This is the load-bearing observation of the design: it removes the external sort, removes the
global merge, and makes the rejoin parallel at the granularity the existing driver already uses.

Note what is *not* required: any global ordering of the evidence. Only two things matter —
batch decomposition is deterministic, and each batch's evidence is contiguous. §8 exploits that
to let threads write wherever they like.

**Verified** (2026-07): `FastqPairedByteReader::fill_buf` (`fastq_byte_reader.rs:225`) calls
`read()` to fill both internal buffers to capacity — looping, with an explicit comment about
`GzDecoder` returning short reads — then `load_lines()` counts newlines in each, rounds each down
to a multiple of 4 (`& !0b11`), takes `min(lines1, lines2)`, and drains exactly that many records.
The reader state is a pure function of the byte stream and the capacity, and it sits behind a
`Mutex`, so *which* thread calls it is irrelevant to the batch sequence. This is read off the
implementation, not proven — hence the test in §13.5.

**Required upstream change (bioreader, ~30 lines):** a batch counter on `FastqPairedByteReader`,
returned from `fill_buf`, surfaced through `PairedFastqReader`, plus a
`read_fastq_paired_end_state_par` variant whose closure receives `(rec1, rec2, read_id, state)`.

**Determinism obligations.** These are correctness, not tuning:
- `buffer_size` must be identical in every shard pass and the rejoin → put it in the manifest; do
  not leave it as the `2usize.pow(24)` literal currently inlined at
  [process_fastq.rs](src/align/process_fastq.rs) call sites.
- Input bytes must be identical (record size + a cheap hash in the manifest).
- Prove it: same FASTQ, `--threads 1` vs `--threads 32`, assert an identical
  `(batch_id, idx) → record` map. Cheap test, and everything downstream rests on it.

## 7. Intermediate format

Records are grouped per `(read, mate)`, ranges in qpos order, and written a batch at a time (§8).
LEB128 varints unless stated.

```
group:
  read_id_delta   varint      // delta vs previous group in this stream; ~1 for most groups
  mate            u8          // 0 | 1
  n_ranges        varint
  range[n_ranges]

range:
  qpos_delta      varint      // delta vs previous range in this group (qpos ascending, unique)
  n_positions     varint      // == vrange.positions.len() -- the global sort key
  flags_dist      u8          // bit0: headered; bits4..7: min_dist (0..=16)
  n_cells         varint
  cells           n_cells × u64   // raw VCell, verbatim
```

The payload is **the raw `VCell`, verbatim**. `Seed::from_flexmer::<K,C,F>(qpos, rpos, value, dist)`
is a pure function of `VD::get(cell.0)` plus `qpos` and `dist`, and `qpos`/`dist` are per-*range*,
not per-seed. The shard never constructs a `Seed` — it copies 8-byte cells and the rejoin builds
the `Seed`s. All per-seed metadata is hoisted into the range header.

- Headerless ranges: `bit0 = 0`, `n_cells == n_positions ≤ 2`; the rejoin uses
  `Seed::from_coremer` (length `C`) rather than `from_flexmer` (length `K` when `dist == 0`).
- Flex-discarded ranges are **not emitted at all** (§4). Under §4a they would come back as
  `n_cells = 0` header-only records.
- A read with no hits in this shard emits nothing.
- `n_positions` stays even though the cutoff barely binds (§5) — it is 1–2 B and it is the only
  thing making the cutoff meaningful.
- Optional: delta-varint the `cells` within a range. They are likely ascending in `(rval, rpos)`
  if the build scans references in order — **verify before relying on it**. ~30 % on the dominant
  field.
- Optional: zstd-1 frame **per batch**, not per file, so the §8 index still points at a frame
  boundary. Write-once/read-once, so level 1 is the right end of the curve.

Put the encoder behind a sink so the on-disk and in-RAM paths (§10) share one implementation:

```rust
trait EvidenceSink { fn put(&mut self, buf: &[u8]); }   // impl for File and for Vec<u8>
```

## 8. Physical layout — one file per shard

```
run.manifest.json          // buffer_size, N, shard boundaries, input hashes, K/C/F/S, params
evidence.s{shard}.bin      // N files
evidence.s{shard}.idx      // flat array: batch_id -> (offset u64, len u32)
```

A thread buffers a whole batch's records in RAM (~3.5 MB per shard per batch), then reserves
space atomically and writes at that offset — **no lock held during IO**:

```rust
let off = shard_tail.fetch_add(buf.len() as u64, Ordering::SeqCst);
file.write_at(&buf, off)?;                  // pwrite; threads write in parallel
idx.push((batch_id, off, buf.len()));       // merged into the flat array at pass end
```

Batches therefore land in the file in **arbitrary order** — whichever thread reserves first — and
that is fine: `idx[batch_id]` names the blob explicitly. Within a blob, records are in `read_id`
order because one thread wrote it walking the batch in order. That is the only ordering the
rejoin needs (§6).

The index is a flat array indexed by dense `batch_id`: ~1800 entries (100 M reads at ~55 k
reads/batch for a 2^24 buffer) × 12 B ≈ **21 KB per shard**. `len = 0` for a batch with no hits.

**No consolidation pass, no sort, and no per-thread files.** Write buffering costs ~3.5 MB × T per
pass (~110 MB at T=32, and only one shard pass is resident at a time).

Rejected alternative: one file per `(shard, thread)`, N × T = 320 files, index carrying a
`thread_id`. It works and needs no write buffering, but it is more files and more machinery for
nothing — the index does the joining either way.

**Do not reuse `OutputBuffer`** ([src/io/output_buffer.rs](src/io/output_buffer.rs)) for this. It
flushes thread-local buffers into a shared `Arc<Mutex<OutputTarget>>`, which interleaves chunks
from different threads into one stream and destroys per-batch contiguity — the one property
everything here rests on.

Crash handling: append the index entry only after the batch's bytes are flushed (and `fsync` the
data before the index at pass end), so a truncated pass yields a short-but-valid index. A pass is
complete iff its index covers every batch id; otherwise re-run that shard. Shards are independent.

## 9. The rejoin

Drive it with the **existing parallel FASTQ reader**, unchanged. A worker that picks up batch `B`:

1. Looks up `idx[s][B] → (offset, len)` for each of the N shards.
2. Reads those N blobs with **`FileExt::read_at`** and opens a cursor over each. Each blob is
   sorted by `read_id`.
3. Walks the batch's reads in order. For read `r`, advance each cursor while
   `cursor.read_id == r`, collecting range records per mate.

```
B=417 -> idx[0] -> (off=0x1A40, len=3.4M)  \
B=417 -> idx[1] -> (off=0x9C20, len=3.1M)   |  N contiguous blobs,
...                                          |  each sorted by read_id
B=417 -> idx[9] -> (off=0x4410, len=3.6M)  /
```

**No heap, no merge tree.** Reads are consumed in increasing `read_id` and every cursor is sorted
on the same key, so N monotone cursors suffice: a linear scan, parallel at batch granularity, with
no coordination between workers.

**Use positioned reads.** Workers read different offsets of the same N files concurrently;
`seek` + `read` is stateful per file descriptor and shared, so two workers will silently corrupt
each other's reads. `read_at` (pread) has no such state and no fd budget. Per-worker fds
(N × T = 320) also work but buy nothing.

Then, per mate:

```
ranges.sort_unstable_by_key(|r| (r.n_positions, r.qpos))     // §3.3: explicit tie-break
for range in ranges:
    if range.headered:
        for cell in range.cells:
            let (val, pos) = VD::get(cell);
            seeds.push(Seed::from_flexmer::<K,C,F>(range.qpos, pos, val, range.min_dist))
        retrieved += 1
    else:
        for cell in range.cells:
            let (val, pos) = VD::get(cell);
            seeds.push(Seed::from_coremer::<K,C,F>(range.qpos, pos, val))
    if retrieved >= max_ranges { break }                      // headered only -- §3.2
glidesort::sort_by_key(&mut seeds, |s| (s.rval, s.rpos))
```

That is `StdSeedExtractor` with `VRange` swapped for `RangeRecord` and the flex test already
applied upstream. Implement it as a second `SeedExtractor`-shaped type (`ReplaySeedExtractor`) and
**factor the shared walk out of `StdSeedExtractor`** so the two cannot drift. The walk must exist
once in the codebase, or the sharded and unsharded modes diverge the first time someone tunes the
heuristic — and with exactness relaxed, nothing will catch it.

The worker then continues into today's `ModularPE` from `AnchorExtractor` onward, unchanged.
Because the rejoin re-reads the FASTQ in the same batch order, the reads needed for extension,
alignment and output are **already in hand** — no read cache, no random access, no seeking into
gzip.

**Merge memory:** ~55 k reads/batch × ~640 B ≈ 35 MB of evidence per in-flight batch; ×32 workers
≈ 1.1 GB, less if the cursors stream rather than slurp.

## 10. The degenerate case: skip the disk entirely

The intermediate only needs to exist when a sample's evidence exceeds RAM. At ~640 B/pair, 32 GB
holds ~50 M pairs — often a whole sample. The same code path with an in-RAM `EvidenceSink` (§7)
gives shards-outer with **no intermediate files at all**: the per-batch blobs stay in a
`Vec<Vec<u8>>` indexed by `batch_id` per shard — the §8 index degenerates into the outer `Vec`'s
own indexing — and the §9 rejoin walks memory ranges instead of file ranges.

Make the sink and the cursor agnostic to this from day one. It is the difference between a usable
single-sample mode and a design that only pays off on a cluster — and the in-RAM path is far
easier to debug, so build it first.

For many samples against a fixed DB, shards-outer amortises the shard load across all samples, but
requires holding **every** sample's evidence until the last shard finishes: `S × ~64 GB`. At 100
samples that is ~6.4 TB of scratch. The amortisation and the intermediate size pull in opposite
directions and there is no clever way out under key-sharding; it is a disk budget decision.

## 10a. Efficiency: the passes are gzip-bound, not IO-bound

The intuition that this design does "a lot of jumping around in files" does not survive contact
with the numbers — but a different cost does, and it is ~10× bigger.

### The rejoin's access pattern is fine

- **The FASTQ is read contiguously.** The rejoin streams it with the existing driver in the same
  batch order; there is no random access into reads at all. That is precisely why §9 re-reads the
  FASTQ instead of caching reads — read `r`'s sequence arrives exactly when its evidence does.
- **The evidence is large sequential reads.** Per batch (~55 k pairs): N=10 `read_at` calls of
  ~3.5 MB. That is **one seek per 3.5 MB — one seek per ~5,500 reads.**

| per batch | |
|---|---|
| evidence IO | 10 × 3.5 MB ≈ **13 ms** on NVMe (~0.1 ms seek amortized over ~1.2 ms transfer) |
| alignment compute | 55 k pairs × (chain + extend + WFA2) → **hundreds of ms**, optimistically |

Across a sample: ~70 GB read in 3.5 MB chunks ≈ **35 s at 2 GB/s**, against a run measured in tens
of minutes — 1–2 %. 32 workers × N files also gives the device a deep queue, which NVMe wants.

**So the rejoin+align phase ≈ today's unsharded aligner, minus the map lookups, plus ~35 s of
streaming reads.**

Two caveats:
- **NVMe scratch is a requirement.** On spinning disk or a network filesystem, 10 seeks/batch ×
  32 workers thrashes and this analysis inverts.
- **`buffer_size` is the IO knob**, not just the determinism constant (§6): it sets blob size
  directly. 2^24 → 3.5 MB blobs; 2^26 → 14 MB blobs, fewer and larger reads.

### The real cost: N+1 × single-threaded inflate

**`GzDecoder::read` runs inside the byte-reader mutex.** `PairedFastqReader::load_batch_par` locks
the shared `FastqPairedByteReader`; `fill_buf → read()` calls `self.file1.read(..)` on the
`GzDecoder` while holding that lock (`fastq_byte_reader.rs:191-206`). Inflate is therefore
single-threaded *and* serialized against every worker.

flate2 inflates at ~150–250 MB/s on one core. A 100 M-pair sample is ~60 GB uncompressed →
**~5 min of serialized inflate per pass**. Shard passes do little compute by comparison (1/N of
the lookups over T threads), so they are **~100 % gzip-bound**, paid N+1 times ≈ **~55 min of pure
decompression** at N=10.

> **Independent finding:** this already throttles the *unsharded* aligner. At 32 threads, 200 MB/s
> of inflate feeds workers that starve unless per-read work exceeds ~50 µs. Worth measuring on main
> regardless of ProjectShard. Note `gzp` (parallel gzip) is in `Cargo.toml` but used nowhere;
> `flate2::read::GzDecoder` is what is wired up at
> [process_fastq.rs](src/align/process_fastq.rs):98,126,307,334.

Fixes, cheapest first:

1. **Decompress once to scratch.** The critical section becomes a 16 MB memcpy from page cache
   (~3 ms) instead of ~80 ms of inflate. No format change; ~60 GB scratch. Do this first.
2. **`gzp`** — parallel inflate needs block-structured input (BGZF/mgzip); plain gzip cannot be
   split without an index. Only helps if you control the input format.
3. **Partition the syncmers in a pass 0.** Parse the FASTQ once, extract syncmers once, and write
   each to the shard file its core-mer belongs to (same batch-blob scheme as §8). Then pass `s`
   streams only its own partition (~1/N of the data) and never touches the FASTQ; the FASTQ is
   inflated exactly twice (pass 0 + rejoin). Kills the N× inflate **and** the N× syncmer
   extraction (hurdle 12).

   Record: `read_id u64 | qpos u16 | Kmer<31> u64` ≈ 18 B raw, ~10–12 B delta-coded within a read.
   At ~23 syncmers/mate: **~46–83 GB**, i.e. comparable to the evidence itself, and each pass reads
   only 1/N of it.

   This is the structurally right answer — the same partition-then-regroup trick applied one stage
   earlier — but it is a second record format and roughly doubles scratch. Build it only if
   measurement says syncmer extraction and inflate actually dominate.

## 11. Sizing

Per mate, with `r` = ranges with a hit surviving the flex filter, `c` = mean cells per range:

```
bytes/mate ≈ r × (4 + c × 8)      (+ ~3 B per (read,mate,shard) group header)
```

Estimates: 150 bp read, closed syncmer density `≈ 2/(C-S+2) = 2/13` → ~23 syncmers/mate; most hit
a comprehensive DB → `r ≈ 20`; `c ≈ 2`.

```
per mate     ≈ 20 × (4 + 16)                    ≈  400 B
per pair                                        ≈  800 B  (~640 B with varint deltas)
per sample   100 M pairs                        ≈ 64–80 GB
```

Total is **independent of N** — each seed is emitted once, by whichever shard owns its key. N only
duplicates group headers (~3 B × N × 2 per pair ≈ 60 B/pair at N=10).

**Measure before trusting any of this.** Instrument `StdSeedExtractor` on a real sample to
histogram `ranges.len()`, `positions.len()`, and `count_at_min_dist`. The intermediate size rides
entirely on `r` and `c`, and the `count_at_min_dist` histogram is what prices §4. This measurement
is the go/no-go for the design.

## 12. Hurdles

| # | Hurdle | Mitigation |
|---|---|---|
| 1 | `--ranges` cutoff is a global rank; shards can't apply it | Emit at range granularity with `n_positions`; rejoin applies the cutoff (§5). Barely binds at N=10 |
| 2 | Flex threshold retroactively raised 16→128 by a global condition | **Drop the recovery path**; shard applies flex=16 locally; re-tune `--max-best-flex` (§4). §4a if measurement objects |
| 3 | `sort_unstable_by_key` tie-break depends on pre-sort order | Sort by `(n_positions, qpos)` — explicit and deterministic (§3.3) |
| 4 | Recovery path duplicates seeds (existing bug) | Fix on main first — it pollutes the baseline you tune against (§3.4) |
| 5 | No read id in bioreader | Add batch counter; `read_id = (batch_id, idx)` (§6) |
| 6 | `buffer_size` inlined at call sites; a change silently corrupts the rejoin | Move to manifest; assert on load (§6) |
| 7 | Batch determinism is assumed, not proven | Test: 1 thread vs 32 → identical id→record map (§6) |
| 8 | ~64–80 GB scratch/sample; `S ×` that for shards-outer over many samples | In-RAM sink when it fits (§10); otherwise a disk budget decision |
| 9 | `stats` (`ranges`, `retrieved_ranges`, `seeds`) accumulate in the extractors → counted N× | Rejoin owns the logical stats; shard passes report only IO/pass stats |
| 10 | Partial/crashed shard pass | Index appended after flush; complete iff it covers all batch ids; shards re-runnable independently (§8) |
| 11 | Two copies of the walk will drift — and with exactness relaxed, nothing catches it | Factor the walk out; `Std`/`Replay` share it (§9). **The main correctness risk in the design** |
| 12 | **Inflate is serialized inside the byte-reader mutex, and paid N+1 times — the dominant cost (§10a)** | Decompress once to scratch (cheap); syncmer partition pass 0 if extraction also dominates |
| 13 | Rejoin needs NVMe scratch; HDD/network FS inverts the §10a analysis | Requirement, not a mitigation. `buffer_size` raises blob size if seeks bite |

## 13. Phasing

Steps 1–3 are worth doing regardless of whether ProjectShard proceeds.

1. **Measure** (§11): histogram `r`, `positions.len()`, `count_at_min_dist`; and what fraction of
   reads trigger recovery (§4). Decides go/no-go and prices §4.
1a. **Measure whether the current aligner is already gzip-bound** (§10a): run a sample from `.gz`
   vs. uncompressed at T=32 and compare wall clock. If they match, inflate is the ceiling — that
   is a bug on main worth fixing before it gets multiplied by N+1 here.
2. **Fix the recovery-path duplication** (§3.4) on main — it is the baseline you will tune against.
3. **Sweep "no recovery, `--max-best-flex` ∈ {16,24,32}" against unsharded** using the existing
   `GOLDSTD_EVAL` harness on a simulated sample. This settles §4 **before any sharding code
   exists**, which is the cheapest possible place to settle it.
4. **Factor the `retrieve_seeds` walk** out of `StdSeedExtractor` so it can be driven by `VRange`
   or `RangeRecord` (hurdle 11).
5. **bioreader: batch id + read_id closure** (§6); prove determinism with the 1-vs-32-thread test.
6. **`flexalign shard`**: slice an existing index into N + manifest (§2). Verify
   `∀k: sharded.get_vrange(k) == unsharded.get_vrange(k)` — a cheap exhaustive test over all 2^30
   keys, and the one place a bit-exact guarantee is both achievable and worth having.
7. **`ShardedDB` + evidence emitter**: shard pass writes §7 records + §8 index.
8. **`ReplaySeedExtractor` + rejoin driver** (§9), in-RAM sink first (§10).
9. **Acceptance**: `GOLDSTD_EVAL` sensitivity/precision of sharded vs unsharded on a simulated
   sample, plus a seed-level diff at the cut line to confirm the only differences are the ones
   §4 predicts. With exactness relaxed, this comparison *is* the spec — there is no
   byte-equality backstop.
10. Disk sink, batch index, parallel rejoin, multi-sample shards-outer.

## 14. Open questions

1. Real `r` (ranges/mate) and `c` (cells/range), and the `count_at_min_dist` histogram (§11).
2. What fraction of reads trigger the `--min-ranges` recovery path? Decides §4 vs §4a outright.
3. Scratch budget per sample, and how many samples per DB in one shards-outer run (§10)?
4. Is the recovery-path duplication (§3.4) known/intentional?
5. Do cells within a value block ascend in `(rval, rpos)`? Decides the ~30 % delta-coding win (§7).
6. Is `--max-range-size 256` at build time a deliberate sensitivity decision? It bounds
   `n_positions`, hence `c`, hence §11.
