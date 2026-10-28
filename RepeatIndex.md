# Repeat index: seeds for the cores the main index drops

Status: implemented on branch `repeat-index` (2026-10-28) as `--repeat-index`, off by default, in
its eager form: every main-index miss asks the repeat index. The lazy trigger (section 3), the MAPQ
cap (section 4) and sharded runs (section 5) are not implemented. Measured on GRCh38 (no-alt
analysis set) with flexalign `a1ff1ad` and one simulated sample of 946,262 pairs
(`datasets::human` of the benchmark).

## Results of the implementation

Correct placements (recall) / precision, as the benchmark scores them. The read pass is the median
of 3 runs on 16 threads, SAM to `/dev/null`.

GRCh38, 1.89 M mates:

| | MAPQ ≥ 0 | MAPQ ≥ 1 | MAPQ ≥ 10 | MAPQ ≥ 30 | read pass | load | peak RSS |
| --- | --- | --- | --- | --- | --- | --- | --- |
| minibwa | 98.26 / 98.26 | 97.52 / 99.96 | 97.31 / 99.99 | 96.33 / 100.00 | | | |
| flexalign `a1ff1ad` | 97.57 / 98.41 | 96.25 / 99.78 | 95.57 / 99.88 | 93.95 / 99.92 | 0.75 s | 0.55 s | 8.4 GB |
| `--repeat-index` (bucket cap 64) | 97.77 / 98.42 | 96.70 / 99.77 | 96.06 / 99.87 | 94.58 / 99.92 | 0.87 s | 0.98 s | 9.8 GB |
| bucket cap 16 | 97.72 / 98.43 | 96.57 / 99.77 | 95.91 / 99.86 | 94.43 / 99.93 | 0.81 s | | |
| bucket cap 128 | 97.78 / 98.40 | 96.74 / 99.75 | 96.11 / 99.86 | 94.65 / 99.92 | 0.95 s | | |

- **Gain:** the repeat index recovers about 3,800 placements at MAPQ ≥ 0 and 12,000 at MAPQ ≥ 30;
  it closes about a quarter of the MAPQ ≥ 30 gap to minibwa.
- **No precision gain:** the overconfident MAPQ on dropped repeats is still there, so the MAPQ cap
  of section 4 is still needed.
- **Repeat file:** 1.84 GB on GRCh38 (81,540 cores, 107.7 M positions), built in 20 s on 22 threads.

Bacteroides (1.99 M mates):
- identical alignments within one mate;
- the read pass is the same within noise;
- the repeat section is tiny: 714 cores.

Metagenomics, 5 samples of the full protal reference (9.95 M mates):

| | MAPQ ≥ 0 | MAPQ ≥ 10 | MAPQ ≥ 30 | read pass (1 run) | peak RSS |
| --- | --- | --- | --- | --- | --- |
| off | 4.066 / 85.12 | 3.425 / 95.57 | 2.562 / 98.45 | 7.31 s | 27.2 GB |
| `--repeat-index` | 4.090 / 85.01 | 3.425 / 95.80 | 2.523 / 98.60 | 7.64 s | 40.4 GB |

On that reference the dropped cores are conserved marker genes, not repeats:
- 1.26 M cores, 688 M positions;
- the repeat file is 13.4 GB, built in 146 s with a 54 GB peak;
- the changes are small and mixed, and memory grows by 13 GB.

Keep it off for marker-gene databases.

### Lookup cost: what made it affordable

A first version binary-searched each core's sorted flanks. That doubled the GRCh38 read pass
(0.75 → 1.40 s) for +5.5 % instructions: about 20 dependent cache misses per lookup, 6.6 M
lookups. Two changes, in order:

1. **Batched lookups:** a read's misses are looked up together, in stages with prefetches, so
   their cache misses overlap.
2. **A per-core directory over the distinct flanks, slotted by a hash of the flank.** The copies
   of a young repeat share most of their flanks, so slots by flank value held thousands of
   entries; by hash, a slot holds about four distinct flanks.

Read pass: 1.40 → 1.18 (value-bit slots) → 1.08 (batched) → 0.87 s (hashed distinct-flank slots).

## Problem

The index build sets every core with more than `--max-range-size` (256) positions to zero
(flexmap `keys.rs`, `FMKeys::build`: `ckmer_count > max_range_size` → `ckmer_count = 0`). At query
time such a core looks exactly like a core absent from the genome. For reads from young repeats
(L1HS–L1PA6, alpha satellite) most cores are dropped. Their remaining rare cores point to *some*
copies, often not the true one.

Effect on GRCh38, compared with minibwa (`map -x sr`) on the same reads:

| | flexalign | minibwa |
| --- | --- | --- |
| placed correctly (MAPQ ≥ 0) | 97.57 % | 98.26 % |
| placed correctly with MAPQ ≥ 30 | 93.95 % | 96.33 % |
| precision at MAPQ ≥ 30 | 99.924 % | 100.000 % |

- 12,158 mates are placed correctly with MAPQ > 0 by minibwa and missed by flexalign.
- 3,962 mates are placed wrongly by flexalign with MAPQ ≥ 1.
- 12,278 mates are correct but get MAPQ 0 from flexalign and MAPQ > 0 from minibwa.
- 68 % of the missed mates come from young L1 and alpha satellite, which hold 4.9 % of all mates.
- For 3,213 of 3,406 confident errors in the anchor dump, the true copy was never a candidate.
- Report: https://claude.ai/artifact/5WPRGDLiBA6UssmdgViJhH

flexalign already does the right selection: among a core's positions, keep those whose stored
16-base flank is nearest to the read's (`resolve_flex`, `src/align/process/seed_extractor.rs`).
The selection never runs for repeat cores, because their positions were never stored. The cap
exists because that selection is a linear scan over all of a core's flanks. The format adds hard
limits too:
- 16-bit block offsets and `MAX_KEY_VALUESSIZE` (4,096) in `FMKeys`;
- 13 count bits (8,191 cells) per key in `HashedKeys`.

Raising the cap to 4,096 recovers about 2,000 of the 12,158 mates and costs about 35 % of
read-pass time.

## Measured: what a flank-sorted repeat section holds

Counted with flexmap's own core selection (canonical 15-mer core of each 31-mer, closed syncmer
S = 4, `mask_coremer` key). The tool is `fb/` in the session scratchpad.

Cores by number of positions in GRCh38:

| positions per core | cores | positions |
| --- | --- | --- |
| 1–255 (in the index today) | 63,642,334 | 509.3 M |
| 256–1,023 | 65,734 | 29.9 M |
| 1,024–4,095 | 12,596 | 23.5 M |
| 4,096–8,191 | 1,923 | 10.9 M |
| 8,192–65,535 | 1,571 | 28.2 M |
| ≥ 65,536 | 94 | 15.4 M |
| **dropped (> 256)** | **81,918** | **107.7 M (+21 % on the 509 M indexed today)** |

Splitting the dropped cores' positions by their full flank, i.e. by the whole 31-mer:

| positions in buckets of size | share of the 107.7 M positions |
| --- | --- |
| 1 (a unique 31-mer) | 34 % |
| ≤ 63 | 60 % |
| ≤ 255 | 69 % |
| ≤ 4,095 | 91 % |
| ≥ 4,096 | 9 % |

Probing the reads, the question is whether the read's own exact (core, flank) bucket contains its
true locus, for at least one of its repeat cores:

| reads | have a dropped core | true locus in a bucket of ≤ 64 | ≤ 256 | any size |
| --- | --- | --- | --- | --- |
| missed (12,158) | 90 % | **56 %** | 68 % | 90 % |
| correct, MAPQ 0 (12,278) | 85 % | 66 % | 73 % | 85 % |
| wrong MAPQ ≥ 1, minibwa no better (1,638) | 75 % | 39 % | 46 % | 74 % |
| control, correct with MAPQ > 0 (20,000) | 47 % | – | – | – |

A half-flank fallback adds about one point: for 1 mismatch in the flank, look up the left 8
bases and the right 8 bases separately. Exact buckets carry nearly all of it. ART reads have few
errors, though; real reads with variants will lean on the fallback more.

## Design

### 1. Build: a repeat section next to the index, built by flexalign

Nothing changes in flexmap or in the main index. After the flexmap build, flexalign makes two
more passes over the reference with the same core selection and key:

1. **Count every key.** `u32` per key, 4^15 × 4 B = 4.3 GB transient; the build already peaks
   at 16.9 GB on GRCh38.
2. **Collect** `(key, flank, ref, pos)` for keys with count > `--max-range-size`, the same cap
   as the main index, so the two sets are exactly complementary. Sort by `(key, flank)`.

Single-threaded this took 3 min 19 s for GRCh38. It parallelizes per chromosome like the main
build.

Layout, as new optional sections of the bundle (`src/database/bundle.rs`):

| section | content | GRCh38 |
| --- | --- | --- |
| repeat keys | open-addressing hash, key → (first entry `u64`, entries `u32`); 2^18 slots | ~4 MB |
| repeat flanks | `u32` flank per entry, sorted per key | 431 MB |
| repeat values | packed `ref_bits + pos_bits` per entry (as `PackedValues`), same order | ~485 MB |

- In total about 0.92 GB, or +18 % on the 5.09 GB index.
- A bundle without the sections loads as today and behaves as today. The format needs a section
  flag, not a version break.
- Optional for keys with more than 4,096 entries (3,588 keys): a sampled pivot table, every 256th
  flank, so the equal-range search touches 2–3 lines instead of about 17. Measure before adding.

### 2. Query: know which cores are repeats

In range extraction, when the main lookup misses, probe the repeat hash. On a hit the core is a
**repeat core**: record its key and the read's flank for it, and count it.

- Today's read pass doesn't change; only main-table misses pay one extra probe.
- Misses are about 16 % of a read's cores, and the 4 MB repeat hash stays in cache.
- Every read now knows how many of its cores were repeats.

### 3. Query: a lazy repeat path

After the normal result for a pair (fast-hit or full path), a mate takes the repeat path when it
has at least one repeat core **and** any of:
- it is unmapped;
- its alignment has NM ≥ 1;
- its MAPQ is below 60.

Estimated on the human run, using every 15-mer as a stand-in for the cores:

| trigger | mates that trigger | missed covered | MAPQ-0 covered | wrong MAPQ ≥ 1 covered |
| --- | --- | --- | --- | --- |
| any repeat core (eager) | 63 % | 93 % | 89 % | 86 % |
| any repeat core & (unmapped, NM ≥ 1 or MAPQ < 60) | 25 % | 93 % | 89 % | 86 % |
| ≥ 25 % repeat cores & (same) | 13 % | 85 % | 77 % | 66 % |

v1 uses the second row. The share threshold stays tunable.

The repeat path, per triggered mate:
1. For each repeat core, find the equal range of the read's flank in the core's sorted flanks.
2. If the bucket has 1 to `--repeat-bucket-max` entries (default 64): emit each entry as a seed,
   as a headered range's best-flank member with distance 0 (`Seed::from_flexmer`).
3. If the bucket is larger: count the core as **unresolved**. If it is empty: v1 skips it; the
   half-flank fallback could come later.
4. Merge the repeat seeds into the mate's normal seeds and run anchoring, alignment and pair
   ranking again, for this pair only.
5. MAPQ now compares the competing copies the repeat seeds found. This is the main precision gain.

For a pair where both mates trigger, keep only repeat seeds within insert size of a candidate of
the partner, before anchoring. Two rare positions 350 bp apart are far more specific than either
alone, and this bounds the work for satellite reads. This can be a v1.1 step.

### 4. MAPQ when copies stay unseen

A mate with unresolved repeat cores (buckets > `--repeat-bucket-max`) still has copies flexalign
didn't look at. Offline rules on today's output show that capping MAPQ from the repeat share alone
costs too much recall. Applying a cap only after the repeat path ran is better calibrated. v1:
cap MAPQ at 3 when unresolved cores make up at least half of the read's cores and NM ≥ 1.
Calibrate on the human set against minibwa's per-band precision.

### 5. Shards

- The repeat section is small, so every pass can load it whole, or it can be split at the same
  key bounds as the main index (its keys are sorted).
- The lazy path needs the pair's final result. In a sharded run that exists only after the
  rejoin, so the align half loads the repeat section (0.92 GB on GRCh38).

## Expected effect, and its limits

- **Upper bound for missed mates:** 6,812 of 12,158 get a seed at the true locus with bucket
  cap 64, and 8,311 with cap 256. Each one still needs to anchor, align and win the ranking.
- **If all of them are placed:** placed-correctly (MAPQ ≥ 0) goes from 97.57 % to about 97.93 %
  (cap 64) or 98.0 % (cap 256). minibwa is at 98.26 %.
- **Precision:** the competing copies become candidates, so overconfident MAPQ should drop for
  the 66–86 % of affected mates the trigger covers.
- **Not addressed:** 717 confident errors in recent segmental duplications. Those regions have
  2–10 copies and no dropped core, so this is a different cause. Also not addressed: buckets
  above the cap (9 % of repeat positions), which mostly fall to the MAPQ cap.
- **Metagenomics:** the marker-gene reference has conserved markers and drops cores too. Check
  how many before claiming a gain or a cost there.

## Cost

| | GRCh38 |
| --- | --- |
| index size | +0.92 GB (+18 %) |
| build time | +2–3 min (two reference passes, parallel) |
| build memory | +4.3 GB transient (key counts) |
| read pass, every read | one cache-resident probe per main-table miss |
| read pass, triggered mates (~25 %) | equal-range searches on ~4–8 repeat cores, ≤ 64 extra seeds each, a second anchoring and alignment |

## Plan

| step | work | estimate |
| --- | --- | --- |
| 1 | Build the repeat section (two passes, sort, bundle sections, load) behind `--repeat-index` | 1 day |
| 2 | Repeat-core detection in range extraction; counts in stats; no behaviour change | 0.5 day |
| 3 | Lazy repeat path: trigger, seeds, re-run for the pair | 1–2 days |
| 4 | MAPQ cap for unresolved cores; calibrate on human | 0.5 day |
| 5 | Validate: human (gap, MAPQ bands, PR curve vs minibwa), bacteroides and metagenomics unchanged or better, read-pass time | 0.5 day |
| 6 | Sharded mode: load the section in the align half | 0.5 day |

The gates are those of the speed work:
- no regression on bacteroides and metagenomics alignments;
- read-pass time within noise on bacteroides;
- on human, a measured trade: time against correct placements.

## Open questions

1. Bucket cap: 64 or 256. 256 adds about 1,500 missed mates, and the seed count per triggered
   mate can then reach 4× the normal path. Decide on the human run's time and precision.
2. Trigger: is NM ≥ 1 too loose for real reads with SNVs? A threshold relative to the read's
   expected error may be needed.
3. Half-flank fallback: worth it on real reads? Measure on SRR8797712 or a GIAB sample before
   adding it.
4. Should the fast-hit path skip mates with a high repeat-core share outright? Its vote margin is
   blind to dropped copies in the same way.
