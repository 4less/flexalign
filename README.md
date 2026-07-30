# flexalign

A short-read aligner for large nucleotide references, built around a flexible k-mer index. It maps
single- or paired-end FASTQ against a reference through a staged seed-and-extend pipeline and emits
PAF or SAM.

Two properties shape how it is used:

- **The index is memory-mapped.** Opening it is near-instant and only the pages actually queried are
  read, so a multi-GB index does not have to be loaded into heap before the first read is aligned.
- **The index can be sliced into shards.** A sharded run streams over the reads once per shard,
  holding one shard resident at a time, so peak memory is bounded by the shard rather than by the
  whole index — at the cost of one pass over the reads per shard.

Alignments go to **stdout**; all diagnostics go to stderr. Never let anything else write to stdout,
or the output will be corrupt.

## Install

flexalign needs the **Rust nightly** toolchain (it uses unstable language features and the
`profile-rustflags` cargo feature) and a host C toolchain — `cc`, `make`, `cmake` — to build
libwfa2's vendored WFA2-lib.

```sh
git clone https://github.com/4less/flexalign.git
cd flexalign
just check-deps        # verify the pinned dependencies resolve, no compile
just build             # -> target/release/flexalign
```

On a cluster, where you want a binary on `$PATH` and a build that ignores any local development
overrides:

```sh
just install-cluster              # -> ~/bin
just install-cluster ~/opt/fa     # or a prefix of your choosing
```

That installs nightly via rustup, checks the C toolchain, and builds and installs `flexalign`. One
binary covers everything, sharded runs included. Load your cluster's compiler modules first if it
provides them that way.

Dependencies (`kmerrs`, `flexmap`, `bioreader`) are pinned to exact git revisions, so a build
anywhere reproduces the same binary. `flexmap` owns the on-disk index format, so leaving it floating
risks a newer library reading an index built by an older one.

## Build a database

There is no separate index command: flexalign builds the index on its first run against a reference
and reuses it afterwards. To build deliberately, run it once with a throwaway read pair.

```sh
printf '@r/1\nACGT\n+\nIIII\n' > /tmp/tiny_1.fq
sed 's|/1|/2|' /tmp/tiny_1.fq   > /tmp/tiny_2.fq

flexalign -r reference.fna -1 /tmp/tiny_1.fq -2 /tmp/tiny_2.fq -t 16 > /dev/null
```

This writes, next to the reference:

| file | what it is |
|---|---|
| `reference.fna.flex.index.blob` | the memory-mapped index (keys + values) — what alignment reads |
| `reference.fna.flex.index.keys` / `.values` | build-time artifacts, also the input to slicing |
| `reference.fna.flex.index.refseq` | reference sequences, concatenated + an offset table |
| `reference.fna.flex.id2ref` / `.ref2id` | reference name ↔ id maps |

Add `--force-build` to rebuild an existing index (needed after a change to the index format or to
the k-mer scheme). Building is the expensive step: on a 16.6 GB marker reference it takes tens of
minutes and ~35–40 GB of RAM, and the resulting files total roughly 75 GB.

The scheme constants (K, C, F, S, …) are compile-time generics fixed at the call site in
`src/flexalign.rs`, not runtime options — changing them means editing those constants and rebuilding
both the binary and the index.

## Build a sharded database

Slicing cuts an existing index into N shards. It loads the full index once, so treat it as a
one-time offline cost like the build itself.

```sh
flexalign -r reference.fna -1 /tmp/tiny_1.fq -2 /tmp/tiny_2.fq --shard-slice 4
```

This writes `reference.fna.flex.index.s4.shard{0..3}.blob` plus `…s4.shards.json`, the manifest that
records each shard's key range.

The shard count is part of every filename, so several slicings of the same reference coexist:

```
reference.fna.flex.index.blob            30.0 GB   unsharded
reference.fna.flex.index.s2.shard0.blob  15.0 GB   2-shard slicing
reference.fna.flex.index.s4.shard0.blob   7.5 GB   4-shard slicing, side by side
```

More shards means lower peak memory and more passes over the reads — each shard costs one pass, so
the wall time grows roughly with the shard count while peak memory falls.

## Align against a prebuilt database

```sh
# paired-end, SAM to stdout (the default: carries CIGAR and NM)
flexalign -r reference.fna -1 reads_R1.fq.gz -2 reads_R2.fq.gz -t 16 > out.sam

# PAF instead (mapping positions only)
flexalign -r reference.fna -1 reads_R1.fq.gz -2 reads_R2.fq.gz -t 16 --paf > out.paf

# single-end — always PAF, see below
flexalign -r reference.fna -1 reads.fq.gz -t 16 > out.paf
```

**SAM is the default**, since the aligner has already computed the CIGAR and NM that PAF would throw
away; `--paf` is the lossy option and therefore the one you ask for. Either format emits only records
that reached gapped alignment above `--min-ani` — reads below it are omitted rather than written as
unmapped.

Single-end input is the exception: that pipeline stops at anchors and never runs gapped alignment, so
it has no CIGAR to write and PAF is the only possible output. flexalign resolves that itself rather
than failing, and says so at `--log-level 3`.

**Several samples in one invocation.** `-1` and `-2` accept multiple files, and `--output` then
names a directory; the output prefix is inferred per input. This matters on a large reference,
because the index is loaded once for the whole batch instead of once per sample:

```sh
flexalign -r reference.fna \
  -1 s1_R1.fq.gz s2_R1.fq.gz s3_R1.fq.gz \
  -2 s1_R2.fq.gz s2_R2.fq.gz s3_R2.fq.gz \
  -t 16 --output results/
```

## Align with a sharded database

Sharded alignment is the **same binary** — there is no separate executable. flexalign discovers the
slicings that exist beside the reference and drives the shard passes itself:

```sh
flexalign -r reference.fna -1 reads_R1.fq.gz -2 reads_R2.fq.gz -t 16 --shards 4 > out.sam
```

If exactly one slicing exists it is detected and used without `--shards`. If several do, flexalign
lists them and exits rather than choosing: the shard count decides the run's memory profile, so
picking one silently would make that depend on whatever happens to be on disk. `--no-shards` ignores
any slicing and uses the whole index; `--shards N` creates the slicing if it does not exist yet.

For several samples, pass them all in one invocation — that puts the samples *inside* the shard loop,
so each shard is memory-mapped once for the whole batch instead of once per sample, and the full index
of the rejoin half is loaded once rather than N times:

```sh
flexalign -r reference.fna \
  -1 s1_R1.fq.gz s2_R1.fq.gz s3_R1.fq.gz \
  -2 s1_R2.fq.gz s2_R2.fq.gz s3_R2.fq.gz \
  -t 16 --shards 4 --output results/
```

Sharded runs are paired-end only, and `-t` applies to them exactly as to an unsharded run.

## Bounding memory

**By default the whole database is read into memory during the load** — index *and* reference, each
in one sequential pass. `mmap` on its own transfers nothing: the multi-GB read still happens, one 4K
fault at a time, in whatever order queries touch it, and lands inside read processing where it looks
like compute rather than like the load it is. Reading it up front moves the same bytes, lets the
kernel read ahead, and makes the reported load figure real.

That leaves the reference resident, which on a machine with free RAM costs nothing — those pages are
reclaimable. Under a memory limit it is what gets a run killed, so `--ref-budget` caps the resident
reference and drops pages once it is exceeded, re-faulting what it needs again:

```sh
flexalign … --ref-budget auto     # read the cgroup limit (container, slurm, systemd), else MemAvailable
flexalign … --ref-budget 8G       # or an explicit size
```

`auto` takes 60% of the limit, leaving room for reads, evidence and output buffers. Under an 8 GB
cgroup cap a 4-shard run completes with a budget where it is OOM-killed without one. The cost is
re-faulting: roughly 1.7–2× on the rejoin phase at a 6 GB budget.

`--lazy-ref` goes further: it skips the up-front reference read entirely and implies
`--ref-budget auto`. Use it when RAM is the binding constraint rather than time — alignment touches
only the few hundred bases around each anchor, so most of a large reference is never needed. The
index is always read up front regardless: k-mer lookups hit its whole extent in essentially random
order, so deferring that read does not avoid it, it only converts one streaming read into millions of
faults.

## Tuning

| flag | default | effect |
|---|---|---|
| `--paf` | off (SAM) | emit PAF instead of SAM |
| `--lazy-ref` | off (eager) | skip the up-front reference read; implies `--ref-budget auto` |
| `--min-ani` | 0.85 | identity floor; also sets the aligner's abort ceiling |
| `--min-query-coverage` | 0.70 | fraction of the *alignable* window that must align |
| `--ranges` | 15 | how many minimizer ranges are looked up |
| `--max-range-size` | 256 | occurrence cap per minimizer |
| `--max-best-flex` | 16 | values kept per range after flank ranking |
| `--extend-top-z` | 32 | anchor pairs extended (the first stage that reads reference) |
| `--align-top-y` | 4 | anchor pairs sent to gapped alignment |
| `--min-ranges` | 4 | floor on ranges, for the recovery pass |

Identity and coverage are measured against the bases that *can* align given where the anchor sits —
not the whole read. A read straddling the end of a reference sequence has nowhere to put its
overhanging bases, and charging it for them would reject a high-identity alignment.

Diagnostics, all off by default: `--dump-anchors` (per-read candidate list as TSV on stderr — use
`-t 1` so lines do not interleave), `--time-log FILE` (per-stage timings as TSV), `--log-level 0..5`.

Two experimental flags are present but measured as regressions on a marker-gene reference, and are
off for that reason: `--end-to-end` (rejects clipping that no reference edge explains) and
`--mapq-from-alignment` (MAPQ from aligned rather than anchor scores).

## Development

```sh
just build          # release
just test           # unit tests + the regression check
just run -- -r ref.fna -1 r1.fq -2 r2.fq
just lint
```

`.cargo/config.toml` (gitignored) can redirect the git dependencies to sibling working copies for
local development. It is deliberately absent from cluster builds — `just install-cluster` moves it
aside — because a path override silently defeats the pinned revisions.

The benchmark that measures this aligner against strobealign, bowtie2, minibwa and protal lives in
its own repository, [flexalign_benchmark](https://github.com/4less/flexalign_benchmark).
