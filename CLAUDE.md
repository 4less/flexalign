# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Overview

`flexalign` is a Rust read aligner/mapper for short reads against a reference of
nucleotide sequences. It builds a flexible k-mer index (via the `flexmap` crate),
maps single- or paired-end FASTQ reads through a staged seed-and-extend pipeline,
and emits alignments (currently PAF; SAM support is in progress on the
`paf_sam_support` branch — see `src/align/sam.rs`).

The codebase depends on several sibling crates pulled from git
(`kmerrs`, `flexmap`, `bioreader`) and on `libwfa2` for gapped alignment.

## Dependency ecosystem (4less sibling crates)

These three crates are maintained alongside flexalign and provide the buildingibwa 
blocks of the pipeline. In `Cargo.toml` they are GitHub git-dependencies; locally
they are linked to working copies (see "Local dependency linking" below).

| Crate         | GitHub (org`4less`) | Branch            | Local path                           | Role in flexalign                                                                                                                                                                                  |
| ------------- | --------------------- | ----------------- | ------------------------------------ | -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `bioreader` | `4less/bioreader`   | `main`          | `../bioreader`                     | Buffered & parallel FASTA/FASTQ readers and sequence record types (`FastaReader`, `RefFastqRecord`, `read_fastq_*_state_par`, `is_gzip`). Drives read input in `align/process_fastq.rs`. |
| `kmerrs`    | `4less/kmerr`       | `main`          | `../kmerr` (crate name `kmerrs`) | k-mer / syncmer / minimizer extraction (`ClosedSyncmer`, `Kmer`, `Minimizer`). Backs the `KmerExtractor` stage.                                                                            |
| `flexmap`   | `4less/flexmap`     | `zeroable_load` | `../flexmap`                       | The flexible k-mer index itself — build, on-disk blob storage, and range lookup (`Flexmap`, `FlexmapBlob`, `VRange`, `build::default_build`). Wrapped by `database::flexmap::DB`.       |

### Local dependency linking

The sibling repos live under `/home/fritscher/git/4less/`. They are linked into
builds via a `paths` override in `.cargo/config.toml` (gitignored, machine-local),
so cargo compiles them from the local working copies instead of fetching GitHub.

Workflow: edit code in `../bioreader` / `../kmerr` / `../flexmap`, build & test
flexalign here to pick the changes up immediately, then **push that dep repo to
GitHub** to publish. The committed `Cargo.toml` keeps pointing at the GitHub
sources, so the canonical build is unaffected; delete `.cargo/config.toml` to build
purely against GitHub. (Verify the override with
`cargo metadata | grep -o '"manifest_path":"[^"]*flexmap[^"]*"'`.)

## Toolchain & Build

This crate **requires Rust nightly** — it uses unstable features (`#![feature(...)]`
in `src/lib.rs` / `src/main.rs`: `test`, `generic_const_exprs`, `unboxed_closures`,
etc.) and the `profile-rustflags` cargo feature in `Cargo.toml`.

```bash
# Build (release is the meaningful target — it sets AVX2/FMA target-features and fat LTO)
cargo +nightly build --release
cargo +nightly build              # debug

# Run
cargo +nightly run --release -- --help
./target/release/flexalign --help

# Tests / a single test
cargo +nightly test
cargo +nightly test it_works                       # by name
cargo +nightly test --package flexalign <name>

# Benchmarks (uses libtest `#[bench]`; gated behind the `benchmarking` feature)
cargo +nightly bench --features benchmarking

cargo +nightly clippy --release
```

Note: `Cargo.lock` is in `.gitignore`, and `/data`, `/target`, `/results`,
`/analysis` are ignored. The `data/` directory holds local test fixtures
(reference FASTA, prebuilt `.flex.*` index files, sample reads).

### Compile-time gold-standard evaluation

`GOLDSTD_EVAL` (in `src/lib.rs`) is a **compile-time constant** read from the
`FLEXALIGN_GOLDSTD_EVAL` env var at build time. When set to `1`, the binary is
built to compare each read's mapped reference against the truth encoded in the
read header (simulated reads) and to emit accuracy stats / false-positive reads.

```bash
FLEXALIGN_GOLDSTD_EVAL=1 cargo +nightly build --release
```

Because it is `const`, the gold-standard branches are compiled out entirely when
unset — changing it requires a rebuild, not just a flag.

## Running

The CLI (defined in `src/options.rs`, `clap` derive) takes paired or single reads
and a reference:

```bash
flexalign -r REFERENCE.fna -1 reads_R1.fq[.gz] -2 reads_R2.fq[.gz] [-o OUTDIR] -t THREADS
flexalign -r REFERENCE.fna -1 reads.fq                              # single-end
```

- On first run against a reference, the index is **built and saved** next to the
  reference as `REFERENCE.flex.index` / `.flex.id2ref` / `.flex.ref2id`; later runs
  load it. Use `--force-build` to rebuild.
- With multiple input files, `-o`/`--output` must be a directory; output prefixes
  are inferred per input.
- **stdout carries the alignment output (PAF).** Never print color/log to stdout —
  `main.rs` notes this and all diagnostics go to stderr via `log` + `simple_logger`.
  Verbosity is `--log-level 0..5` (Off..Trace).
- Key tuning args map directly onto pipeline stages: `--ranges`, `--max-range-size`,
  `--max-best-flex`, `--extend-top-x`, `--align-top-y`, `--min-ranges`.

## Architecture

Entry: `main.rs` → `flexalign::run` (`src/flexalign.rs`). `run` constructs/loads the
database, then calls `process_fastq_wrapper_modular` (`src/align/process_fastq.rs`),
which reads FASTQ in parallel (`bioreader`, `--threads`) and drives one pipeline
instance per worker.

The index `DB<K, C, F, S, L, CELLS_PER_BODY, HEADER_THRESHOLD>` is heavily
**const-generic**; these parameters are fixed at the `run` call site in
`src/flexalign.rs` (currently K=31, C=15, F=16, S=4, …) and threaded through the
whole pipeline as generics. Changing the k-mer/syncmer scheme means changing these
constants, not config.

### The modular pipeline

The core design is a **trait-per-stage pipeline** assembled generically. Traits are
declared in `src/align/common.rs`; standard implementations live under
`src/align/process/`. The pipeline structs `Modular` (single-end) and `ModularPE`
(paired-end) in `src/align/modular_workflow.rs` hold one impl per stage and run them
in order:

1. `KmerExtractor` — extract relevant k-mers/syncmers from a read
   (`process/kmer_extractor.rs`, uses `kmerrs::ClosedSyncmer`).
2. `RangeExtractor` — look up index ranges for each core-mer
   (`process/range_extractor.rs`).
3. `SeedExtractor` — turn ranges into concrete seeds, bounded by
   `max_best_flex` / `max_range_size` / `min_ranges` / `ranges`
   (`process/seed_extractor.rs`).
4. `AnchorExtractor` / `PairedAnchorExtractor` — group seeds into anchors
   (`process/anchor_extractor.rs`, with `anchor_extractor_revised.rs`,
   `anchor_sorter.rs`, `seeds_to_anchors.rs` as alternatives/helpers).
5. `PairedAnchorExtender` — extend top-`extend_top_x` anchors against the reference
   using Hamming distance (`process/anchor_extender.rs`). This begins the
   reference-comparison phase; before it, no query/reference sequence comparison
   happens.
6. `Align` + `Heuristic` — gapped alignment of the top-`align_top_y` anchors via WFA2
   (`process/alignment.rs`, `LIBWFA2Alignment`; ANI-based abort scoring).
7. `PAFOutput` / `SAMOutput` — write results (`process/output.rs`). Output is wrapped
   in an `Or<PO, SO>` (`src/align/common.rs`) so PAF and/or SAM can be selected; SAM
   is currently stubbed as `NoSAMOutput`.

MAPQ is computed as a pseudo-mapq from the score gap between best and second-best
anchor (`StdPairedAnchorMAPQ`). Anchor data structures are in
`src/align/data_structures/` (`anchor.rs`). Per-read timing/accuracy counters
accumulate into `Stats` (`src/align/stats.rs`), threaded through every stage.

To add or swap a pipeline stage: define/extend a trait in `src/align/common.rs`,
add an impl under `src/align/process/`, and wire it into the `Modular`/`ModularPE`
construction in `src/align/process_fastq.rs`.

### Database layer

`src/database/` — `common.rs` defines the `FlexalignDatabase` trait
(`get_rname` / `get_reference` / `get_vrange` / `build` / `save` / `load`) and
`DBPaths` (the `.flex.*` file naming). `flexmap.rs` implements it over
`flexmap::FlexmapBlob`, persisted with `savefile`. `build.rs` handles index
construction. `GLOBAL_VERSION` (`src/lib.rs`) is the on-disk format version checked
on save/load.

### Supporting modules

- `src/io/output_buffer.rs` — buffered, thread-shared output sink (`OutputTarget`
  = stdout or file, behind `Arc<Mutex<…>>`).
- `src/distance/` — Hamming/indel distance helpers (`indel_detection.rs`).
- `src/align/eval.rs`, `process/evaluate.rs` — gold-standard accuracy evaluation
  (only meaningful when `GOLDSTD_EVAL` is compiled in).
- `src/utils.rs`, `src/misc.rs` — output-prefix inference and misc helpers.

## Working in this codebase notes

- `process_fastq.rs` contains both an older `process_fastq_wrapper` (uses
  `workflow::Standard`) and the current `process_fastq_wrapper_modular`. `run` calls
  the modular one; prefer it for changes.
- `modular_workflow.rs` carries large blocks of commented-out debugging/exploration
  code and several `eprintln!`/`println!` left in the alignment loop — be aware some
  of these print to stdout and can corrupt PAF/SAM output.
