# Uniqueness flags, inactive: no cost

Question: does flexmap's uniqueness-flag support (flexmap `uniqueness-flags`, 9c1c618) slow
flexalign down when no index uses it? Answer: no. Same index files, same alignments, same speed.

## Setup (2026-10-28, 12-core desktop, desktop apps running)

- `flexalign_base`: flexalign main (7ea91cf) against flexmap main (2fde119).
- `flexalign_flags`: the same flexalign against flexmap 9c1c618 (this branch only bumps the pin).
- Same toolchain (`cargo +nightly build --release`), same target dir.
- Reference: the first genome of the 60 largest Bacteroides species of GTDB r232 (258 MB, 10,294
  contigs). Reads: 1 M simulated 150 bp pairs, 1% substitutions (`make_data.py`).
- Each binary builds its own index; runs are interleaved (base/flags, then flags/base) and
  warm-cache.

## Results

| check | result |
| --- | --- |
| `.flex.index.compact` (620 MB), `.hashed` (660 MB), `.refseq`, `.id2ref` | byte-identical |
| `.flex.ref2id` | same entries, different order (a serialised HashMap; differs run to run) |
| PAF, 1 M pairs, 8 threads | identical after sorting (threads emit in any order) |
| wall, 1 M pairs, 8 threads, 9 runs each (`times_t8.tsv`) | 5.96 s vs 5.96 s median (+0.0%) |
| wall, 250 k pairs, 1 thread, 7 runs each (`times_t1.tsv`) | 8.01 s vs 8.13 s median (+1.5%, ranges overlap: 7.93-8.39 / 7.95-8.29) |

Per stage, from flexalign's own timers, 1 thread, 250 k pairs (median of 9 runs compact, 7 hashed;
`stages_t1*.tsv`). The two stages the change touches are key lookup ("getting ranges") and value
reads ("getting range headers"):

| stage | compact base | compact flags | Δ | hashed base | hashed flags | Δ |
| --- | --- | --- | --- | --- | --- | --- |
| getting ranges | 1.350 s | 1.334 s | -1.2% | 1.027 s | 1.033 s | +0.6% |
| getting range headers | 0.911 s | 0.898 s | -1.4% | 0.894 s | 0.888 s | -0.6% |
| getting kmers (untouched) | 1.854 s | 1.849 s | -0.3% | 1.844 s | 1.843 s | -0.0% |
| hamming screen (untouched) | 1.069 s | 1.045 s | -2.2% | 1.035 s | 1.024 s | -1.0% |

Every difference lies within the run-to-run spread, including on stages the change does not touch.
`perf stat` was not available (`perf_event_paranoid` = 4), so instruction counts were not compared.

## Why no cost is expected

An index without flags has `flag_bits = 0`, so a packed value's stride and read width are both
`ref_bits + pos_bits`, as before; the read masks with the value width instead of the stride (the
same instructions). `PackedRange` keeps its size (two `u32` fields became `u32` + 2 × `u16`). Blob
cells: `VD::get` already masked the reference id, so flags in bits 62-63 are invisible to it. The
file is still written as packed version 1.

## Re-run

    python3 make_data.py ~/.gtdb_downloader/r232/g__Bacteroides.map.tsv <dir>/data 1000000
    # binaries in <dir>/bin/flexalign_{base,flags}, references linked in <dir>/idx_{base,flags}/
    bench.sh <dir> 9 8 <dir>/data/reads_R1.fq <dir>/data/reads_R2.fq t8
    stages.sh <dir> 9 <R1 250k> <R2 250k>
    EXTRA="--index-format hashed" TAG=_hashed stages.sh <dir> 7 <R1> <R2>
    python3 summarize.py <dir>/times_*.tsv
