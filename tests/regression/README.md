# Regression indicator

A lightweight way to see whether a change **improves, worsens, or keeps neutral**
the alignment results, without committing the multi-GB reference/index.

## How it works

`run.sh` aligns a fixed subset of `data/paired_dat{1,2}.fq` (default: first 1000
read pairs) single-threaded for a given `--profile`, normalizes stdout to just the
sorted PAF records, and diffs it against a committed baseline.

- **Committed:** `run.sh`, `baseline.<profile>.paf` (the current golden output, ~160 KB),
  and `history/` (every superseded baseline + a log of all updates).
- **Not committed:** the reference FASTA, the `.flex.*` index, and the reads — they
  live under `data/` (gitignored) and must be present on the machine. The subset is
  derived deterministically with `head`, so no read data needs to be committed.
- Output is single-threaded (`-t 1`) and sorted, so it is deterministic and
  order-independent. Non-PAF lines (logs, stray debug prints) are filtered out.

## Usage

```bash
cargo +nightly build --release

# Compare current output against the baseline (exit 0 = neutral, 1 = changed):
./tests/regression/run.sh

# A different profile or a larger subset:
PROFILE=mid ./tests/regression/run.sh
NREADS=5000 ./tests/regression/run.sh

# After an *intended* change (you've reviewed the diff and it's an improvement),
# refresh the baseline so future runs compare against the new expected output.
# This ARCHIVES the old baseline and logs the update (optionally with a note):
./tests/regression/run.sh --update "switched to revised anchor extractor"
just regression-update "switched to revised anchor extractor"        # via just
```

## Workflow for a change

1. Run `./tests/regression/run.sh` and confirm `NEUTRAL` on a clean tree (baseline valid).
2. Make your change, rebuild.
3. Run it again:
   - `NEUTRAL` → results unchanged.
   - `CHANGED: +N / -M` → inspect `.work/diff.<profile>.txt`. Fewer/"better" mappings
     vs more/"worse" is a judgement call for this dataset; the diff is the evidence.
4. If the change is intended, commit it **together with** a `--update`d baseline so the
   new expectation is recorded. Reviewers see the result delta in the baseline diff.

## History of past results

Every `--update` keeps a full record, so no old result is ever lost:

- `history/baseline.<profile>.<timestamp>.paf` — a snapshot of each superseded baseline.
- `history/log.tsv` — one row per update: UTC timestamp, profile, git SHA, record
  count, `+added`/`-removed` vs the previous baseline, and your note. View it with
  `column -t -s$'\t' tests/regression/history/log.tsv`.

The active baseline is `baseline.<profile>.paf`; the full set of results over time is
that file plus everything under `history/`. To compare the current build against an
older result, diff its archived snapshot against fresh output, e.g.:

```bash
diff <(LC_ALL=C sort tests/regression/history/baseline.default.<timestamp>.paf) \
     tests/regression/.work/current.default.paf
```

> The baseline only captures *that* output changed, not whether it's biologically
> more correct. For true accuracy, build with `FLEXALIGN_GOLDSTD_EVAL=1` (the reads'
> headers encode their true origin) and compare the gold-standard stats this emits.
> This harness is the fast, always-on indicator; gold-standard eval is the deeper check.
