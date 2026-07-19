#!/usr/bin/env bash
#
# Regression indicator for flexalign alignment output.
#
# Runs the release binary on a fixed, deterministically-derived read subset for a
# given --profile, normalizes stdout to just the PAF records (sorted by read id,
# stray non-PAF lines dropped), and diffs it against a committed baseline.
#
#   ./tests/regression/run.sh            compare current output vs baseline
#   ./tests/regression/run.sh --update   (re)write the baseline from current output
#   PROFILE=mid ./tests/regression/run.sh        check a different profile
#   NREADS=5000 ./tests/regression/run.sh        use a larger subset
#
# Exit codes:  0 = neutral (matches baseline)   1 = changed (regression indicator)
#              2 = setup error
#
# Inputs (reference, index, reads) live under data/ and are machine-local
# (data/ is gitignored). Only the small normalized baseline + this script are
# committed. The read subset is derived with `head`, so it is reproducible given
# the same data/ files.
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
REG="$ROOT/tests/regression"
DATA="$ROOT/data"
REF="$DATA/GCA_000006155.2_ASM615v2_genomic.fna"
BIN="$ROOT/target/release/flexalign"

PROFILE="${PROFILE:-default}"
NREADS="${NREADS:-1000}"          # number of read *pairs*
NLINES=$((NREADS * 4))
BASELINE="$REG/baseline.$PROFILE.paf"
WORK="$REG/.work"
CURRENT="$WORK/current.$PROFILE.paf"
HIST="$REG/history"
LOG="$HIST/log.tsv"
UPDATE=0
NOTE=""
[ "${1:-}" = "--update" ] && { UPDATE=1; NOTE="${2:-}"; }

git_sha() { git -C "$ROOT" rev-parse --short HEAD 2>/dev/null || echo nogit; }

[ -x "$BIN" ] || { echo "Build first:  cargo +nightly build --release" >&2; exit 2; }
[ -f "$REF" ] || { echo "Missing reference: $REF" >&2; exit 2; }
[ -f "$DATA/paired_dat1.fq" ] || { echo "Missing reads: $DATA/paired_dat1.fq" >&2; exit 2; }

mkdir -p "$WORK"
head -n "$NLINES" "$DATA/paired_dat1.fq" > "$WORK/r1.fq"
head -n "$NLINES" "$DATA/paired_dat2.fq" > "$WORK/r2.fq"

# -t 1 => deterministic record order. Keep only PAF records (they contain tabs),
# then sort so the check is independent of output order and free of log/debug noise.
"$BIN" -r "$REF" -1 "$WORK/r1.fq" -2 "$WORK/r2.fq" -t 1 --profile "$PROFILE" --log-level 1 \
    2> "$WORK/stderr.$PROFILE.log" \
  | grep -P '\t' | LC_ALL=C sort > "$CURRENT"

records=$(wc -l < "$CURRENT")

if [ "$UPDATE" = 1 ] || [ ! -f "$BASELINE" ]; then
    mkdir -p "$HIST"
    [ -f "$LOG" ] || printf 'timestamp_utc\tprofile\tgit_sha\trecords\tadded_vs_prev\tremoved_vs_prev\tnote\n' > "$LOG"

    ts=$(date -u +%Y%m%dT%H%M%SZ)
    added=0
    removed=0

    # Archive the outgoing baseline (the one being replaced) so old results are
    # never lost, and record how the new output differs from it.
    if [ -f "$BASELINE" ]; then
        cp "$BASELINE" "$HIST/baseline.$PROFILE.$ts.paf"
        if ! diff -u "$BASELINE" "$CURRENT" > "$WORK/upd.$PROFILE.diff"; then
            added=$(grep -c '^+[^+]' "$WORK/upd.$PROFILE.diff" || true)
            removed=$(grep -c '^-[^-]' "$WORK/upd.$PROFILE.diff" || true)
        fi
        echo "Previous baseline archived: history/baseline.$PROFILE.$ts.paf"
    fi

    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
        "$ts" "$PROFILE" "$(git_sha)" "$records" "$added" "$removed" "$NOTE" >> "$LOG"

    cp "$CURRENT" "$BASELINE"
    echo "Baseline updated: $BASELINE ($records records, +$added/-$removed vs previous) [profile=$PROFILE]"
    echo "History log: history/log.tsv"
    exit 0
fi

if diff -u "$BASELINE" "$CURRENT" > "$WORK/diff.$PROFILE.txt"; then
    echo "NEUTRAL: output matches baseline ($records records) [profile=$PROFILE]"
    exit 0
else
    added=$(grep -c '^+[^+]' "$WORK/diff.$PROFILE.txt" || true)
    removed=$(grep -c '^-[^-]' "$WORK/diff.$PROFILE.txt" || true)
    echo "CHANGED [profile=$PROFILE]: +$added new / -$removed gone vs baseline ($records records now)"
    echo "Full diff: $WORK/diff.$PROFILE.txt"
    echo "If the change is intended/an improvement, refresh with: $0 --update"
    exit 1
fi
