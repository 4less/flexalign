#!/usr/bin/env bash
# Per-stage times from flexalign's log, interleaved runs of both binaries, single-threaded.
# usage: stages.sh <dir> <reps> <R1> <R2>
set -euo pipefail
dir=$1; reps=$2; r1=$3; r2=$4
out=$dir/stages_t1${TAG:-}.tsv
echo -e "rep\tbinary\tstage\tseconds" > "$out"
for rep in $(seq 1 "$reps"); do
  if (( rep % 2 )); then order="base flags"; else order="flags base"; fi
  for v in $order; do
    "$dir/bin/flexalign_$v" -r "$dir/idx_$v/reference.fna" -1 "$r1" -2 "$r2" -t 1 ${EXTRA:-} 2>&1 >/dev/null \
      | sed 's/\x1b\[[0-9;]*m//g' | grep -E '^Time for [a-z ]+\.{3,}[0-9.]+(ns|µs|ms|s)$' \
      | sed -E 's/^Time for ([a-z ]+)\.+([0-9.]+)(ns|µs|ms|s)$/\1\t\2\t\3/' \
      | awk -F'\t' -v r=$rep -v b=$v 'BEGIN{OFS="\t"} {m=($3=="s")?1:($3=="ms")?1e-3:($3=="µs")?1e-6:1e-9; print r,b,$1,$2*m}' >> "$out"
  done
done
