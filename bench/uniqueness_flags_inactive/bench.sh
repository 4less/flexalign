#!/usr/bin/env bash
# Interleaved timing of two flexalign binaries on the same reads: base, flags, base, flags, ...
# usage: bench.sh <dir> <reps> <threads> <R1> <R2> <tag>
set -euo pipefail
dir=$1; reps=$2; t=$3; r1=$4; r2=$5; tag=$6
out=$dir/times_$tag.tsv
echo -e "rep\tbinary\tthreads\twall_s\tuser_s\tsys_s\tmaxrss_kb" > "$out"
for v in base flags; do cat "$dir"/idx_$v/*.compact "$dir"/idx_$v/*.refseq > /dev/null; done
for rep in $(seq 1 "$reps"); do
  # Alternate which binary goes first, so neither always runs right after the other.
  if (( rep % 2 )); then order="base flags"; else order="flags base"; fi
  for v in $order; do
    /usr/bin/time -f "%e\t%U\t%S\t%M" -o "$dir/.t" "$dir/bin/flexalign_$v" -r "$dir/idx_$v/reference.fna" \
      -1 "$r1" -2 "$r2" -t "$t" > "$dir/out_$v/$tag.paf" 2> "$dir/out_$v/$tag.log"
    echo -e "$rep\t$v\t$t\t$(cat "$dir/.t")" >> "$out"
  done
done
cmp -s <(sort "$dir/out_base/$tag.paf") <(sort "$dir/out_flags/$tag.paf") && echo "$tag: PAF identical (sorted; threads emit in any order)" || echo "$tag: PAF DIFFERS"
