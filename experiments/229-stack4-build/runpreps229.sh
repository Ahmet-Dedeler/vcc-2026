#!/bin/bash
# 229 — the preps, STRICTLY SERIAL and never alongside a build (`vcc prep` peaks at 22.33 GiB).
# Dry-run first (mandate falsifier 5), then the real package.
set -uo pipefail
cd "$(dirname "$0")/../.."
HERE=experiments/229-stack4-build
G=data/challenge/gene_names.csv
P=data/challenge/pert_counts.csv

for tag in ka-h1union-g028 stack6-g035 stack4-g035; do
  h5=submissions/229-$tag.h5ad
  vcc=submissions/229-$tag.vcc
  echo "=== $tag dry-run $(date -u +%FT%TZ) ==="
  /usr/bin/time -l vcc prep "$h5" -g "$G" --perts "$P" --dry-run \
      > "$HERE/work/prep-$tag-dryrun.log" 2>&1
  echo "  dry-run exit=$?"
  echo "=== $tag prep $(date -u +%FT%TZ) ==="
  /usr/bin/time -l vcc prep "$h5" -g "$G" --perts "$P" -o "$vcc" \
      > "$HERE/work/prep-$tag.log" 2>&1
  echo "  prep exit=$?"
  if [ -f "$vcc" ]; then
    echo -n "  $tag vcc sha256 "; shasum -a 256 "$vcc" | awk '{print $1}'
    echo -n "  $tag vcc bytes  "; stat -f %z "$vcc"
  fi
  grep -E "maximum resident set size" "$HERE/work/prep-$tag.log" | tail -1
done
echo "=== PREPS DONE $(date -u +%FT%TZ) ==="
