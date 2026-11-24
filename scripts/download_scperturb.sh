#!/usr/bin/env bash
# scPerturb (Peidli et al., Nat Methods 2024): 50 harmonized single-cell perturbation
# datasets, ~40 GB. Standardized .obs schema across studies, which is what makes them
# usable as a single cross-context training corpus.
#
# Source: https://zenodo.org/records/10044268  ·  http://projects.sanderlab.org/scperturb/
# Smallest-first so cheap datasets land early.
set -uo pipefail
DEST="data/external/scperturb"
LIST="scripts/scperturb_files.tsv"
mkdir -p "$DEST"

while IFS=$'\t' read -r key size url; do
  out="$DEST/$key"
  if [ -f "$out" ]; then
    have=$(stat -f%z "$out" 2>/dev/null || echo 0)
    [ "$have" = "$size" ] && { echo "have  $key"; continue; }
    echo "resize $key ($have != $size), refetching"
  fi
  printf 'get   %-52s %6.2f GB\n' "$key" "$(echo "$size/1000000000" | bc -l)"
  if curl -fL --retry 5 --retry-delay 5 -C - -s -S -o "$out.part" "$url"; then
    mv "$out.part" "$out"
  else
    echo "FAIL  $key"
  fi
done < "$LIST"
echo "---"; du -sh "$DEST"
