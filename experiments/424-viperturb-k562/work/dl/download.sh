#!/bin/bash
# 424: fetch VIPerturb-seq genome-wide bins + manifest from Zenodo 18460279 (CC BY 4.0)
set -u
base=https://zenodo.org/api/records/18460279/files
for f in genome_wide_manifest.txt genome_wide_binC.RDS genome_wide_binA.RDS genome_wide_binB.RDS; do
  echo "$(date -u +%FT%TZ) start $f"
  curl -sS -L --retry 5 --retry-delay 10 -C - -o "$f" "$base/$f/content"
  echo "$(date -u +%FT%TZ) done $f rc=$? size=$(stat -f%z "$f")"
done
md5 genome_wide_manifest.txt genome_wide_bin*.RDS
echo ALLDONE
