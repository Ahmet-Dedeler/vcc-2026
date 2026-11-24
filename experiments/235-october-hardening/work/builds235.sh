#!/bin/bash
# 235 known-answer rows 1 and 2, through 222's byte-exact build script.
set -euo pipefail
cd "$(dirname "$0")/../../.."
W=experiments/235-october-hardening/work
experiments/222-h1-union-build/build222.sh submissions/235-ka-cd4f034-f034.h5ad > "$W/build-ka.log" 2>&1
VCC_H1_UNION_G=0.28 experiments/222-h1-union-build/build222.sh submissions/235-h1union-g028.h5ad > "$W/build-arm.log" 2>&1
/usr/bin/time -l uv run vcc prep submissions/235-h1union-g028.h5ad \
  -g data/challenge/gene_names.csv --perts data/challenge/pert_counts.csv \
  -o submissions/235-h1union-g028.vcc > "$W/prep.log" 2>&1
shasum -a 256 submissions/235-ka-cd4f034-f034.h5ad submissions/235-h1union-g028.h5ad \
  submissions/235-h1union-g028.vcc > "$W/hashes.txt"
wc -c submissions/235-ka-cd4f034-f034.h5ad submissions/235-h1union-g028.h5ad \
  submissions/235-h1union-g028.vcc >> "$W/hashes.txt"
cat "$W/hashes.txt"
