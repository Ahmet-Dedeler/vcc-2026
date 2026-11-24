#!/bin/bash
# 248 -- `experiments/222-h1-union-build/build222.sh` (= 210's = 202's = 191's flags, byte for byte)
# with the shipped VCC_H1_UNION_G=0.28 pinned and TWO passthroughs: VCC_UNION2_SOURCES, VCC_UNION2_G.
#
#   alpha 0.982 | cd4 f 0.34 | h1 union g 0.28 | min_cpm 10 | tail_spread 1.0 | keep_top 8000
#   keep_spread 0.25 | overdispersion 4.0 | whitener_id union_pruned | emitter rc_callset
#
# usage: [VCC_UNION2_SOURCES=dld1_lowmoi VCC_UNION2_G=0.35] build248.sh <out.h5ad>
#   both unset -> KNOWN-ANSWER ROW 1: sha256 8df4dadd5a63b76c7ca3e3d57a511e7e52000b79d016d140897c55dd24c3ecf6
#                 on 12,648,807,800 bytes (the shipped h1union), or the default path is not inert.
#   VCC_UNION2_SOURCES=dld1_lowmoi VCC_UNION2_G=0.35 -> the rider, submissions/248-cum1-dld1-g035.h5ad
set -euo pipefail
cd "$(dirname "$0")/../.."
OUT="$1"
export VCC_ALPHA=0.982
export VCC_CD4_TARGET_F=0.34
export VCC_H1_UNION_G=0.28
export VCC_UNION2_SOURCES="${VCC_UNION2_SOURCES:-}"
export VCC_UNION2_G="${VCC_UNION2_G:-}"
unset VCC_CD4_SIZING VCC_CD4_TARGET_G VCC_UNION_SOURCES VCC_CD4_DEFLATE_K
exec /usr/bin/time -l uv run python scripts/make_submission_pm_cd4.py \
  --whitener data/donors/whitener_union_pruned.npz \
  --whitener-metadata data/donors/whitener_union_pruned.json \
  --whitener-id union_pruned \
  --min-cpm 10 \
  --tail-spread 1.0 \
  --keep-top 8000 \
  --keep-spread 0.25 \
  --overdispersion 4.0 \
  --emitter rc_callset \
  --max-rss-gib 24 \
  -o "$OUT"
