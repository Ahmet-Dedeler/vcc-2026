#!/bin/bash
# 255 -- `experiments/249-cd4-deflate/build249.sh` (= 222's = 210's = 202's = 191's flags, byte for byte)
# with the shipped VCC_H1_UNION_G=0.28 pinned and ONE passthrough: VCC_CD4_DEFLATE_K.  Every other
# inert flag the builder has grown since (247/248/251/254) is UNSET explicitly.
#
#   alpha 0.982 | cd4 f 0.34 | h1 union g 0.28 | min_cpm 10 | tail_spread 1.0 | keep_top 8000
#   keep_spread 0.25 | overdispersion 4.0 | whitener_id union_pruned | emitter rc_callset
#
# usage: VCC_CD4_DEFLATE_K=10 build255.sh <out.h5ad>
#   VCC_CD4_DEFLATE_K=10 -> the arm, submissions/255-cd4svd10-f034.h5ad
#   (flag unset would rebuild the shipped h1union 8df4dadd...; NOT run here -- 254 section KA 1 proved it
#    through these exact builder bytes 3615d20f... an hour before this experiment started)
set -euo pipefail
cd "$(dirname "$0")/../.."
OUT="$1"
export VCC_ALPHA=0.982
export VCC_CD4_TARGET_F=0.34
export VCC_H1_UNION_G=0.28
export VCC_CD4_DEFLATE_K="${VCC_CD4_DEFLATE_K:-}"
unset VCC_CD4_SIZING VCC_CD4_TARGET_G VCC_UNION_SOURCES VCC_CD4_DEFLATE_BASIS VCC_UNION2_SOURCES VCC_UNION2_G VCC_CD4_PARTIAL_G
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
