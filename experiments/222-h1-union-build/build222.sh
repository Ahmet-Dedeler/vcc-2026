#!/bin/bash
# 222 — `cd4f034`'s configuration with ONE variable: the coverage-aware sizing rule and the
# H1 union residual it unlocks.
#
# This is `experiments/210-cd4-build/build210.sh` with exactly ONE line added --
# `export VCC_H1_UNION_G="${VCC_H1_UNION_G:-}"` -- so the arm's dose is the only thing that
# distinguishes it and the KNOWN-ANSWER ROW re-runs the incumbent through this same script
# with the flag empty.  Every other flag, every other environment variable and the emitter
# are byte-identical to 210's, which is byte-identical to 202's and 191's:
#   alpha 0.982 | cd4 f 0.34 | min_cpm 10 | tail_spread 1.0 | keep_top 8000
#   keep_spread 0.25 | overdispersion 4.0 | whitener_id union_pruned | emitter rc_callset
#
# usage: build222.sh <out.h5ad>
#
#   build222.sh submissions/222-ka-cd4f034-f034.h5ad
#                                  KNOWN-ANSWER ROW 1: with VCC_H1_UNION_G unset its sha256
#                                  must be
#                                  7036e89f8a03c82d69e79a5c6690b275adae52dbb3ff72bb6ab23e7ca622bf7a
#                                  on 12,648,807,800 bytes, or the flag-off path is not the
#                                  path that produced `cd4f034` and nothing after it reads.
#
#   VCC_H1_UNION_G=0.28 build222.sh submissions/222-h1union-g028.h5ad      the arm.
set -euo pipefail
cd "$(dirname "$0")/../.."
OUT="$1"
export VCC_ALPHA=0.982
export VCC_CD4_TARGET_F=0.34
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
