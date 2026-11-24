#!/bin/bash
# 229 — `cd4f034`'s configuration with ONE variable: which banks the coverage-aware second
# residual unions, and its dose.
#
# This is `experiments/222-h1-union-build/build222.sh` with exactly ONE line added --
# `export VCC_UNION_SOURCES="${VCC_UNION_SOURCES:-}"` -- so the arm's source list is the only
# thing that distinguishes it from the shipped `h1union`, and the two KNOWN-ANSWER ROWS re-run
# the incumbent and the shipped arm through this same script.  Every other flag, every other
# environment variable and the emitter are byte-identical to 222's, which is 210's, which is
# 202's and 191's:
#   alpha 0.982 | cd4 f 0.34 | min_cpm 10 | tail_spread 1.0 | keep_top 8000
#   keep_spread 0.25 | overdispersion 4.0 | whitener_id union_pruned | emitter rc_callset
#
# usage: build229.sh <out.h5ad>
#
#   build229.sh submissions/229-ka-cd4f034-f034.h5ad
#                                  KNOWN-ANSWER ROW 1: with both flags unset its sha256 must be
#                                  7036e89f8a03c82d69e79a5c6690b275adae52dbb3ff72bb6ab23e7ca622bf7a
#                                  on 12,648,807,800 bytes, or the flag-off path is not the path
#                                  that produced `cd4f034` and nothing after it reads.
#
#   VCC_H1_UNION_G=0.28 build229.sh submissions/229-ka-h1union-g028.h5ad
#                                  KNOWN-ANSWER ROW 2: VCC_UNION_SOURCES unset selects the H1
#                                  pair, so this must be
#                                  8df4dadd5a63b76c7ca3e3d57a511e7e52000b79d016d140897c55dd24c3ecf6
#                                  and its `.vcc`
#                                  f13ce40bca908af7cb85b8e84887af473aae4588410c238529139b7d1f245b1e
#                                  on 2,824,785,920 bytes -- the shipped `h1union` package.
#
#   VCC_H1_UNION_G=0.35 \
#   VCC_UNION_SOURCES=dld1_lowmoi,jiang_5line_noK562,igvf_cardiomyocyte,myllia_a549 \
#     build229.sh submissions/229-stack4-g035.h5ad                              the arm.
set -euo pipefail
cd "$(dirname "$0")/../.."
OUT="$1"
export VCC_ALPHA=0.982
export VCC_CD4_TARGET_F=0.34
export VCC_UNION_SOURCES="${VCC_UNION_SOURCES:-}"
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
