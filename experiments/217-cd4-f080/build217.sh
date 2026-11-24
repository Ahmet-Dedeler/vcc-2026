#!/bin/bash
# 217 — `cd4f034`'s configuration, built locally, with ONE variable: the CD4 dose.
#
# This is `experiments/210-cd4-build/build210.sh` with exactly ONE executable line changed —
# the default of `VCC_CD4_TARGET_F` becomes 0.80 instead of 0.34 — so the arm's dose is the
# default and the KNOWN-ANSWER ROW can re-run the incumbent's 0.34 through this same script.
# Every other flag, every other environment variable and the emitter are byte-identical to
# 210's, which is byte-identical to 202's, which is byte-identical to 191's:
#   alpha 0.982 | min_cpm 10 | tail_spread 1.0 | keep_top 8000 | keep_spread 0.25
#   overdispersion 4.0 | whitener_id union_pruned | spread_map block | effect_scale 1.0
#   no shrink | no mass_fix | donor weights all 1.0 | emitter rc_callset
#
# Set `f`, never `fraction` (207 §4 fn 1): the builder solves f = 0.80 to fraction ~ 0.758 on
# the challenge panel.  `VCC_CD4_FRACTION` is the legacy path and stays at its 0.05 default.
# `VCC_ALPHA` is NOT re-solved — 207 §6.2 measured alpha* = 0.982 exactly, 211 re-solved it at
# 0.9808/0.9879 (inside noise); re-solving would make this two variables.
#
# usage: build217.sh <arm> [output.h5ad]
#
#   VCC_CD4_TARGET_F=0.34 build217.sh rc_callset submissions/217-ka-cd4f034-f034.h5ad
#                                  KNOWN-ANSWER ROW: its sha256 must be
#                                  7036e89f8a03c82d69e79a5c6690b275adae52dbb3ff72bb6ab23e7ca622bf7a
#                                  on 12,648,807,800 bytes, or this environment differs from
#                                  210's and nothing after it is readable.
#
#   build217.sh rc_callset submissions/217-cd4f080.h5ad         the arm, at the default f = 0.80.
#
# (The output path does not reach the emitted bytes: 191 and 202 wrote `rc_hybrid` to two
#  different paths and got the same sha256, e0b0677f…; 210 re-proved it on `rc_callset`.)
set -euo pipefail
cd "$(dirname "$0")/../.."
ARM="$1"
OUT="${2:-submissions/217-${ARM}.h5ad}"
EXTRA=()
case "$ARM" in
  a15_local) ;;
  rc_hybrid) EXTRA=(--emitter rc_hybrid) ;;
  rc_callset) EXTRA=(--emitter rc_callset) ;;
  *) echo "unknown arm $ARM" >&2; exit 2 ;;
esac
export VCC_ALPHA=0.982
export VCC_CD4_TARGET_F="${VCC_CD4_TARGET_F:-0.80}"
exec /usr/bin/time -l uv run python scripts/make_submission_pm_cd4.py \
  --whitener data/donors/whitener_union_pruned.npz \
  --whitener-metadata data/donors/whitener_union_pruned.json \
  --whitener-id union_pruned \
  --min-cpm 10 \
  --tail-spread 1.0 \
  --keep-top 8000 \
  --keep-spread 0.25 \
  --overdispersion 4.0 \
  --max-rss-gib 24 \
  ${EXTRA[@]+"${EXTRA[@]}"} \
  -o "$OUT"
