#!/bin/bash
# 202 — a15's configuration, built locally, with ONE variable: the emitter.
#
# Byte-for-byte the same invocation 191 used (`experiments/191-rc-hybrid-build/build191.sh`)
# and 201 reused, with `rc_hybrid` replaced by `rc_callset`.  a15's configuration is read off
# `submissions/a15-cd4f010.vcc.meta.json` field by field:
#   alpha 0.982 | cd4_target_f 0.1 | min_cpm 10 | tail_spread 1.0 | keep_top 8000
#   keep_spread 0.25 | overdispersion 4.0 | whitener_id union_pruned | spread_map block
#   effect_scale 1.0 | no shrink | no mass_fix | donor weights all 1.0
# Everything else is the builder's shipped default, which is what a15 shipped too.
#
# usage: build202.sh rc_hybrid          191's arm, rebuilt through the PATCHED builder.
#                                       KNOWN-ANSWER ROW 1: its sha256 must be
#                                       e0b0677f984a94a1c4a7af2add196908cdab90d7acb463eb1d6e0397c2d870a8
#                                       on 13,202,460,120 bytes, or the patch moved the
#                                       incumbent and the build stops.
#        build202.sh rc_callset         the arm.
set -euo pipefail
cd "$(dirname "$0")/../.."
ARM="$1"
EXTRA=()
case "$ARM" in
  a15_local) ;;
  rc_hybrid) EXTRA=(--emitter rc_hybrid) ;;
  rc_callset) EXTRA=(--emitter rc_callset) ;;
  *) echo "unknown arm $ARM" >&2; exit 2 ;;
esac
export VCC_ALPHA=0.982
export VCC_CD4_TARGET_F=0.1
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
  -o "submissions/202-${ARM}.h5ad"
