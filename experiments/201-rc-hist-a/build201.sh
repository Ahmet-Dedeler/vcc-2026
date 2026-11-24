#!/bin/bash
# 201 — a15's configuration, built locally, with ONE variable: the emitter.
#
# Byte-for-byte the same invocation 191 used (`experiments/191-rc-hybrid-build/build191.sh`),
# with `rc_hybrid` replaced by `rc_hist_a`.  a15's configuration is read off
# `submissions/a15-cd4f010.vcc.meta.json` field by field:
#   alpha 0.982 | cd4_target_f 0.1 | min_cpm 10 | tail_spread 1.0 | keep_top 8000
#   keep_spread 0.25 | overdispersion 4.0 | whitener_id union_pruned | spread_map block
#   effect_scale 1.0 | no shrink | no mass_fix | donor weights all 1.0
# Everything else is the builder's shipped default, which is what a15 shipped too.
#
# usage: build201.sh rc_hist_a          (--emitter rc_hist_a, nothing else)
#        build201.sh rc_hybrid          (191's arm, for a rebuild check)
#        build201.sh a15_local          (the shipped quasi emitter)
set -euo pipefail
cd "$(dirname "$0")/../.."
ARM="$1"
EXTRA=()
case "$ARM" in
  a15_local) ;;
  rc_hybrid) EXTRA=(--emitter rc_hybrid) ;;
  rc_hist_a) EXTRA=(--emitter rc_hist_a) ;;
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
  -o "submissions/201-${ARM}.h5ad"
