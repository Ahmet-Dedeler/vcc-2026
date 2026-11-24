#!/usr/bin/env bash
# 438: rehearse the day-of command end to end on 238's essential-free pessimistic bundle (efp), timing it.
cd "$(dirname "$0")"
export PYTHONDONTWRITEBYTECODE=1
R=~/Code/virtual-cell-challenge-2026/experiments
echo "$(date -u +%FT%TZ) start dayof rehearsal (efp bundle, prod fold operators)" >> work/chain438.log
/usr/bin/time -l uv run python dayof438.py --bundle $R/238-essential-free-panel/work/ef-bundle \
  --op $R/244-october-rehearsal/work/run1/whitener_union_pruned_def.npz \
  --cd4-root $R/238-essential-free-panel/work/cd4/efp/cd4_compact --tag dayof_efp --synthetic > work/logs/dayof_efp.log 2>&1
echo "$(date -u +%FT%TZ) end dayof rehearsal rc=$?" >> work/chain438.log
