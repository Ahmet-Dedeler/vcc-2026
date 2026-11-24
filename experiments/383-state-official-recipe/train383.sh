#!/usr/bin/env bash
# 383 — official-recipe training launcher (mirrors 106 run_training.sh, model=state).
# Usage: train383.sh <max_steps> [run_name]
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
HERE="$ROOT/experiments/383-state-official-recipe"
STATE_DIR="$ROOT/data/external/state-vcc/state"
SUPPORT_SRC="$HERE/data/competition_support_set"
SUPPORT_DIR="$HERE/data/competition_support_set-log1p"
RUN_ROOT="$HERE/runs"
LOCAL_TOML="$RUN_ROOT/starter.local.toml"

MAX_STEPS="${1:?max_steps required}"
NAME="${2:-official-state}"
CKPT_EVERY="${CKPT_EVERY:-20000}"

if [[ ! -f "$SUPPORT_DIR/starter.toml" ]]; then
  echo "missing prepared support view: $SUPPORT_DIR/starter.toml" >&2
  exit 2
fi

# Memory gate per protocol: free >= 40%, and yield to the PI build.
FREE_PCT="$(memory_pressure | grep -o 'System-wide memory free percentage: [0-9]*%' | grep -o '[0-9]*')"
if [[ "$FREE_PCT" -lt 40 ]]; then
  echo "memory gate: free ${FREE_PCT}% < 40%, refusing" >&2
  exit 3
fi
if pgrep -f def_pipeline >/dev/null; then
  echo "PI build running, refusing" >&2
  exit 3
fi

mkdir -p "$RUN_ROOT"
python3 "$HERE/prepare_local_config.py" \
  "$SUPPORT_SRC/starter.toml" "$SUPPORT_DIR" "$LOCAL_TOML"
cd "$STATE_DIR"

export STATE_USE_MPS=1
export PYTORCH_ENABLE_MPS_FALLBACK=0
export NUMBA_NUM_THREADS=4 OMP_NUM_THREADS=4

exec python3 "$ROOT/scripts/run_with_memory_guard.py" \
  --max-rss-gib 40 \
  --min-free-percent 15 \
  --poll-seconds 5 \
  -- \
  "$STATE_DIR/.venv/bin/state" tx train \
  "data.kwargs.toml_config_path=$LOCAL_TOML" \
  data.kwargs.num_workers=8 \
  data.kwargs.pin_memory=false \
  data.kwargs.batch_col=batch_var \
  data.kwargs.pert_col=target_gene \
  data.kwargs.cell_type_key=cell_type \
  data.kwargs.control_pert=non-targeting \
  "data.kwargs.perturbation_features_file=$SUPPORT_DIR/ESM2_pert_features.pt" \
  "training.max_steps=$MAX_STEPS" \
  "training.ckpt_every_n_steps=$CKPT_EVERY" \
  model=state \
  use_wandb=false \
  "output_dir=$RUN_ROOT" \
  "name=$NAME"
