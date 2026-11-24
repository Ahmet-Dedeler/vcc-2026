#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
STATE_DIR="$ROOT/data/external/state-vcc/state"
FROZEN_SUPPORT_DIR="$ROOT/data/external/state-vcc/competition_support_set"
SUPPORT_DIR="$ROOT/data/external/state-vcc/competition_support_set-log1p"
RUN_ROOT="$ROOT/experiments/106-state-baseline/runs"
LOCAL_TOML="$RUN_ROOT/starter.local.toml"

PROFILE="${1:-state}"
MAX_STEPS="${2:-20000}"

case "$PROFILE" in
  state|state_sm) ;;
  *)
    echo "profile must be state or state_sm" >&2
    exit 2
    ;;
esac

if [[ ! -f "$SUPPORT_DIR/starter.toml" ]]; then
  echo "missing prepared support view: $SUPPORT_DIR/starter.toml" >&2
  exit 2
fi

mkdir -p "$RUN_ROOT"
python3 "$ROOT/experiments/106-state-baseline/prepare_local_config.py" \
  "$FROZEN_SUPPORT_DIR/starter.toml" "$SUPPORT_DIR" "$LOCAL_TOML"
cd "$STATE_DIR"

export STATE_USE_MPS=1
export PYTORCH_ENABLE_MPS_FALLBACK=0

exec python3 "$ROOT/scripts/run_with_memory_guard.py" \
  --max-rss-gib 44 \
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
  training.ckpt_every_n_steps=20000 \
  "model=$PROFILE" \
  use_wandb=false \
  "output_dir=$RUN_ROOT" \
  "name=official-$PROFILE"
