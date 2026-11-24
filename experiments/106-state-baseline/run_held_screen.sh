#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
HERE="$ROOT/experiments/106-state-baseline"
STATE_PY="$ROOT/data/external/state-vcc/state/.venv/bin/python"
GUARD="$ROOT/scripts/run_with_memory_guard.py"
FROZEN_CHECKPOINT="$HERE/frozen/state_sm-step20000"
RUN_ID="${EXP106_RUN_ID:-$(date -u +%Y%m%dT%H%M%SZ)}"
RESUME="${EXP106_RESUME:-0}"
if [[ ! "$RUN_ID" =~ ^[A-Za-z0-9][A-Za-z0-9._-]{0,127}$ ]]; then
  echo "invalid EXP106_RUN_ID: $RUN_ID" >&2
  exit 2
fi
RUN_ROOT="$HERE/work/runs/$RUN_ID"
WORK="$RUN_ROOT/a14"
CONTROL_DIR="$RUN_ROOT/controls"
INCUMBENT_DIR="$RUN_ROOT/incumbent-a14"
LOCK_DIR="$HERE/work/.held-screen.lock"

mkdir -p "$HERE/work"
if ! mkdir "$LOCK_DIR" 2>/dev/null; then
  echo "another experiment-106 held screen holds $LOCK_DIR" >&2
  exit 2
fi
cleanup() {
  rmdir "$LOCK_DIR" 2>/dev/null || true
}
trap cleanup EXIT

assert_no_other_heavy() {
  if pgrep -f \
    '[c]ell-eval2|[s]tate tx train|[s]tate_effect_adapter.py|[p]repare_held_pruned_operators.py|[b]uild_union_whitener.py|[r]un_with_memory_guard.py' \
    >/dev/null; then
    echo "another heavy VCC process is active; refusing to continue." >&2
    exit 2
  fi
}

if pgrep -f "state tx train.*official-state_sm" >/dev/null; then
  echo "STATE training is still active; refusing to start the held screen." >&2
  exit 2
fi
assert_no_other_heavy
if [[ -e "$RUN_ROOT" && "$RESUME" != "1" ]]; then
  echo "evaluation run already exists: $RUN_ROOT" >&2
  exit 2
fi
mkdir -p "$RUN_ROOT"
echo "experiment-106 evaluation run: $RUN_ID (resume=$RESUME)"

"$STATE_PY" "$HERE/freeze_checkpoint.py"

mkdir -p "$CONTROL_DIR"
for spec in "hct116:HCT116" "hek293t:HEK293T"; do
  surface="${spec%%:*}"
  line="${spec##*:}"
  controls="$CONTROL_DIR/orion-${surface}-controls.h5ad"
  if [[ -f "$controls" && -f "${controls}.manifest.json" ]]; then
    echo "reusing control-only artifact: $controls"
  elif [[ -e "$controls" || -e "${controls}.manifest.json" ]]; then
    echo "incomplete control-only artifact: $controls" >&2
    exit 2
  else
    uv run python "$GUARD" \
      --max-rss-gib 24 \
      --min-free-percent 15 \
      --poll-seconds 5 \
      -- \
      "$ROOT/.venv/bin/python" "$HERE/prepare_held_controls.py" \
        --line "$line" \
        --output "$controls"
  fi
done

assert_no_other_heavy
uv run python "$GUARD" \
  --max-rss-gib 24 \
  --min-free-percent 15 \
  --poll-seconds 5 \
  -- \
  "$ROOT/.venv/bin/python" "$HERE/prepare_held_pruned_operators.py" \
    --max-rss-gib 20

for surface in hct116 hek293t; do
  if [[ -d "$INCUMBENT_DIR/$surface" ]]; then
    echo "reusing held incumbent artifact: $INCUMBENT_DIR/$surface"
  else
    uv run python "$GUARD" \
      --max-rss-gib 24 \
      --min-free-percent 15 \
      --poll-seconds 5 \
      -- \
      uv run python "$HERE/freeze_a14_endpoints.py" \
        --surface "$surface" \
        --controls-h5ad "$CONTROL_DIR/orion-${surface}-controls.h5ad" \
        --output-dir "$INCUMBENT_DIR/$surface"
  fi
done

assert_no_other_heavy
if [[ -d "$WORK/hct116" ]]; then
  echo "reusing STATE inference artifact: $WORK/hct116"
else
  STATE_USE_MPS=1 PYTORCH_ENABLE_MPS_FALLBACK=0 \
    uv run python "$GUARD" \
      --max-rss-gib 24 \
      --min-free-percent 15 \
      --poll-seconds 5 \
      -- \
      "$STATE_PY" "$HERE/state_effect_adapter.py" \
        --controls-h5ad "$CONTROL_DIR/orion-hct116-controls.h5ad" \
        --incumbent-npz "$INCUMBENT_DIR/hct116/frozen.npz" \
        --incumbent-manifest "$INCUMBENT_DIR/hct116/manifest.json" \
        --context HCT116 \
        --checkpoint-manifest "$FROZEN_CHECKPOINT/checkpoint-manifest.json" \
        --run-dir "$FROZEN_CHECKPOINT" \
        --output-dir "$WORK/hct116" \
        --evaluation-run-id "$RUN_ID"
fi

assert_no_other_heavy
if [[ -d "$WORK/hek293t" ]]; then
  echo "reusing STATE inference artifact: $WORK/hek293t"
else
  STATE_USE_MPS=1 PYTORCH_ENABLE_MPS_FALLBACK=0 \
    uv run python "$GUARD" \
      --max-rss-gib 24 \
      --min-free-percent 15 \
      --poll-seconds 5 \
      -- \
      "$STATE_PY" "$HERE/state_effect_adapter.py" \
        --controls-h5ad "$CONTROL_DIR/orion-hek293t-controls.h5ad" \
        --incumbent-npz "$INCUMBENT_DIR/hek293t/frozen.npz" \
        --incumbent-manifest "$INCUMBENT_DIR/hek293t/manifest.json" \
        --context HEK293T \
        --checkpoint-manifest "$FROZEN_CHECKPOINT/checkpoint-manifest.json" \
        --run-dir "$FROZEN_CHECKPOINT" \
        --output-dir "$WORK/hek293t" \
        --evaluation-run-id "$RUN_ID"
fi

if [[ ! -f "$WORK/common-targets.txt" ]]; then
  "$ROOT/.venv/bin/python" "$HERE/score_state_blends.py" common-targets \
    --frozen-dir "$WORK/hct116" \
    --frozen-dir "$WORK/hek293t" \
    --output "$WORK/common-targets.txt"
fi

for surface in hct116 hek293t; do
  if [[ ! -f "$WORK/$surface/blend-metrics.json" ]]; then
    "$ROOT/.venv/bin/python" "$HERE/score_state_blends.py" score \
      --surface "$surface" \
      --frozen-dir "$WORK/$surface" \
      --target-allowlist "$WORK/common-targets.txt"
  fi
done

if [[ ! -f "$WORK/summary.json" ]]; then
  "$ROOT/.venv/bin/python" "$HERE/score_state_blends.py" summarize \
    --surface-dir "$WORK/hct116" \
    --surface-dir "$WORK/hek293t" \
    --output "$WORK/summary.json"
fi

echo "experiment-106 summary: $WORK/summary.json"
