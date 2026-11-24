#!/usr/bin/env bash
# 318 prologue — build `base` + `ptshuf2` ON the worker, run every known-answer gate, keep only
# the arm for scoring. `base` is NOT scored: it is 221's rc_callset (gate 5 asserts the
# sha256) and 317 already recorded its six-member row on a fresh worker from the same object.
set -Eeuo pipefail
export PATH="/root/.local/bin:$PATH"
export UV_CACHE_DIR="${WORK}/uv-cache"
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-8}"
export OPENBLAS_NUM_THREADS="${OPENBLAS_NUM_THREADS:-8}"
export MKL_NUM_THREADS="${MKL_NUM_THREADS:-8}"
mkdir -p "${WORK}/prologue-out" "${PREDS_DIR}"

# the 1.29 GB operator and the three staged npz, uploaded once for 317's sweep
aws s3 sync s3://<your-s3-bucket>/input/cloud-score/ladder5-20260904/payload/ "$(pwd)/" --only-show-errors
ls -la data/donors/ experiments/170-gamma-mse/work/ 2>&1 | head -20

uv sync --frozen

VCC_ALPHA=0.982 VCC_MASS_FIX=minnorm VCC_WHITEN_GAMMA=0.65 VCC_WHITEN_LAM=0.0 \
VCC_CD4_TARGET_F=0.10 \
uv run python experiments/318-per-target-magnitude/build318.py \
  --real "$REAL" --out-dir "$PREDS_DIR" \
  --report "${WORK}/prologue-out/build318-ptshuf2.json" \
  --only "base,ptshuf2"

sha256sum "${PREDS_DIR}"/*.h5ad >"${WORK}/prologue-out/candidate-sha256.txt"
cat "${PREDS_DIR}"/../prologue-out/candidate-sha256.txt 2>/dev/null || cat "${WORK}/prologue-out/candidate-sha256.txt"
rm -f "${PREDS_DIR}/base.h5ad"
