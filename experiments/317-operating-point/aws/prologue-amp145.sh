#!/usr/bin/env bash
# 317 prologue — build `base` + `amp145` ON the worker, run every known-answer gate, keep
# only the arm for scoring (`base` is already scored: 221's rc_callset, KA-1).
set -Eeuo pipefail
export PATH="/root/.local/bin:$PATH"
export UV_CACHE_DIR="${WORK}/uv-cache"
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-8}"
export OPENBLAS_NUM_THREADS="${OPENBLAS_NUM_THREADS:-8}"
export MKL_NUM_THREADS="${MKL_NUM_THREADS:-8}"
mkdir -p "${WORK}/prologue-out" "${PREDS_DIR}"

# The 1.29 GB operator and the three staged npz land straight from S3 (uploaded once for the
# whole sweep) instead of through eight per-run payload uploads.
aws s3 sync s3://<your-s3-bucket>/input/cloud-score/ladder5-20260904/payload/ "$(pwd)/" --only-show-errors
ls -la data/donors/ experiments/170-gamma-mse/work/ 2>&1 | head -20

uv sync --frozen

VCC_ALPHA=0.982 VCC_MASS_FIX=minnorm VCC_WHITEN_GAMMA=0.65 VCC_WHITEN_LAM=0.0 \
VCC_CD4_TARGET_F=0.10 \
uv run python experiments/317-operating-point/build317.py \
  --real "$REAL" --out-dir "$PREDS_DIR" \
  --report "${WORK}/prologue-out/build317-amp145.json" \
  --only "base,amp145"

sha256sum "${PREDS_DIR}"/*.h5ad >"${WORK}/prologue-out/candidate-sha256.txt"
cat "${WORK}/prologue-out/candidate-sha256.txt"
# gate 5 has already compared it; the scorer only needs the arm.
rm -f "${PREDS_DIR}/base.h5ad"
