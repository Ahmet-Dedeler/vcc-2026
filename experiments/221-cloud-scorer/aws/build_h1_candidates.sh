#!/usr/bin/env bash
# 221 prologue — build the four board-anchored h1-test candidates ON the worker.
#
# `scripts/cloud_score_worker.sh` runs this in $REPO with WORK, REAL and PREDS_DIR exported,
# BEFORE it scores. Everything heavy (the 1.29 GB operator, the 2.7 GB control block, four
# 52,400 x 18,080 emissions) therefore happens in the cloud, which is the whole point of the
# mandate — the laptop never goes above the 8 GB ceiling.
set -Eeuo pipefail

export PATH="/root/.local/bin:$PATH"
export UV_CACHE_DIR="${WORK}/uv-cache"
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-8}"
export OPENBLAS_NUM_THREADS="${OPENBLAS_NUM_THREADS:-8}"
export MKL_NUM_THREADS="${MKL_NUM_THREADS:-8}"
mkdir -p "${WORK}/prologue-out" "${PREDS_DIR}"

echo "[prologue] repo: $(pwd)"
ls -la data/donors/ experiments/170-gamma-mse/work/ 2>&1 | head -30

uv sync --frozen

VCC_ALPHA=0.982 VCC_MASS_FIX=minnorm VCC_WHITEN_GAMMA=0.65 VCC_WHITEN_LAM=0.0 \
VCC_CD4_TARGET_F=0.10 \
uv run python experiments/221-cloud-scorer/build_h1_candidates.py \
  --real "$REAL" \
  --out-dir "$PREDS_DIR" \
  --report "${WORK}/prologue-out/build221.json"

sha256sum "${PREDS_DIR}"/*.h5ad >"${WORK}/prologue-out/candidate-sha256.txt"
cat "${WORK}/prologue-out/candidate-sha256.txt"
