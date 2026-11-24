#!/bin/bash
set -Eeuo pipefail

BUCKET="<your-s3-bucket>"
PREFIX="exp057-kolf-union-whiten/20260822T210455Z"
ROOT="/work/repo/experiments/057-kolf-union-whiten"
INPUTS="${ROOT}/inputs"
OUT="${ROOT}/out"
LOG="${ROOT}/cloud.log"
PYTHON="/work/cd4/repo/.venv/bin/python"

mkdir -p "${INPUTS}" "${OUT}"
exec > >(tee -a "${LOG}") 2>&1

finish() {
  status=$?
  trap - EXIT
  printf '%s\n' "${status}" > "${ROOT}/exit-code.txt"
  aws s3 cp "${LOG}" "s3://${BUCKET}/${PREFIX}/output/cloud.log" \
    --only-show-errors || true
  aws s3 cp "${ROOT}/exit-code.txt" \
    "s3://${BUCKET}/${PREFIX}/output/exit-code.txt" \
    --only-show-errors || true
  aws s3 cp "${OUT}/" "s3://${BUCKET}/${PREFIX}/output/out/" \
    --recursive --only-show-errors || true
  exit "${status}"
}
trap finish EXIT

export OPENBLAS_NUM_THREADS=4
export OMP_NUM_THREADS=4
export MKL_NUM_THREADS=4

aws s3 cp "s3://${BUCKET}/${PREFIX}/input/run.py" "${ROOT}/run.py" \
  --only-show-errors
aws s3 cp "s3://${BUCKET}/${PREFIX}/input/gene_names.csv" \
  "${INPUTS}/gene_names.csv" --only-show-errors
aws s3 cp "s3://${BUCKET}/${PREFIX}/input/pert_counts.csv" \
  "${INPUTS}/pert_counts.csv" --only-show-errors
aws s3 cp "s3://${BUCKET}/${PREFIX}/input/exp055-results.json" \
  "${INPUTS}/exp055-results.json" --only-show-errors
aws s3 cp "s3://${BUCKET}/${PREFIX}/input/kolf2/" "${INPUTS}/kolf2/" \
  --recursive --only-show-errors

for donor in orion_hct116 orion_hek293t k562_gwps
do
  mkdir -p "${INPUTS}/${donor}"
  aws s3 cp "s3://${BUCKET}/input/clean-hct/${donor}/" \
    "${INPUTS}/${donor}/" --recursive --only-show-errors
done

cd "${ROOT}"
/usr/bin/time -v "${PYTHON}" run.py \
  --output-dir "${OUT}" \
  --basis-per-donor 6000 \
  --operator-seed 0 \
  --gammas 0 0.55 0.65 0.75 0.85 \
  --query-folds 5 \
  --chunk-rows 64 \
  --max-rss-gib 28.0

"${PYTHON}" - <<'PY'
import hashlib
import json
from pathlib import Path

root = Path("/work/repo/experiments/057-kolf-union-whiten")
files = {}
for path in sorted((root / "out").glob("*")):
    if not path.is_file():
        continue
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while block := handle.read(8 << 20):
            digest.update(block)
    files[path.name] = {"bytes": path.stat().st_size, "sha256": digest.hexdigest()}
(root / "artifact-manifest.json").write_text(json.dumps(files, indent=2) + "\n")
PY

aws s3 cp "${ROOT}/artifact-manifest.json" \
  "s3://${BUCKET}/${PREFIX}/output/artifact-manifest.json" \
  --only-show-errors
