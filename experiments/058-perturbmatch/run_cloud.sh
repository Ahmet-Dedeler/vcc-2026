#!/usr/bin/env bash
set -euo pipefail

BUCKET="${BUCKET:-<your-s3-bucket>}"
WORK="${WORK:-/work/exp058}"
PYTHON="${PYTHON:-/work/cd4/repo/.venv/bin/python}"
INPUT="s3://${BUCKET}/input/exp058"
OUTPUT="s3://${BUCKET}/output/exp058"

mkdir -p "${WORK}/code" "${WORK}/results"
for file in extract.py evaluate_lodo.py test_perturbmatch.py pert_counts.csv gene_names.csv
do
  aws s3 cp "${INPUT}/${file}" "${WORK}/code/${file}" --only-show-errors
done

cd "${WORK}/code"

if [[ ! -f "${WORK}/perturbmatch/meta.json" ]]
then
  /usr/bin/time -v "${PYTHON}" extract.py \
    --panel-csv "${WORK}/code/pert_counts.csv" \
    --output "${WORK}/perturbmatch" \
    2>&1 | tee "${WORK}/results/extract.log"
fi

tar -C "${WORK}" -czf "${WORK}/results/perturbmatch-compact.tar.gz" perturbmatch
sha256sum "${WORK}/results/perturbmatch-compact.tar.gz" \
  > "${WORK}/results/perturbmatch-compact.tar.gz.sha256"

if [[ ! -f "${WORK}/results/results.json" ]]
then
  /usr/bin/time -v "${PYTHON}" evaluate_lodo.py \
    --k562 /work/cd4/donors/k562_gwps \
    --hct116 /work/cd4/donors/orion_hct116 \
    --hek293t /work/cd4/donors/orion_hek293t \
    --perturbmatch "${WORK}/perturbmatch" \
    --panel-csv "${WORK}/code/pert_counts.csv" \
    --genes-csv "${WORK}/code/gene_names.csv" \
    --max-rss-gib 12 \
    --output "${WORK}/results/results.json" \
    2>&1 | tee "${WORK}/results/evaluate.log"
fi

for file in \
  extract.log \
  evaluate.log \
  results.json \
  results.csv \
  per_target_pds.csv \
  perturbmatch-compact.tar.gz \
  perturbmatch-compact.tar.gz.sha256
do
  aws s3 cp "${WORK}/results/${file}" "${OUTPUT}/${file}" --only-show-errors
done
