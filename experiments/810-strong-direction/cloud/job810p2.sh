#!/usr/bin/env bash
# 810 -- job810p.sh (ea9e743f...) with ONE change: the fold's part rows come from the input mirror (synced to
# experiments/810-strong-direction/work/parts/ by labworker; a server-side S3 copy of output/lab/810parts/, same keys and sizes)
# because the instance role may not read s3://.../output/lab/ (run 20261016T081954Z-s810-price-a: 403 on HeadObject, 22/22
# jobs failed before any pricing). Pricing code unchanged (price810.py da001560...).
set -Eeuo pipefail
cd /work/repo
fid="$1"; stem="${fid/\//__}"
for d in experiments/*/; do mkdir -p "${d}work"; done
export GUARD801_ALLOW_MKDIR=1 KA804_OUT=/work/ka804 P810_OUT=/work/p810 PARTS810=/work/repo/experiments/810-strong-direction/work/parts
export OPENBLAS_NUM_THREADS=2 OMP_NUM_THREADS=2 MKL_NUM_THREADS=2
uv run python experiments/810-strong-direction/price810.py "$fid"
D=experiments/810-strong-direction/work/price
mkdir -p "$D/bulk/$stem" "$D/meta"
cp /work/p810/bulk/$stem/*.npz "$D/bulk/$stem/"
cp /work/p810/meta/$stem.json "$D/meta/"
echo "$fid priced"
