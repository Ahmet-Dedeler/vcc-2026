#!/usr/bin/env bash
# 810 -- price810 on one fold (EC2 only; run by labworker from /work/repo; ONLY after the judge's GO on PROPOSAL-2.md).
#   job810p.sh <panel>/<held>: fetch the fold's part rows from s3://$BUCKET/output/lab/810parts/, price inc / 707 / sgu0 (KA) /
#   sgu1 / sgu1all (price810.py), copy the STANDARD-FORMAT arrays into experiments/810-strong-direction/work/price/ (collected).
#   Needs sync810.txt plus experiments/804-known-answers/work/run (804's inc / 707 arrays: the KA).
set -Eeuo pipefail
cd /work/repo
fid="$1"; stem="${fid/\//__}"
for d in experiments/*/; do mkdir -p "${d}work"; done
export GUARD801_ALLOW_MKDIR=1 KA804_OUT=/work/ka804 P810_OUT=/work/p810 PARTS810=/work/parts810
export OPENBLAS_NUM_THREADS=2 OMP_NUM_THREADS=2 MKL_NUM_THREADS=2
mkdir -p /work/parts810
aws s3 cp "s3://${BUCKET:-<your-s3-bucket>}/output/lab/810parts/${stem}-parts.npz" /work/parts810/ --only-show-errors
uv run python experiments/810-strong-direction/price810.py "$fid"
D=experiments/810-strong-direction/work/price
mkdir -p "$D/bulk/$stem" "$D/meta"
cp /work/p810/bulk/$stem/*.npz "$D/bulk/$stem/"
cp /work/p810/meta/$stem.json "$D/meta/"
rm -f /work/parts810/${stem}-parts.npz
echo "$fid priced"
