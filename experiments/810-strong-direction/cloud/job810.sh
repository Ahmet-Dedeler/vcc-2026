#!/usr/bin/env bash
# 810 -- one fold's decomposition on EC2 (run by labworker from /work/repo; never on the Mac).
#   job810.sh <panel>/<held>   L folds wait for 801's prep marker of their panel; ABC folds compose under 801's frame lock.
#   compose (parts, KA vs the registry dump) then realise (engine801 incumbent; ABC: KA-706); outputs copied into
#   experiments/810-strong-direction/work/decomp/ (collected) and the full part rows to s3://$BUCKET/output/lab/810parts/.
set -Eeuo pipefail
cd /work/repo
fid="$1"; panel="${fid%%/*}"; stem="${fid/\//__}"
C=experiments/801-lab-folds/work/cloud
D=experiments/810-strong-direction/work/decomp
mkdir -p "$C" "$D"
for d in experiments/*/; do mkdir -p "${d}work"; done
export GUARD801_ALLOW_MKDIR=1 OUT810=/work/out810 PARTS810=/work/parts810
export OPENBLAS_NUM_THREADS=2 OMP_NUM_THREADS=2 MKL_NUM_THREADS=2
if [[ "$panel" == ABC ]]; then
  flock $C/frames.lock uv run python experiments/810-strong-direction/decomp810.py compose "$fid"
else
  i=0
  while [[ ! -f $C/$panel.prep.ok ]]; do [[ -f $C/$panel.prep.fail ]] && { echo "prep $panel failed"; exit 8; }; sleep 10; i=$((i+1)); [[ $i -lt 240 ]] || exit 4; done
  uv run python experiments/810-strong-direction/decomp810.py compose "$fid"
fi
uv run python experiments/810-strong-direction/decomp810.py realise "$fid"
cp /work/out810/${stem}-* "$D/"
aws s3 cp /work/parts810/${stem}-parts.npz "s3://${BUCKET:-<your-s3-bucket>}/output/lab/810parts/${stem}-parts.npz" --only-show-errors
rm -f /work/parts810/${stem}-parts.npz
echo "$fid done"
