#!/bin/bash
# run-kit811.sh <TAG> -- 811's run-aonly.sh (git blob 33f5558cd1648495eebcd6f3d0d990ee47af1e7a) with the panel TAG as a parameter
# (801 follow-up 7): extend 811's banks to <TAG>'s targets (kit811ext.py = ext811a with KIT811_PANELS), then price <TAG>'s 3 cancer
# folds + the L1/k562_gwps known-answer fold with kit811.py (= aonly811 with the fold assert opened to the kit tag). Outputs
# experiments/811-lolo-model/work/run-m811kit-<TAG>/. For the lead to launch through 811 on the day.
set -uo pipefail
TAG="${1:?TAG}"
for d in experiments/*/; do mkdir -p "${d}work"; done
export GUARD801_ALLOW_MKDIR=1 OPENBLAS_NUM_THREADS=2 OMP_NUM_THREADS=2 MKL_NUM_THREADS=2 AONLY811_OUT=/work/kit811 KIT811_PANELS="$TAG"
uv run python experiments/801-lab-folds/kit/kit811ext.py || { echo "kit811ext FAILED"; exit 3; }
FOLDS="L1/k562_gwps $TAG/k562_gwps $TAG/orion_hct116 $TAG/orion_hek293t"
echo $FOLDS | tr ' ' '\n' | xargs -P 4 -I{} sh -c 'uv run python experiments/801-lab-folds/kit/kit811.py {} > /work/k811-$(echo {} | tr / _).log 2>&1; echo "{} rc=$?"' | tee /work/k811-rc.txt
O=experiments/811-lolo-model/work/run-m811kit-$TAG
mkdir -p $O/logs && cp -R /work/kit811/. $O/ && cp /work/k811-*.log /work/k811-rc.txt $O/logs/
grep -h "KA \|_l25 " /work/k811-*.log | head -100
! grep -v "rc=0" /work/k811-rc.txt
