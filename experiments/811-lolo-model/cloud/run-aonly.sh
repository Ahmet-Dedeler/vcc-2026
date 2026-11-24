#!/bin/bash
# 811 A-only instrument read (EC2 only; labworker runs it in /work/repo): extend the banks (ext811a), then price 18 A folds
# + the L1/k562_gwps known-answer fold with aonly811.py, 10 at a time.
set -uo pipefail
for d in experiments/*/; do mkdir -p "${d}work"; done
export GUARD801_ALLOW_MKDIR=1 OPENBLAS_NUM_THREADS=2 OMP_NUM_THREADS=2 MKL_NUM_THREADS=2 AONLY811_OUT=/work/aonly811
uv run python experiments/811-lolo-model/ext811a.py || { echo "ext811a FAILED"; exit 3; }
FOLDS="L1/k562_gwps $(for k in 1 2 3 4 5 6; do for h in k562_gwps orion_hct116 orion_hek293t; do printf 'A%s/%s ' $k $h; done; done)"
echo $FOLDS | tr ' ' '\n' | xargs -P 10 -I{} sh -c 'uv run python experiments/811-lolo-model/aonly811.py {} > /work/a811-$(echo {} | tr / _).log 2>&1; echo "{} rc=$?"' | tee /work/a811-rc.txt
mkdir -p experiments/811-lolo-model/work/run-m811aonly/logs && cp -R /work/aonly811/. experiments/811-lolo-model/work/run-m811aonly/ && cp /work/a811-*.log /work/a811-rc.txt experiments/811-lolo-model/work/run-m811aonly/logs/
grep -h "KA \|_l25 " /work/a811-*.log | head -100
! grep -v "rc=0" /work/a811-rc.txt
