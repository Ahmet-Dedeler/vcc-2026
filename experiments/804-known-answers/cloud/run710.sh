#!/bin/bash
# 804: the 710 row on all 22 bulk folds (11 at a time), then grade every row present (EC2 only; labworker runs it in /work/repo).
set -uo pipefail
for d in experiments/*/; do mkdir -p "${d}work"; done
export GUARD801_ALLOW_MKDIR=1 OPENBLAS_NUM_THREADS=2 OMP_NUM_THREADS=2 MKL_NUM_THREADS=2 KA804_OUT=/work/ka804
FOLDS="$(for k in 1 2 3 4 5 6; do for h in k562_gwps orion_hct116 orion_hek293t; do echo L$k/$h; done; done; for h in k562_gwps orion_hct116 orion_hek293t kolf2; do echo ABC/$h; done)"
echo "$FOLDS" | xargs -P 11 -I{} sh -c 'uv run python experiments/804-known-answers/price710.py {} > /work/p710-$(echo {} | tr / _).log 2>&1; echo "{} rc=$?"'
cat /work/p710-*.log | grep -E "done|Error|Traceback|KA" | head -60
mkdir -p experiments/804-known-answers/work/run && cp -R /work/ka804/. experiments/804-known-answers/work/run/
n=$(ls experiments/804-known-answers/work/run/bulk/*/710.npz | wc -l); echo "710 arrays: $n"
uv run python experiments/804-known-answers/grade804.py --surfaces W,ABC > experiments/804-known-answers/work/grade804.txt 2>&1
cat experiments/804-known-answers/work/grade804.txt
[ "$n" -eq 22 ]
