#!/bin/bash
# 804: row 710 on the three ABCp folds (EC2 only; labworker runs it in /work/repo).
set -uo pipefail
for d in experiments/*/; do mkdir -p "${d}work"; done
export GUARD801_ALLOW_MKDIR=1 OPENBLAS_NUM_THREADS=4 OMP_NUM_THREADS=4 MKL_NUM_THREADS=4 KA804_OUT=/work/ka804p
printf 'ABCp/k562_gwps\nABCp/orion_hct116\nABCp/orion_hek293t\n' | xargs -P 3 -I{} sh -c 'uv run python experiments/804-known-answers/price710p.py {} > /work/p710p-$(echo {} | tr / _).log 2>&1; echo "{} rc=$?"'
cat /work/p710p-*.log | grep -E "done|Error|Traceback|KA|cos_agg" | head -40
mkdir -p experiments/804-known-answers/work/run710p && cp -R /work/ka804p/. experiments/804-known-answers/work/run710p/
n=$(ls experiments/804-known-answers/work/run710p/bulk/*/710.npz 2>/dev/null | wc -l); echo "710 arrays: $n"
[ "$n" -eq 3 ]
