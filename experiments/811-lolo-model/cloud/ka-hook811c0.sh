#!/bin/bash
# 811 C0 hook: KA-same-construction / KA-scale on engine folds (EC2 only): kasame811c0.py on ABCp x3 + L1/k562.
set -uo pipefail
for d in experiments/*/; do mkdir -p "${d}work"; done
export GUARD801_ALLOW_MKDIR=1 OPENBLAS_NUM_THREADS=2 OMP_NUM_THREADS=2 MKL_NUM_THREADS=2 KASAME811_OUT=/work/kasame811
O=experiments/811-lolo-model/work/hook811/kac0; mkdir -p $O; rc=0
for f in ABCp/k562_gwps ABCp/orion_hct116 ABCp/orion_hek293t L1/k562_gwps; do
  uv run python experiments/811-lolo-model/kasame811c0.py $f > $O/kasame-$(echo $f | tr / _).log 2>&1 && echo "KA-same-C0 $f: PASS" | tee -a $O/summary.txt || { echo "KA-same-C0 $f: FAIL" | tee -a $O/summary.txt; rc=1; }
done
cp /work/kasame811/c0-*.json $O/ 2>/dev/null
exit $rc
