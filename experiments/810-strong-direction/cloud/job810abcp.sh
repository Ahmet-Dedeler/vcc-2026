#!/usr/bin/env bash
# 810 -- the ABCp information read for one A/B/C fold (EC2 only; LEAD-DECISION-1.md): parts with 062's operator (KA vs 518's dump
# and the GO run's parts), parts with the DAY operator (KA vs 801's ABCp dump), then price810.py (frozen 9c33b24f...) on ABCp/<held>.
set -Eeuo pipefail
cd /work/repo
h="$1"
for d in experiments/*/; do mkdir -p "${d}work"; done
export GUARD801_ALLOW_MKDIR=1 KA804_OUT=/work/ka804 P810_OUT=/work/p810 PARTS810=/work/parts810 PLOLO810=/work/plolo810
export OPENBLAS_NUM_THREADS=2 OMP_NUM_THREADS=2 MKL_NUM_THREADS=2
D=experiments/810-strong-direction/work/abcp
mkdir -p "$D/parts" "$D/bulk" "$D/meta"
uv run python experiments/810-strong-direction/decompabcp810.py "$h" 062
uv run python experiments/810-strong-direction/decompabcp810.py "$h" day
cp /work/parts810/ABC__$h-parts.json /work/parts810/ABCp__$h-parts.json "$D/parts/"
uv run python experiments/810-strong-direction/price810.py "ABCp/$h"
mkdir -p "$D/bulk/ABCp__$h"
cp /work/p810/bulk/ABCp__$h/*.npz "$D/bulk/ABCp__$h/"
cp /work/p810/meta/ABCp__$h.json "$D/meta/"
echo "ABCp/$h priced"
