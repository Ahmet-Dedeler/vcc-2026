#!/usr/bin/env bash
# 810 -- the LABC instrument read (EC2 only; lab lead 10-16 ~10:55Z). One job for the three LABC folds:
#   inputs arranged in the EC2 checkout (never on the Mac): decomp810's L-panel path reads work/<panel>/{bundle,cd4/cd4_compact,ops};
#   for LABC those are data/challenge's bundle files, A/B/C's production CD4 banks (byte-identical to 801's LABC banks, 801 §14)
#   and ABCp's day operators (the operators 801's LABC dumps used). decomp810's KA (parts sum == the LABC dump's fin < 1e-12) is
#   what proves the arrangement right. Then price810 (frozen 9c33b24f...), labc810.py's KA gate, check800 aggregate --no-wear x3.
set -Eeuo pipefail
cd /work/repo
for d in experiments/*/; do mkdir -p "${d}work"; done
export GUARD801_ALLOW_MKDIR=1 KA804_OUT=/work/ka804 P810_OUT=/work/p810 PARTS810=/work/parts810 OUT810=/work/out810
export OPENBLAS_NUM_THREADS=3 OMP_NUM_THREADS=3 MKL_NUM_THREADS=3
W=experiments/801-lab-folds/work/LABC
mkdir -p $W/bundle $W/cd4 $W/ops
cp data/challenge/context_A.h5ad data/challenge/context_B.h5ad data/challenge/context_C.h5ad data/challenge/manifest.json $W/bundle/
cmp data/challenge/pert_counts.csv $W/bundle/pert_counts.csv && cmp data/challenge/gene_names.csv $W/bundle/gene_names.csv
cp -R data/external/cd4t/cloud-compact/cd4_compact $W/cd4/cd4_compact
cp experiments/801-lab-folds/work/ABCp/ops/held_*.npz experiments/801-lab-folds/work/ABCp/ops/held_*.json $W/ops/
D=experiments/810-strong-direction/work/labc
mkdir -p $D/bulk $D/meta $D/parts
for h in k562_gwps orion_hct116 orion_hek293t; do
  uv run python experiments/810-strong-direction/decomp810.py compose "LABC/$h"
  cp /work/out810/LABC__$h-comp.json $D/parts/
  uv run python experiments/810-strong-direction/price810.py "LABC/$h"
  mkdir -p $D/bulk/LABC__$h
  cp /work/p810/bulk/LABC__$h/*.npz $D/bulk/LABC__$h/
  cp /work/p810/meta/LABC__$h.json $D/meta/
done
uv run python experiments/810-strong-direction/labc810.py $D
for a in sgu1 sgu1all 707; do
  uv run python experiments/800-checker/check800.py aggregate $D --cand $a --base inc --family 810-labc-read --cls transfer --no-wear > $D/agg-$a.txt 2>&1
done
echo "LABC read done"
