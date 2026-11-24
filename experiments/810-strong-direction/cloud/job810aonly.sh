#!/usr/bin/env bash
# 810 -- the A-only instrument read (EC2 only; lab lead 10-16 12:44Z). Two subcommands, run by labworker from /work/repo:
#   prep <A_k>          the panel's inputs in the EC2 checkout, the way 801's `cloud801.sh prep` makes them (199's def-bundle
#                       controls, sha256-checked; cd4_801.py's CD4 extraction; frames801.py under the shared frame lock), but
#                       without 801's S3 probe; operators from 810's server-side copy of output/lab/ops801/<A_k>
#                       (input/lab/repo/experiments/810-strong-direction/work/aops/<A_k>), each checked against 801's registry.
#   fold <A_k> <held>   after its panel's prep: decomp810 compose (KA: parts sum == 801's dump708 fin < 1e-12, h1c1 == fPP),
#                       then price810.py (frozen 9c33b24f...) on <A_k>/<held>; arrays to work/aonly/{bulk,meta,parts}.
set -Eeuo pipefail
cd /work/repo
for d in experiments/*/; do mkdir -p "${d}work"; done
export GUARD801_ALLOW_MKDIR=1 KA804_OUT=/work/ka804 OPENBLAS_NUM_THREADS=2 OMP_NUM_THREADS=2 MKL_NUM_THREADS=2
W=experiments/801-lab-folds/work
C=/work/aonly-markers
D=experiments/810-strong-direction/work/aonly
mkdir -p $C $D/bulk $D/meta $D/parts
sub="$1"; L="$2"
[[ "$L" =~ ^A[1-6]$ ]] || { echo "only A1-A6"; exit 2; }
if [[ "$sub" == prep ]]; then
  trap 'touch '"$C/$L"'.fail' ERR
  B=$W/$L/bundle
  for x in D E F; do cp experiments/199-def-runbook/work/def-bundle/context_$x.h5ad $B/; done
  ( cd $B && sha256sum context_D.h5ad context_E.h5ad context_F.h5ad ) | while read -r s f; do
    grep -q "^$s  experiments/199-def-runbook/work/def-bundle/$f\$" $W/manifest/local_sha256_defbundle.txt || { echo "$L: $f differs"; exit 7; }
  done
  uv run python experiments/801-lab-folds/cd4_801.py "$L"
  mkdir -p $W/$L/ops
  cp experiments/810-strong-direction/work/aops/$L/held_*.npz experiments/810-strong-direction/work/aops/$L/held_*.json $W/$L/ops/
  for h in k562_gwps orion_hct116 orion_hek293t; do
    got="$(sha256sum $W/$L/ops/held_$h.npz | cut -d' ' -f1)"; want="$(jq -r ".[\"$L/$h\"].operator_sha256" experiments/810-strong-direction/work/aops/registry-aonly.json)"
    [[ "$got" == "$want" ]] || { echo "$L/$h operator sha $got != registry $want"; exit 11; }
  done
  flock /work/frames.lock uv run python experiments/801-lab-folds/frames801.py "$L"
  touch $C/$L.ok
  echo "$L prepped"
else
  h="$3"; i=0
  while [[ ! -f $C/$L.ok ]]; do [[ -f $C/$L.fail ]] && { echo "prep $L failed"; exit 8; }; sleep 10; i=$((i+1)); [[ $i -lt 240 ]] || exit 4; done
  stem="${L}__$h"
  export OUT810=/work/out810-$stem PARTS810=/work/parts810-$stem P810_OUT=/work/p810-$stem
  uv run python experiments/810-strong-direction/decomp810.py compose "$L/$h"
  cp $OUT810/$stem-comp.json $D/parts/
  uv run python experiments/810-strong-direction/price810.py "$L/$h"
  mkdir -p $D/bulk/$stem
  cp $P810_OUT/bulk/$stem/*.npz $D/bulk/$stem/
  cp $P810_OUT/meta/$stem.json $D/meta/
  rm -rf $PARTS810
  echo "$L/$h priced"
fi
