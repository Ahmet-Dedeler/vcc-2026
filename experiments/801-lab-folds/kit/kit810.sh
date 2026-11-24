#!/usr/bin/env bash
# kit810.sh -- 810's A-only job (experiments/810-strong-direction/cloud/job810aonly.sh, git blob c9586229616791983afd4e04548b708bfb006af2)
# with the panel TAG as a parameter (801 follow-up 7: the panel kit's output feeds it). Changes, nothing else: any kit tag
# (not only A1-A6); operators from 801's mirror of the kit's fits (experiments/801-lab-folds/work/<TAG>/ops, copied there by
# kit801-collect.sh), each checked against the kit's registry fragment (work/registry-<TAG>.json); outputs under
# experiments/810-strong-direction/work/kit-<TAG>/. price810.py / decomp810.py are called unchanged.
#   prep <TAG>          |   fold <TAG> <held>
set -Eeuo pipefail
cd /work/repo
for d in experiments/*/; do mkdir -p "${d}work"; done
export GUARD801_ALLOW_MKDIR=1 KA804_OUT=/work/ka804 OPENBLAS_NUM_THREADS=2 OMP_NUM_THREADS=2 MKL_NUM_THREADS=2
W=experiments/801-lab-folds/work
C=/work/kit810-markers
D=experiments/810-strong-direction/work/kit-$2
mkdir -p $C $D/bulk $D/meta $D/parts
sub="$1"; L="$2"
[[ "$L" =~ ^[A-Za-z][A-Za-z0-9]{0,15}$ && "$L" != L7 && "$L" != L8 ]] || { echo "bad tag $L"; exit 2; }
if [[ "$sub" == prep ]]; then
  trap 'touch '"$C/$L"'.fail' ERR
  B=$W/$L/bundle
  for x in D E F; do cp experiments/199-def-runbook/work/def-bundle/context_$x.h5ad $B/; done
  ( cd $B && sha256sum context_D.h5ad context_E.h5ad context_F.h5ad ) | while read -r s f; do
    grep -q "^$s  experiments/199-def-runbook/work/def-bundle/$f\$" $W/manifest/local_sha256_defbundle.txt || { echo "$L: $f differs"; exit 7; }
  done
  uv run python experiments/801-lab-folds/cd4_801.py "$L"
  mkdir -p $W/$L/ops
  ls $W/$L/ops/held_*.npz >/dev/null   # the kit mirror already placed them here
  for h in k562_gwps orion_hct116 orion_hek293t; do
    got="$(sha256sum $W/$L/ops/held_$h.npz | cut -d' ' -f1)"; want="$(jq -r ".[\"$L/$h\"].operator_sha256" $W/registry-$L.json)"
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
