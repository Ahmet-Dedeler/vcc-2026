#!/bin/bash
# run-kit811c.sh <TAG> -- 811's composite arms (LAB-DEF-RUNBOOK Amendment 1) on a kit panel's folds. EC2 only (labworker runs it in
# /work/repo). A copy of 801's kit/run-kit811.sh (sha256 9ccd12f114872b63…) with the pricing step swapped for
# composite811.py; 801's kit811ext.py (bank extension to <TAG>'s targets) is called UNCHANGED.
#   1. kit811ext: banks -> banks_ext for <TAG>'s targets (its own KA: L1's CD4 compact byte-identical to 801's);
#   2. ka804 references for <TAG> (if the kit produced experiments/801-lab-folds/work/ka-<tag, lower case>) -> 811's work/refcomp/<TAG>__<held>/;
#   3. composite811.py --mode day --ka on L1/k562_gwps (KA-705/707/L0/811) and ABCp/k562_gwps (KA-705/707/L0), and on <TAG>'s
#      three cancer folds (--ka if their ka804 refs exist, so KA-705/707 gate them too; else KA-L0 only via --ka is skipped);
#   4. outputs -> experiments/811-lolo-model/work/run-comp811kit-<TAG>/. Any rc != 0 fails the run.
set -uo pipefail
TAG="${1:?TAG}"
for d in experiments/*/; do mkdir -p "${d}work"; done
export GUARD801_ALLOW_MKDIR=1 OPENBLAS_NUM_THREADS=2 OMP_NUM_THREADS=2 MKL_NUM_THREADS=2 KIT811_PANELS="$TAG"
export COMPOSITE811_OUT=/work/comp811kit COMPOSITE811_BANKS=experiments/811-lolo-model/work/banks_ext
uv run python experiments/801-lab-folds/kit/kit811ext.py || { echo "kit811ext FAILED"; exit 3; }
KAREF=experiments/801-lab-folds/work/ka-$(echo "$TAG" | tr "A-Z" "a-z")/bulk     # kit801 names ka dirs by the lower-case tag
TAGKA=""
if [ -d "$KAREF" ]; then
  for h in k562_gwps orion_hct116 orion_hek293t; do
    mkdir -p experiments/811-lolo-model/work/refcomp/${TAG}__$h
    cp $KAREF/${TAG}__$h/705.npz $KAREF/${TAG}__$h/707.npz experiments/811-lolo-model/work/refcomp/${TAG}__$h/ || exit 4
  done
  TAGKA="--ka"
fi
echo "ka804 refs for $TAG: ${TAGKA:-none (KA folds L1 + ABCp gate the run)}"
{ echo "L1/k562_gwps --ka"; echo "ABCp/k562_gwps --ka"; for h in k562_gwps orion_hct116 orion_hek293t; do echo "$TAG/$h $TAGKA"; done; } \
  | xargs -P 5 -I{} sh -c 'f=$(echo {} | cut -d" " -f1); uv run python experiments/811-lolo-model/composite811.py {} --mode day > /work/c811-$(echo $f | tr / _).log 2>&1; echo "$f rc=$?"' \
  | tee /work/c811-rc.txt
O=experiments/811-lolo-model/work/run-comp811kit-$TAG
mkdir -p $O/logs && cp -R /work/comp811kit/. $O/ && cp /work/c811-*.log /work/c811-rc.txt $O/logs/
grep -h "KA \|B0 \|B1 \|C0 \|C1 " /work/c811-*.log | head -80
! grep -v "rc=0" /work/c811-rc.txt
