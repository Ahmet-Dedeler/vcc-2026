#!/bin/bash
# 424: price one VIPerturb table variant end to end (4 PLOLO folds, A/B/C energies, pricing), serially.
set -u
cd "$(dirname "$0")"
v=$1
for f in orion_hct116 orion_hek293t kolf2 k562_gwps; do
  VCC424_VIP=$v /usr/bin/time -l uv run python comp424.py $f > work/comp424_$v-$f.log 2>&1 || { echo "comp $f failed"; exit 1; }
done
VCC424_VIP=$v /usr/bin/time -l uv run python abc424.py > work/abc424_$v.log 2>&1 || { echo "abc failed"; exit 1; }
VCC424_VIP=$v uv run python analyse424.py > work/analyse424_$v.stdout 2>&1 || { echo "analyse failed"; exit 1; }
echo "done $v"
