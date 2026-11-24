#!/bin/bash
# 229 — the FOURTH build: `stack6`, the buildable form of "the new sources ride ALONGSIDE the
# H1 pair" (NEXT-BUILD-SPEC §4).  Waits for the three-build chain to clear, because two builds
# must never overlap on this machine.
#
# Measured, three cancer PLOLO folds, marginal over the incumbent `h1union`, at g = 0.35:
#   ARM6  +0.00268 (1.41x the shape-only gate), 3/3 cosine, 3/3 pds, F1 and F3 pass
#   ARM4  -0.00102 (the literally pre-registered source list, which DROPS the H1 pair)
set -uo pipefail
cd "$(dirname "$0")/../.."
HERE=experiments/229-stack4-build

while pgrep -f "make_submission_pm_cd4.py --whitener" > /dev/null; do sleep 60; done
echo "=== stack6g035 $(date -u +%FT%TZ) ==="
env VCC_H1_UNION_G=0.35 \
    VCC_UNION_SOURCES=gse295214_h1,h1_2025,dld1_lowmoi,jiang_5line_noK562,igvf_cardiomyocyte,myllia_a549 \
    "$HERE/build229.sh" submissions/229-stack6-g035.h5ad > "$HERE/work/build-stack6g035.log" 2>&1
echo "stack6g035 exit=$?"
if [ -f submissions/229-stack6-g035.h5ad ]; then
  echo -n "stack6g035 sha256 "; shasum -a 256 submissions/229-stack6-g035.h5ad | awk '{print $1}'
  echo -n "stack6g035 bytes  "; stat -f %z submissions/229-stack6-g035.h5ad
fi
grep -E "maximum resident set size" "$HERE/work/build-stack6g035.log" | tail -1
echo "=== DONE $(date -u +%FT%TZ) ==="
