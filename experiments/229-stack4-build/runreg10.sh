#!/bin/bash
# 229 — INBOX-2's known-answer row: REGISTRY MEMBERSHIP CANNOT AFFECT A BUILD.
# Rebuild the recommended package with the ten-entry registry (the four October sources added,
# none of them selected).  Its .h5ad must hash to the six-entry build's
# e2b9238723fea123d1130c88dbfcf7518d5a5fa44b9a27988abfc459de1936cd, or registry membership is
# not inert and pre-registering unpriceable sources is not safe.
set -uo pipefail
cd "$(dirname "$0")/../.."
HERE=experiments/229-stack4-build
while pgrep -f "make_submission_pm_cd4.py --whitener" > /dev/null; do sleep 30; done
echo "=== reg10 $(date -u +%FT%TZ) ==="
env VCC_H1_UNION_G=0.35 \
    VCC_UNION_SOURCES=gse295214_h1,h1_2025,dld1_lowmoi,jiang_5line_noK562,igvf_cardiomyocyte,myllia_a549 \
    "$HERE/build229.sh" submissions/229-stack6-g035-reg10.h5ad > "$HERE/work/build-reg10.log" 2>&1
echo "reg10 exit=$?"
echo -n "reg10 sha256 "; shasum -a 256 submissions/229-stack6-g035-reg10.h5ad | awk '{print $1}'
echo -n "reg10 bytes  "; stat -f %z submissions/229-stack6-g035-reg10.h5ad
echo "required     e2b9238723fea123d1130c88dbfcf7518d5a5fa44b9a27988abfc459de1936cd / 12657196408"
grep -E "maximum resident set size" "$HERE/work/build-reg10.log" | tail -1
echo "=== DONE $(date -u +%FT%TZ) ==="
