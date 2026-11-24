#!/bin/bash
# 630 capacity pre-check: restore 446bj2z's package (carries zs), run the module on 60 targets/context, then delete the
# restored h5ad (derived; the .vcc stays). Every step logs; the delete runs only if the restore produced the file.
set -u
R=~/Code/virtual-cell-challenge-2026
W=$R/experiments/630-rank-neutral-up/work
H5=$W/restored_446bj2z.h5ad
cd $R
echo "$(date -u +%FT%TZ) restore start"
[ -f "$H5" ] || uv run python scripts/restore_prediction_h5ad.py submissions/446bj2z-cap15-partial50-clean4q80-pertscale3-ballast-sharpen175-stacklodo05load-ctxlevelm-mde010coarse-zr5-zs250.vcc "$H5"
echo "$(date -u +%FT%TZ) restore rc $? size $(du -h "$H5" 2>/dev/null | cut -f1)"
uv run python experiments/630-rank-neutral-up/capacity630.py "$H5" 60 "${VARIANT:-b_k100}"
echo "$(date -u +%FT%TZ) capacity rc $?"
[ -f "$H5" ] && rm -f "$H5" && echo "$(date -u +%FT%TZ) removed restored h5ad"
echo "$(date -u +%FT%TZ) END"
