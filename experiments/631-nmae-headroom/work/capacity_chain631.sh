#!/bin/bash
# 631 A/B/C capacity for the DNMAG rider: restore 446bj2z, run 630 v3 then DN d_m200 on 60 targets/context, delete the
# restored h5ad. One gate602 ticket per heavy run; the delete runs only if the restore made the file. END last.
set -u
R=~/Code/virtual-cell-challenge-2026
W=$R/experiments/631-nmae-headroom/work
H5=$W/restored_446bj2z.h5ad
G=experiments/602-council-gate/gate602.sh
cd $R
echo "$(date -u +%FT%TZ) restore start"
[ -f "$H5" ] || bash $G 631-restore 4 -- uv run python scripts/restore_prediction_h5ad.py submissions/446bj2z-cap15-partial50-clean4q80-pertscale3-ballast-sharpen175-stacklodo05load-ctxlevelm-mde010coarse-zr5-zs250.vcc "$H5"
echo "$(date -u +%FT%TZ) restore rc $? size $(du -h "$H5" 2>/dev/null | cut -f1)"
bash $G 631-capacity 6 -- uv run python experiments/631-nmae-headroom/capacity631.py "$H5" 60 b_k100_q20r3p10 d_m200
echo "$(date -u +%FT%TZ) capacity rc $?"
[ -f "$H5" ] && rm -f "$H5" && echo "$(date -u +%FT%TZ) removed restored h5ad"
echo "$(date -u +%FT%TZ) END"
