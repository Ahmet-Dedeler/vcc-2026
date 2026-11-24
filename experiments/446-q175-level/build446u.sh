#!/usr/bin/env bash
# 446u = 433k with q 1.5 -> 1.75 and NO per-context scale: a MEASUREMENT build (never submitted). Its render energy
# (cliffcheck414) gives the q 1.75 re-level factors the way 417b's artifact gave 417k's. Config diff vs 433k must be
# exactly bulk_ctx_scale + bulk_sharpen_q. The .vcc is deleted after the cliff check has been recorded.
REPO=~/Code/virtual-cell-challenge-2026; cd "$REPO" || exit 1
E=experiments/428-zero-gene-audit; X=experiments/446-q175-level; LOG="$X/work/build446u.log"
log(){ printf '%s %s\n' "$(date -u +%FT%TZ)" "$*" | tee -a "$LOG"; }
B=63dfd7d0e402738854072e57c81777a71fd773212877a6f26540a363b2aa34d3; EM=edb5b8dcbee519683b9d50c3556d4470ca6290594d4d1611be3d55c4fb7d41c7
REF="$REPO/submissions/433k-cap15-partial50-clean4q80-pertscale3-ballast-sharpen150-stackload-ctxlevel-mde010coarse-zr5.h5ad.meta.json"
D="$REPO/experiments/417-bulk-sharpen/work/dump_h1x2_cd4x15.npz"; DS=$(shasum -a 256 "$D" | cut -d' ' -f1)
log "446u start (q 1.75, unlevelled, measurement only)"
( unset VCC_DE_MAG VCC_RC_AUDIT_DUMP VCC_BULK_CTX_SCALE
  export VCC_BULK_LOAD="$D:$DS" VCC_DE_MDE=0.1,0.1,0.3,coarse VCC_DE_ZERO_REPAIR=-5
  bash $E/run_arm_generic.sh 446u 1.75 "$REF" the-autoclave-read-the-manual \
    446u-measure-q175-unlevelled "$B" "$EM" bulk_ctx_scale,bulk_sharpen_q ) >> "$X/work/build446u.nohup" 2>&1
log "446u: $(grep -E '^\[446u\] (PASS|FAIL|CLIFF)|audit|package sha|DONE$|stopping' $E/work/build_446u/run446u.log | tr '\n' ' ' | cut -c1-700)"
rm -f "$REPO/submissions/446u-measure-q175-unlevelled.vcc"
log "measurement package deleted (never for submission); h5ad kept for re-level until factors are recorded"
log "END"
