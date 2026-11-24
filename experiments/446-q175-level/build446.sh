#!/usr/bin/env bash
# 446 = 433k with q 1.5 -> 1.75 and q 1.75's OWN re-level (factors from ntclevel417.py on the 446u measurement build).
# usage: bash build446.sh "A=..,B=..,C=.."      Config diff vs 433k must be exactly bulk_ctx_scale + bulk_sharpen_q.
REPO=~/Code/virtual-cell-challenge-2026; cd "$REPO" || exit 1
SCALE="${1:?VCC_BULK_CTX_SCALE string}"
E=experiments/428-zero-gene-audit; X=experiments/446-q175-level; LOG="$X/work/build446.log"
log(){ printf '%s %s\n' "$(date -u +%FT%TZ)" "$*" | tee -a "$LOG"; }
B=63dfd7d0e402738854072e57c81777a71fd773212877a6f26540a363b2aa34d3; EM=edb5b8dcbee519683b9d50c3556d4470ca6290594d4d1611be3d55c4fb7d41c7
REF="$REPO/submissions/433k-cap15-partial50-clean4q80-pertscale3-ballast-sharpen150-stackload-ctxlevel-mde010coarse-zr5.h5ad.meta.json"
D="$REPO/experiments/417-bulk-sharpen/work/dump_h1x2_cd4x15.npz"; DS=$(shasum -a 256 "$D" | cut -d' ' -f1)
log "446 start (q 1.75 levelled: $SCALE)"
( unset VCC_DE_MAG VCC_RC_AUDIT_DUMP
  export VCC_BULK_LOAD="$D:$DS" VCC_DE_MDE=0.1,0.1,0.3,coarse VCC_DE_ZERO_REPAIR=-5 VCC_BULK_CTX_SCALE="$SCALE"
  bash $E/run_arm_generic.sh 446 1.75 "$REF" the-sequencer-joined-a-choir \
    446-cap15-partial50-clean4q80-pertscale3-ballast-sharpen175-stackload-ctxlevel-mde010coarse-zr5 "$B" "$EM" bulk_ctx_scale,bulk_sharpen_q ) >> "$X/work/build446.nohup" 2>&1
log "446: $(grep -E '^\[446\] (PASS|FAIL|CLIFF)|audit|package sha|DONE$|stopping' $E/work/build_446/run446.log | tr '\n' ' ' | cut -c1-700)"
rm -f $E/work/build_446/the-sequencer-joined-a-choir.h5ad
log "END"
