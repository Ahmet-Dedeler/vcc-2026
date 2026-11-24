#!/usr/bin/env bash
# 446b = 446 re-levelled after ntclevel measured 446 at rho 0.978/0.965/0.969 (target 433k: 0.991/0.983/0.987; rho ~ f^1.86).
# usage: bash build446.sh "A=..,B=..,C=.."      Config diff vs 433k must be exactly bulk_ctx_scale + bulk_sharpen_q.
REPO=~/Code/virtual-cell-challenge-2026; cd "$REPO" || exit 1
SCALE="${1:?VCC_BULK_CTX_SCALE string}"
E=experiments/428-zero-gene-audit; X=experiments/446-q175-level; LOG="$X/work/build446b.log"
log(){ printf '%s %s\n' "$(date -u +%FT%TZ)" "$*" | tee -a "$LOG"; }
B=63dfd7d0e402738854072e57c81777a71fd773212877a6f26540a363b2aa34d3; EM=edb5b8dcbee519683b9d50c3556d4470ca6290594d4d1611be3d55c4fb7d41c7
REF="$REPO/submissions/433k-cap15-partial50-clean4q80-pertscale3-ballast-sharpen150-stackload-ctxlevel-mde010coarse-zr5.h5ad.meta.json"
D="$REPO/experiments/417-bulk-sharpen/work/dump_h1x2_cd4x15.npz"; DS=$(shasum -a 256 "$D" | cut -d' ' -f1)
log "446b start (q 1.75 levelled: $SCALE)"
( unset VCC_DE_MAG VCC_RC_AUDIT_DUMP
  export VCC_BULK_LOAD="$D:$DS" VCC_DE_MDE=0.1,0.1,0.3,coarse VCC_DE_ZERO_REPAIR=-5 VCC_BULK_CTX_SCALE="$SCALE"
  bash $E/run_arm_generic.sh 446b 1.75 "$REF" the-thermocycler-got-a-haircut \
    446b-cap15-partial50-clean4q80-pertscale3-ballast-sharpen175-stackload-ctxlevel-mde010coarse-zr5 "$B" "$EM" bulk_ctx_scale,bulk_sharpen_q ) >> "$X/work/build446b.nohup" 2>&1
log "446b: $(grep -E '^\[446b\] (PASS|FAIL|CLIFF)|audit|package sha|DONE$|stopping' $E/work/build_446b/run446b.log | tr '\n' ' ' | cut -c1-700)"
rm -f $E/work/build_446b/the-thermocycler-got-a-haircut.h5ad
log "END"
