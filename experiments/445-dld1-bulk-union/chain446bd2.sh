#!/usr/bin/env bash
# 446bd2 = 446bd re-levelled to 446b's matched rho. ntclevel read 446bd at rho 0.9649/0.9573/0.9620 vs 446b's
# 0.9908/0.9827/0.9873 (union 2 adds a residual and the global preserve rescales the base: -2.6 % render energy).
# Factors f_446b * (rho_446b/rho_446bd)^(1/1.86). Config diff vs 446b: bulk_load, bulk_load_sha256, bulk_ctx_scale.
REPO=~/Code/virtual-cell-challenge-2026; cd "$REPO" || exit 1
E=experiments/428-zero-gene-audit; X=experiments/445-dld1-bulk-union; LOG="$X/work/chain446bd2.log"
log(){ printf '%s %s\n' "$(date -u +%FT%TZ)" "$*" | tee -a "$LOG"; }
until grep -q " END$" experiments/446-q175-level/work/chain446bj2.log 2>/dev/null; do
  [ "$(date -u +%Y%m%d%H%M)" -ge 202609250900 ] && { log "chain446bj2 never ended by 09:00Z"; exit 1; }; sleep 30; done
B=63dfd7d0e402738854072e57c81777a71fd773212877a6f26540a363b2aa34d3; EM=edb5b8dcbee519683b9d50c3556d4470ca6290594d4d1611be3d55c4fb7d41c7
REF="$REPO/submissions/446b-cap15-partial50-clean4q80-pertscale3-ballast-sharpen175-stackload-ctxlevel-mde010coarse-zr5.h5ad.meta.json"
DD="$REPO/$X/work/dump_h1x2_cd4x15_dld1u2g020.npz"; DDS=$(shasum -a 256 "$DD" | cut -d' ' -f1)
[ "$DDS" = "575b4fdae8635cb0dc4d2e39fc74cce6e43c313d17046a1776e20e43c2e6fccc" ] || { log "DLD-1 dump sha $DDS unexpected"; exit 1; }
SCALE="A=1.1162,B=1.1680,C=1.1885"
log "446bd2 start ($SCALE)"
( unset VCC_DE_MAG VCC_RC_AUDIT_DUMP
  export VCC_BULK_LOAD="$DD:$DDS" VCC_DE_MDE=0.1,0.1,0.3,coarse VCC_DE_ZERO_REPAIR=-5 VCC_BULK_CTX_SCALE="$SCALE"
  bash $E/run_arm_generic.sh 446bd2 1.75 "$REF" the-microscope-took-up-birdwatching \
    446bd2-cap15-partial50-clean4q80-pertscale3-ballast-sharpen175-stackdld1u2load-ctxlevelm-mde010coarse-zr5 "$B" "$EM" bulk_ctx_scale,bulk_load,bulk_load_sha256 ) >> "$X/work/build446bd2.nohup" 2>&1
log "446bd2: $(grep -E '^\[446bd2\] (PASS|FAIL|CLIFF)|audit|package sha|DONE$|stopping' $E/work/build_446bd2/run446bd2.log | tr '\n' ' ' | cut -c1-700)"
rm -f $E/work/build_446bd2/the-microscope-took-up-birdwatching.h5ad
log "END"
