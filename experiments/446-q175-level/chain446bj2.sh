#!/usr/bin/env bash
# 446bj2 = 446bj re-levelled to 446b's matched rho. ntclevel read 446bj at rho 0.9978/0.9916/0.9935 vs 446b's
# 0.9908/0.9827/0.9873 (the LODO key adds 0.6-0.9 % render energy at q 1.75; 187's claims put B at ~0.999 of its cliff).
# Factors f_446b * (rho_446b/rho_446bj)^(1/1.86). Config diff vs 446b: bulk_load, bulk_load_sha256, bulk_ctx_scale.
REPO=~/Code/virtual-cell-challenge-2026; cd "$REPO" || exit 1
E=experiments/428-zero-gene-audit; X=experiments/446-q175-level; LOG="$X/work/chain446bj2.log"
log(){ printf '%s %s\n' "$(date -u +%FT%TZ)" "$*" | tee -a "$LOG"; }
until grep -q " END$" experiments/445-dld1-bulk-union/work/chain446bd.log 2>/dev/null; do
  [ "$(date -u +%Y%m%d%H%M)" -ge 202609250800 ] && { log "chain446bd never ended by 08:00Z"; exit 1; }; sleep 30; done
B=63dfd7d0e402738854072e57c81777a71fd773212877a6f26540a363b2aa34d3; EM=edb5b8dcbee519683b9d50c3556d4470ca6290594d4d1611be3d55c4fb7d41c7
REF="$REPO/submissions/446b-cap15-partial50-clean4q80-pertscale3-ballast-sharpen175-stackload-ctxlevel-mde010coarse-zr5.h5ad.meta.json"
DL="$REPO/experiments/417-bulk-sharpen/work/dump_h1x2_cd4x15_lodo05.npz"; DLS=$(shasum -a 256 "$DL" | cut -d' ' -f1)
[ "$DLS" = "00ccfcba1b67fa24ff743ee0c64067e01158147ffe49c5b8ca37cef67600b667" ] || { log "LODO dump sha $DLS unexpected"; exit 1; }
SCALE="A=1.0962,B=1.1461,C=1.1681"
log "446bj2 start ($SCALE)"
( unset VCC_DE_MAG VCC_RC_AUDIT_DUMP
  export VCC_BULK_LOAD="$DL:$DLS" VCC_DE_MDE=0.1,0.1,0.3,coarse VCC_DE_ZERO_REPAIR=-5 VCC_BULK_CTX_SCALE="$SCALE"
  bash $E/run_arm_generic.sh 446bj2 1.75 "$REF" the-centrifuge-joined-a-book-club \
    446bj2-cap15-partial50-clean4q80-pertscale3-ballast-sharpen175-stacklodo05load-ctxlevelm-mde010coarse-zr5 "$B" "$EM" bulk_ctx_scale,bulk_load,bulk_load_sha256 ) >> "$X/work/build446bj2.nohup" 2>&1
log "446bj2: $(grep -E '^\[446bj2\] (PASS|FAIL|CLIFF)|audit|package sha|DONE$|stopping' $E/work/build_446bj2/run446bj2.log | tr '\n' ' ' | cut -c1-700)"
rm -f $E/work/build_446bj2/the-centrifuge-joined-a-book-club.h5ad
log "END"
