#!/usr/bin/env bash
# After 445's dumps end: 446b (q 1.75 at 433k's matched rho, corrected factors), then its two riders (each one variable vs 446b):
#   446bj = 446b with the stack+LODO-key dump: vs 446b exactly bulk_load, bulk_load_sha256
#   446bz = 446b + VCC_DE_MAG=dn_k250_zs:     vs 446b exactly de_mag
REPO=~/Code/virtual-cell-challenge-2026; cd "$REPO" || exit 1
E=experiments/428-zero-gene-audit; X=experiments/446-q175-level; LOG="$X/work/chain446b.log"
log(){ printf '%s %s\n' "$(date -u +%FT%TZ)" "$*" | tee -a "$LOG"; }
until grep -q " END$" experiments/445-dld1-bulk-union/work/chain445dump.log 2>/dev/null; do
  [ "$(date -u +%Y%m%d%H%M)" -ge 202609250600 ] && { log "445 dumps never ended by 06:00Z"; exit 1; }; sleep 30; done
log "445 dumps ended; building 446b"
bash $X/build446b.sh "A=1.1004,B=1.1517,C=1.1720"
grep -q " DONE$" $E/work/build_446b/run446b.log || { log "446b did not finish clean -- no riders"; exit 1; }
B=63dfd7d0e402738854072e57c81777a71fd773212877a6f26540a363b2aa34d3; EM=edb5b8dcbee519683b9d50c3556d4470ca6290594d4d1611be3d55c4fb7d41c7
REF="$REPO/submissions/446b-cap15-partial50-clean4q80-pertscale3-ballast-sharpen175-stackload-ctxlevel-mde010coarse-zr5.h5ad.meta.json"
SCALE="A=1.1004,B=1.1517,C=1.1720"
D="$REPO/experiments/417-bulk-sharpen/work/dump_h1x2_cd4x15.npz"; DS=$(shasum -a 256 "$D" | cut -d' ' -f1)
DL="$REPO/experiments/417-bulk-sharpen/work/dump_h1x2_cd4x15_lodo05.npz"; DLS=$(shasum -a 256 "$DL" | cut -d' ' -f1)
[ "$DLS" = "00ccfcba1b67fa24ff743ee0c64067e01158147ffe49c5b8ca37cef67600b667" ] || { log "LODO dump sha $DLS unexpected"; exit 1; }
log "446bj start"
( unset VCC_DE_MAG VCC_RC_AUDIT_DUMP
  export VCC_BULK_LOAD="$DL:$DLS" VCC_DE_MDE=0.1,0.1,0.3,coarse VCC_DE_ZERO_REPAIR=-5 VCC_BULK_CTX_SCALE="$SCALE"
  bash $E/run_arm_generic.sh 446bj 1.75 "$REF" the-spectrophotometer-wrote-a-sonnet \
    446bj-cap15-partial50-clean4q80-pertscale3-ballast-sharpen175-stacklodo05load-ctxlevel-mde010coarse-zr5 "$B" "$EM" bulk_load,bulk_load_sha256 ) >> "$X/work/build446bj.nohup" 2>&1
log "446bj: $(grep -E '^\[446bj\] (PASS|FAIL|CLIFF)|audit|package sha|DONE$|stopping' $E/work/build_446bj/run446bj.log | tr '\n' ' ' | cut -c1-700)"
rm -f $E/work/build_446bj/the-spectrophotometer-wrote-a-sonnet.h5ad
log "446bz start"
( unset VCC_RC_AUDIT_DUMP
  export VCC_BULK_LOAD="$D:$DS" VCC_DE_MDE=0.1,0.1,0.3,coarse VCC_DE_ZERO_REPAIR=-5 VCC_BULK_CTX_SCALE="$SCALE" VCC_DE_MAG=dn_k250_zs
  bash $E/run_arm_generic.sh 446bz 1.75 "$REF" the-vortex-mixer-started-journaling \
    446bz-cap15-partial50-clean4q80-pertscale3-ballast-sharpen175-stackload-ctxlevel-mde010coarse-zr5-zs250 "$B" "$EM" de_mag ) >> "$X/work/build446bz.nohup" 2>&1
log "446bz: $(grep -E '^\[446bz\] (PASS|FAIL|CLIFF)|audit|package sha|DONE$|stopping' $E/work/build_446bz/run446bz.log | tr '\n' ' ' | cut -c1-700)"
rm -f $E/work/build_446bz/the-vortex-mixer-started-journaling.h5ad
log "END"
