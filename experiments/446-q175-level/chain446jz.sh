#!/usr/bin/env bash
# After 445's dumps end: the two 09-26 riders on 446's bulk, each one variable vs 446 (the prebuilt 433kz/433kj cover the
# branch where 446 does not go):
#   446j = 446 with the stack+LODO-key dump (417j's, clip 0.5/3): vs 446 exactly bulk_load, bulk_load_sha256
#   446z = 446 + VCC_DE_MAG=dn_k250_zs                          : vs 446 exactly de_mag
# 433j's energies matched 433h's to 0.1 %, so 446's re-level holds for 446j (its CLIFF line re-checks it).
REPO=~/Code/virtual-cell-challenge-2026; cd "$REPO" || exit 1
E=experiments/428-zero-gene-audit; X=experiments/446-q175-level; LOG="$X/work/chain446jz.log"
log(){ printf '%s %s\n' "$(date -u +%FT%TZ)" "$*" | tee -a "$LOG"; }
until grep -q " END$" experiments/445-dld1-bulk-union/work/chain445dump.log 2>/dev/null; do
  [ "$(date -u +%Y%m%d%H%M)" -ge 202609250600 ] && { log "445 dumps never ended by 06:00Z"; exit 1; }; sleep 30; done
grep -q " DONE$" $E/work/build_446/run446.log || { log "446 did not finish clean -- no riders"; exit 1; }
B=63dfd7d0e402738854072e57c81777a71fd773212877a6f26540a363b2aa34d3; EM=edb5b8dcbee519683b9d50c3556d4470ca6290594d4d1611be3d55c4fb7d41c7
REF="$REPO/submissions/446-cap15-partial50-clean4q80-pertscale3-ballast-sharpen175-stackload-ctxlevel-mde010coarse-zr5.h5ad.meta.json"
SCALE="$(uv run python -c "import json; s=json.load(open('$REF'))['config']['bulk_ctx_scale']; print(','.join(f'{k}={v}' for k,v in s.items()))")"
D="$REPO/experiments/417-bulk-sharpen/work/dump_h1x2_cd4x15.npz"; DS=$(shasum -a 256 "$D" | cut -d' ' -f1)
DL="$REPO/experiments/417-bulk-sharpen/work/dump_h1x2_cd4x15_lodo05.npz"; DLS=$(shasum -a 256 "$DL" | cut -d' ' -f1)
[ "$DLS" = "00ccfcba1b67fa24ff743ee0c64067e01158147ffe49c5b8ca37cef67600b667" ] || { log "LODO dump sha $DLS unexpected"; exit 1; }
log "446j start (scale $SCALE)"
( unset VCC_DE_MAG VCC_RC_AUDIT_DUMP
  export VCC_BULK_LOAD="$DL:$DLS" VCC_DE_MDE=0.1,0.1,0.3,coarse VCC_DE_ZERO_REPAIR=-5 VCC_BULK_CTX_SCALE="$SCALE"
  bash $E/run_arm_generic.sh 446j 1.75 "$REF" the-spectrophotometer-wrote-a-sonnet \
    446j-cap15-partial50-clean4q80-pertscale3-ballast-sharpen175-stacklodo05load-ctxlevel-mde010coarse-zr5 "$B" "$EM" bulk_load,bulk_load_sha256 ) >> "$X/work/build446j.nohup" 2>&1
log "446j: $(grep -E '^\[446j\] (PASS|FAIL|CLIFF)|audit|package sha|DONE$|stopping' $E/work/build_446j/run446j.log | tr '\n' ' ' | cut -c1-700)"
rm -f $E/work/build_446j/the-spectrophotometer-wrote-a-sonnet.h5ad
log "446z start"
( unset VCC_RC_AUDIT_DUMP
  export VCC_BULK_LOAD="$D:$DS" VCC_DE_MDE=0.1,0.1,0.3,coarse VCC_DE_ZERO_REPAIR=-5 VCC_BULK_CTX_SCALE="$SCALE" VCC_DE_MAG=dn_k250_zs
  bash $E/run_arm_generic.sh 446z 1.75 "$REF" the-vortex-mixer-started-journaling \
    446z-cap15-partial50-clean4q80-pertscale3-ballast-sharpen175-stackload-ctxlevel-mde010coarse-zr5-zs250 "$B" "$EM" de_mag ) >> "$X/work/build446z.nohup" 2>&1
log "446z: $(grep -E '^\[446z\] (PASS|FAIL|CLIFF)|audit|package sha|DONE$|stopping' $E/work/build_446z/run446z.log | tr '\n' ' ' | cut -c1-700)"
rm -f $E/work/build_446z/the-vortex-mixer-started-journaling.h5ad
log "END"
