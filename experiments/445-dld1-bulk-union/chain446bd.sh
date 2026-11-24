#!/usr/bin/env bash
# After chain446b ends: 445r on 446b's bulk = 446b with VCC_BULK_LOAD swapped to the DLD-1 union-2 dump (575b4fda).
# Config diff vs 446b must be exactly bulk_load, bulk_load_sha256.
REPO=~/Code/virtual-cell-challenge-2026; cd "$REPO" || exit 1
E=experiments/428-zero-gene-audit; X=experiments/445-dld1-bulk-union; LOG="$X/work/chain446bd.log"
log(){ printf '%s %s\n' "$(date -u +%FT%TZ)" "$*" | tee -a "$LOG"; }
until grep -q " END$" experiments/446-q175-level/work/chain446b.log 2>/dev/null; do
  [ "$(date -u +%Y%m%d%H%M)" -ge 202609250800 ] && { log "chain446b never ended by 08:00Z"; exit 1; }; sleep 30; done
grep -q " DONE$" $E/work/build_446b/run446b.log || { log "446b did not finish clean -- no 446bd"; exit 1; }
B=63dfd7d0e402738854072e57c81777a71fd773212877a6f26540a363b2aa34d3; EM=edb5b8dcbee519683b9d50c3556d4470ca6290594d4d1611be3d55c4fb7d41c7
REF="$REPO/submissions/446b-cap15-partial50-clean4q80-pertscale3-ballast-sharpen175-stackload-ctxlevel-mde010coarse-zr5.h5ad.meta.json"
DD="$REPO/$X/work/dump_h1x2_cd4x15_dld1u2g020.npz"; DDS=$(shasum -a 256 "$DD" | cut -d' ' -f1)
[ "$DDS" = "575b4fdae8635cb0dc4d2e39fc74cce6e43c313d17046a1776e20e43c2e6fccc" ] || { log "DLD-1 dump sha $DDS unexpected"; exit 1; }
log "446bd start"
( unset VCC_DE_MAG VCC_RC_AUDIT_DUMP
  export VCC_BULK_LOAD="$DD:$DDS" VCC_DE_MDE=0.1,0.1,0.3,coarse VCC_DE_ZERO_REPAIR=-5 VCC_BULK_CTX_SCALE="A=1.1004,B=1.1517,C=1.1720"
  bash $E/run_arm_generic.sh 446bd 1.75 "$REF" the-pipette-learned-to-whistle \
    446bd-cap15-partial50-clean4q80-pertscale3-ballast-sharpen175-stackdld1u2load-ctxlevel-mde010coarse-zr5 "$B" "$EM" bulk_load,bulk_load_sha256 ) >> "$X/work/build446bd.nohup" 2>&1
log "446bd: $(grep -E '^\[446bd\] (PASS|FAIL|CLIFF)|audit|package sha|DONE$|stopping' $E/work/build_446bd/run446bd.log | tr '\n' ' ' | cut -c1-700)"
rm -f $E/work/build_446bd/the-pipette-learned-to-whistle.h5ad
log "END"
