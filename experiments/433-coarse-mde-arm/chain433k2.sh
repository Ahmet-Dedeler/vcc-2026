#!/usr/bin/env bash
# 09-26 candidates on 433k's bulk (q 1.5 + 417k's re-level), built before the 09-25 reads so either branch is ready:
#   433kz = 433k + VCC_DE_MAG=dn_k250_zs            : vs 433k exactly de_mag
#   433kj = 433k with the stack+LODO-key dump (417j's): vs 433k exactly bulk_load, bulk_load_sha256
# First arms built through the NEW driver (5463fe28) and builder (63dfd7d0, T6 guard). Sequential; exits if a gate fails.
REPO=~/Code/virtual-cell-challenge-2026; cd "$REPO" || exit 1
E=experiments/428-zero-gene-audit; X=experiments/433-coarse-mde-arm; LOG="$X/work/chain433k2.log"
log(){ printf '%s %s\n' "$(date -u +%FT%TZ)" "$*" | tee -a "$LOG"; }
B=63dfd7d0e402738854072e57c81777a71fd773212877a6f26540a363b2aa34d3; EM=edb5b8dcbee519683b9d50c3556d4470ca6290594d4d1611be3d55c4fb7d41c7
REF="$REPO/submissions/433k-cap15-partial50-clean4q80-pertscale3-ballast-sharpen150-stackload-ctxlevel-mde010coarse-zr5.h5ad.meta.json"
D="$REPO/experiments/417-bulk-sharpen/work/dump_h1x2_cd4x15.npz"; DS=$(shasum -a 256 "$D" | cut -d' ' -f1)
DL="$REPO/experiments/417-bulk-sharpen/work/dump_h1x2_cd4x15_lodo05.npz"; DLS=$(shasum -a 256 "$DL" | cut -d' ' -f1)
[ "$DLS" = "00ccfcba1b67fa24ff743ee0c64067e01158147ffe49c5b8ca37cef67600b667" ] || { log "LODO dump sha $DLS unexpected"; exit 1; }
log "433kz start"
( unset VCC_RC_AUDIT_DUMP
  export VCC_BULK_LOAD="$D:$DS" VCC_DE_MDE=0.1,0.1,0.3,coarse VCC_DE_ZERO_REPAIR=-5 VCC_BULK_CTX_SCALE="A=1.0464,B=1.092,C=1.114" VCC_DE_MAG=dn_k250_zs
  bash $E/run_arm_generic.sh 433kz 1.5 "$REF" the-lab-coat-started-a-podcast \
    433kz-cap15-partial50-clean4q80-pertscale3-ballast-sharpen150-stackload-ctxlevel-mde010coarse-zr5-zs250 "$B" "$EM" de_mag ) >> "$X/work/build433kz.nohup" 2>&1
log "433kz: $(grep -E '^\[433kz\] (PASS|FAIL|CLIFF)|audit|package sha|DONE$|stopping' $E/work/build_433kz/run433kz.log | tr '\n' ' ' | cut -c1-700)"
rm -f $E/work/build_433kz/the-lab-coat-started-a-podcast.h5ad
grep -q " DONE$" $E/work/build_433kz/run433kz.log || { log "433kz did not finish clean -- stopping"; exit 1; }
log "433kj start"
( unset VCC_DE_MAG VCC_RC_AUDIT_DUMP
  export VCC_BULK_LOAD="$DL:$DLS" VCC_DE_MDE=0.1,0.1,0.3,coarse VCC_DE_ZERO_REPAIR=-5 VCC_BULK_CTX_SCALE="A=1.0464,B=1.092,C=1.114"
  bash $E/run_arm_generic.sh 433kj 1.5 "$REF" the-petri-dish-went-viral \
    433kj-cap15-partial50-clean4q80-pertscale3-ballast-sharpen150-stacklodo05load-ctxlevel-mde010coarse-zr5 "$B" "$EM" bulk_load,bulk_load_sha256 ) >> "$X/work/build433kj.nohup" 2>&1
log "433kj: $(grep -E '^\[433kj\] (PASS|FAIL|CLIFF)|audit|package sha|DONE$|stopping' $E/work/build_433kj/run433kj.log | tr '\n' ' ' | cut -c1-700)"
rm -f $E/work/build_433kj/the-petri-dish-went-viral.h5ad
log "END"
