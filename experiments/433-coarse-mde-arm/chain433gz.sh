#!/usr/bin/env bash
# After the 432 worker releases the builder (its claude process has exited AND no builder has run for 10 min):
#   433g = 428g with coarse (417g = 417a + re-level A1.002/B1.0408/C1.0648; coarse MDE; zero repair): vs 428g exactly de_mde + emitter sha.
#          The 09-25 fallback if 417e's stack reads badly.
#   433z = 433h + VCC_DE_MAG=dn_k250_zs (429's zs pass; 435: stacks on coarse, nmae -0.00183 on H1): vs 433h exactly de_mag.
# The wait loop is capped by a timestamp (21:00Z 09-24), never an iteration count.
REPO=~/Code/virtual-cell-challenge-2026; cd "$REPO" || exit 1
E=experiments/428-zero-gene-audit; X=experiments/433-coarse-mde-arm; LOG="$X/work/chain433gz.log"
log(){ printf '%s %s\n' "$(date -u +%FT%TZ)" "$*" | tee -a "$LOG"; }
WPID="$(cat experiments/432-def-rehearsal-new-knobs/claude-run.pid)"
log "WAITING for worker 432 (pid $WPID) to exit and the builder to idle 10 min (pid $$)"
idle=0
while true; do
  [ "$(date -u +%Y%m%d%H%M)" -ge 202609242100 ] && { log "gave up at 21:00Z 09-24"; exit 1; }
  if ps -p "$WPID" >/dev/null 2>&1 || ps -axo comm,args | awk '$1 ~ /[Pp]ython/ && /make_submission_pm_cd4|def_pipeline/' | grep -q .; then idle=0; else idle=$((idle+1)); fi
  [ $idle -ge 20 ] && break
  sleep 30
done
B=2c90d9e1cc2c19fd427da35ffe6a6d1306d3097382ff7a343233423cb710573f; EM=edb5b8dcbee519683b9d50c3556d4470ca6290594d4d1611be3d55c4fb7d41c7
[ "$(shasum -a 256 scripts/make_submission_pm_cd4.py | cut -d' ' -f1)" = "$B" ] || { log "builder is not $B -- stopping"; exit 1; }
log "433g start"
( unset VCC_DE_MAG VCC_RC_AUDIT_DUMP VCC_BULK_LOAD
  export VCC_DE_MDE=0.1,0.1,0.3,coarse VCC_DE_ZERO_REPAIR=-5 VCC_BULK_CTX_SCALE="A=1.002,B=1.0408,C=1.0648"
  bash $E/run_arm_generic.sh 433g 1.25 "$REPO/submissions/428g-cap15-partial50-clean4q80-pertscale3-ballast-sharpen125-ctxlevel-mde010-zr5.h5ad.meta.json" \
    the-autoclave-wrote-a-sonnet 433g-cap15-partial50-clean4q80-pertscale3-ballast-sharpen125-ctxlevel-mde010coarse-zr5 "$B" "$EM" de_mde,de_sign_emitter_sha256 ) >> "$X/work/build433g.nohup" 2>&1
log "433g: $(grep -E '^\[433g\] (PASS|FAIL|CLIFF)|audit|package sha|DONE$|stopping' $E/work/build_433g/run433g.log | tr '\n' ' ' | cut -c1-700)"
rm -f $E/work/build_433g/the-autoclave-wrote-a-sonnet.h5ad
D="$REPO/experiments/417-bulk-sharpen/work/dump_h1x2_cd4x15.npz"; DS=$(shasum -a 256 "$D" | cut -d' ' -f1)
log "433z start"
( unset VCC_RC_AUDIT_DUMP
  export VCC_BULK_LOAD="$D:$DS" VCC_DE_MDE=0.1,0.1,0.3,coarse VCC_DE_ZERO_REPAIR=-5 VCC_BULK_CTX_SCALE="A=1.0088,B=1.05,C=1.0736" VCC_DE_MAG=dn_k250_zs
  bash $E/run_arm_generic.sh 433z 1.25 "$REPO/submissions/433h-cap15-partial50-clean4q80-pertscale3-ballast-sharpen125-stackload-mde010coarse-zr5-ctxlevel.h5ad.meta.json" \
    the-vortexer-joined-a-choir 433z-cap15-partial50-clean4q80-pertscale3-ballast-sharpen125-stackload-mde010coarse-zr5-ctxlevel-zs250 "$B" "$EM" de_mag ) >> "$X/work/build433z.nohup" 2>&1
log "433z: $(grep -E '^\[433z\] (PASS|FAIL|CLIFF)|audit|package sha|DONE$|stopping' $E/work/build_433z/run433z.log | tr '\n' ' ' | cut -c1-700)"
rm -f $E/work/build_433z/the-vortexer-joined-a-choir.h5ad
log "END"
