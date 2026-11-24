#!/usr/bin/env bash
# 433j = 433h with the stack dumped under 420's LODO key at clip(0.5, 3) (417j's dump): vs 433h exactly bulk_load +
# bulk_load_sha256. 420 priced the key step at about +0.003; 417j (the same step on 417h) is built but unsent.
# Waits for chain433gz.sh's END (the 433g/433z builds); wait capped by timestamp (22:00Z 09-24).
REPO=~/Code/virtual-cell-challenge-2026; cd "$REPO" || exit 1
E=experiments/428-zero-gene-audit; X=experiments/433-coarse-mde-arm; LOG="$X/work/chain433j.log"
log(){ printf '%s %s\n' "$(date -u +%FT%TZ)" "$*" | tee -a "$LOG"; }
log "WAITING for chain433gz END (pid $$)"
until grep -q " END$" "$X/work/chain433gz.log" 2>/dev/null; do
  grep -qE "gave up|stopping" "$X/work/chain433gz.log" 2>/dev/null && { log "chain433gz did not finish clean -- stopping"; exit 1; }
  [ "$(date -u +%Y%m%d%H%M)" -ge 202609242200 ] && { log "gave up at 22:00Z 09-24"; exit 1; }; sleep 60; done
B=2c90d9e1cc2c19fd427da35ffe6a6d1306d3097382ff7a343233423cb710573f; EM=edb5b8dcbee519683b9d50c3556d4470ca6290594d4d1611be3d55c4fb7d41c7
[ "$(shasum -a 256 scripts/make_submission_pm_cd4.py | cut -d' ' -f1)" = "$B" ] || { log "builder is not $B -- stopping"; exit 1; }
D="$REPO/experiments/417-bulk-sharpen/work/dump_h1x2_cd4x15_lodo05.npz"; DS=$(shasum -a 256 "$D" | cut -d' ' -f1)
[ "$DS" = "00ccfcba1b67fa24ff743ee0c64067e01158147ffe49c5b8ca37cef67600b667" ] || { log "LODO dump sha $DS unexpected -- stopping"; exit 1; }
log "433j start"
( unset VCC_DE_MAG VCC_RC_AUDIT_DUMP
  export VCC_BULK_LOAD="$D:$DS" VCC_DE_MDE=0.1,0.1,0.3,coarse VCC_DE_ZERO_REPAIR=-5 VCC_BULK_CTX_SCALE="A=1.0088,B=1.05,C=1.0736"
  bash $E/run_arm_generic.sh 433j 1.25 "$REPO/submissions/433h-cap15-partial50-clean4q80-pertscale3-ballast-sharpen125-stackload-mde010coarse-zr5-ctxlevel.h5ad.meta.json" \
    the-pipette-tip-ran-for-office 433j-cap15-partial50-clean4q80-pertscale3-ballast-sharpen125-stacklodo05load-mde010coarse-zr5-ctxlevel "$B" "$EM" bulk_load,bulk_load_sha256 ) >> "$X/work/build433j.nohup" 2>&1
log "433j: $(grep -E '^\[433j\] (PASS|FAIL|CLIFF)|audit|package sha|DONE$|stopping' $E/work/build_433j/run433j.log | tr '\n' ' ' | cut -c1-700)"
rm -f $E/work/build_433j/the-pipette-tip-ran-for-office.h5ad
log "END"
