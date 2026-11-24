#!/usr/bin/env bash
# q 1.75 levelled on 433k's config (the next step of the dose curve that paid +0.0046 at q 1.5 on 09-25):
#   433q  = 433k with q 1.75 and NO re-level (the energy source for the re-level; packaged but not for submission)
#   433qL = 433k with q 1.75 + its own EV-optimal re-level (ctxlevel417.py on 433q's energies at 433k's board raw mse)
# Both: vs 433k the config diff is exactly bulk_sharpen_q + bulk_ctx_scale.
REPO=~/Code/virtual-cell-challenge-2026; cd "$REPO" || exit 1
E=experiments/428-zero-gene-audit; X=experiments/433-coarse-mde-arm; LOG="$X/work/chain433q.log"
log(){ printf '%s %s\n' "$(date -u +%FT%TZ)" "$*" | tee -a "$LOG"; }
B=63dfd7d0e402738854072e57c81777a71fd773212877a6f26540a363b2aa34d3; EM=edb5b8dcbee519683b9d50c3556d4470ca6290594d4d1611be3d55c4fb7d41c7
REF="$REPO/submissions/433k-cap15-partial50-clean4q80-pertscale3-ballast-sharpen150-stackload-ctxlevel-mde010coarse-zr5.h5ad.meta.json"
D="$REPO/experiments/417-bulk-sharpen/work/dump_h1x2_cd4x15.npz"; DS=$(shasum -a 256 "$D" | cut -d' ' -f1)
log "433q start (q 1.75, no re-level)"
( unset VCC_DE_MAG VCC_RC_AUDIT_DUMP VCC_BULK_CTX_SCALE
  export VCC_BULK_LOAD="$D:$DS" VCC_DE_MDE=0.1,0.1,0.3,coarse VCC_DE_ZERO_REPAIR=-5
  bash $E/run_arm_generic.sh 433q 1.75 "$REF" the-spectrophotometer-took-a-nap \
    433q-cap15-partial50-clean4q80-pertscale3-ballast-sharpen175-stackload-mde010coarse-zr5 "$B" "$EM" bulk_sharpen_q,bulk_ctx_scale ) >> "$X/work/build433q.nohup" 2>&1
log "433q: $(grep -E '^\[433q\] (PASS|FAIL|CLIFF)|audit|package sha|DONE$|stopping' $E/work/build_433q/run433q.log | tr '\n' ' ' | cut -c1-600)"
grep -q "^\[433q\] PASS" $E/work/build_433q/run433q.log || { log "433q config diff did not pass -- stopping"; exit 1; }
LV="$(uv run python experiments/417-bulk-sharpen/ctxlevel417.py the-spectrophotometer-took-a-nap.h5ad 0.648523 2>&1)"
echo "$LV" >> "$LOG"
CTX="$(echo "$LV" | grep -o 'VCC_BULK_CTX_SCALE=[A-Z0-9=.,]*' | tail -1 | cut -d= -f2-)"
log "re-level string for q 1.75: ${CTX:-<none>}"
[ -n "$CTX" ] || { log "no re-level string -- stopping"; exit 1; }
rm -f $E/work/build_433q/the-spectrophotometer-took-a-nap.h5ad
log "433qL start (q 1.75 + $CTX)"
( unset VCC_DE_MAG VCC_RC_AUDIT_DUMP
  export VCC_BULK_LOAD="$D:$DS" VCC_DE_MDE=0.1,0.1,0.3,coarse VCC_DE_ZERO_REPAIR=-5 VCC_BULK_CTX_SCALE="$CTX"
  bash $E/run_arm_generic.sh 433qL 1.75 "$REF" the-spectrophotometer-woke-up \
    433qL-cap15-partial50-clean4q80-pertscale3-ballast-sharpen175-stackload-ctxlevel-mde010coarse-zr5 "$B" "$EM" bulk_sharpen_q,bulk_ctx_scale ) >> "$X/work/build433qL.nohup" 2>&1
log "433qL: $(grep -E '^\[433qL\] (PASS|FAIL|CLIFF)|audit|package sha|DONE$|stopping' $E/work/build_433qL/run433qL.log | tr '\n' ' ' | cut -c1-600)"
rm -f $E/work/build_433qL/the-spectrophotometer-woke-up.h5ad
log "END"
