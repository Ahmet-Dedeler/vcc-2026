#!/usr/bin/env bash
# After 448's rehearsal ends and restores the builder: 446br2 = 446b with VCC_DE_ZERO_REPAIR=-2 (config diff exactly de_zero_repair).
REPO=~/Code/virtual-cell-challenge-2026; cd "$REPO" || exit 1
E=experiments/428-zero-gene-audit; X=experiments/449-repair-floor; LOG="$X/work/chain449.log"
log(){ printf '%s %s\n' "$(date -u +%FT%TZ)" "$*" | tee -a "$LOG"; }
until grep -q " END$" experiments/448-def-rehearsal-q175/work/run448.log 2>/dev/null && ! ps -axo comm,args | awk '$1 ~ /^bash$/ && /run448.sh/' | grep -q .; do
  [ "$(date -u +%Y%m%d%H%M)" -ge 202609251200 ] && { log "448 never ended by 12:00Z"; exit 1; }; sleep 30; done
B=b6c3f13141e707012588cc354f62758889db9c38c3757f9765ebea37a8a733f8; EM=edb5b8dcbee519683b9d50c3556d4470ca6290594d4d1611be3d55c4fb7d41c7
[ "$(shasum -a 256 scripts/make_submission_pm_cd4.py | cut -d' ' -f1)" = "$B" ] || { log "builder not restored to b6c3f131 -- stopping"; exit 1; }
REF="$REPO/submissions/446b-cap15-partial50-clean4q80-pertscale3-ballast-sharpen175-stackload-ctxlevel-mde010coarse-zr5.h5ad.meta.json"
D="$REPO/experiments/417-bulk-sharpen/work/dump_h1x2_cd4x15.npz"; DS=$(shasum -a 256 "$D" | cut -d' ' -f1)
log "446br2 start (floor -2)"
( unset VCC_DE_MAG VCC_RC_AUDIT_DUMP
  export VCC_BULK_LOAD="$D:$DS" VCC_DE_MDE=0.1,0.1,0.3,coarse VCC_DE_ZERO_REPAIR=-2 VCC_BULK_CTX_SCALE="A=1.1004,B=1.1517,C=1.1720"
  bash $E/run_arm_generic.sh 446br2 1.75 "$REF" the-pipette-rack-learned-origami \
    446br2-cap15-partial50-clean4q80-pertscale3-ballast-sharpen175-stackload-ctxlevel-mde010coarse-zr2 "$B" "$EM" de_zero_repair ) >> "$X/work/build446br2.nohup" 2>&1
log "446br2: $(grep -E '^\[446br2\] (PASS|FAIL|CLIFF)|audit|package sha|DONE$|stopping' $E/work/build_446br2/run446br2.log | tr '\n' ' ' | cut -c1-700)"
A="$E/work/build_446br2/the-pipette-rack-learned-origami.h5ad"
[ -s "$A" ] && uv run python $X/audit_extremes.py "$A" > "$X/work/audit-446br2-extremes.txt" 2>&1 && log "extremes after floor -2: $(head -1 $X/work/audit-446br2-extremes.txt)"
rm -f "$A"
log "END"
