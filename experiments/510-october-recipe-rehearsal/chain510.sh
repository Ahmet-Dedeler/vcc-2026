#!/usr/bin/env bash
# 510: the October recipe as decided on 2026-09-27, rehearsed on efr (238's ef-bundle-r) with the same-line call FORCED
# the 452b way (D=K562, E=HCT116, F=A549: in-mix banks on two contexts, a sparse one on the third). 487's same-line setting:
# args 5-6 = 0.5 and "" (w 0.5, no pt, programme on). Plus every other candidate switch: DEF_DAY_LODO=1 +
# DEF_DAY_LODO_HI_PCT=96.64697522769474, DEF_DAY_SHARPEN_PT=stercile:2.5,2.0,1.75, DEF_DAY_PROGRAMME_SOURCE=h1_2025_full,
# DEF_DAY_PRIOR_GATE=le1. No known answer exists for this configuration: the checks are DONE, postbuild PASS, every final
# rho <= 1, the builder restored, and the audits recorded (same-line w per context, gate counts, cap, terciles).
# Serial after chain497e (process exit); needs 497's def_day.sh (939e32ec) at HEAD. Never submits.
set -uo pipefail
REPO=~/Code/virtual-cell-challenge-2026; cd "$REPO" || exit 1
X=experiments/510-october-recipe-rehearsal; W=$X/work; LOG=$W/chain510.log; DD=scripts/def_day.sh; BU=scripts/make_submission_pm_cd4.py
EFR=experiments/238-essential-free-panel/work/ef-bundle-r
log(){ printf '%s %s\n' "$(date -u +%FT%TZ)" "$*" >> "$LOG"; }
echo $$ > $W/chain510.pid
log "armed (pid $$): waiting for chain497e's process to exit, 497's def_day.sh at HEAD, builder == HEAD"
P497E="$(cat experiments/497-october-switches/work/chain497e.pid 2>/dev/null)"
until ! ps -p "${P497E:-1}" -o args= 2>/dev/null | /usr/bin/grep -q chain497e \
      && ! ps -axo comm,args | awk '$1 ~ /^bash$/ && /chain4(7[69]b?|9[5-8])\.sh|chain496p\.sh|chain497e\.sh|run_ka_4|run_arm_generic\.sh|def_day\.sh|run_dump4/' | /usr/bin/grep -q . \
      && [ "$(shasum -a 256 $DD | cut -c1-16)" = e29e793a7964df03 ] \
      && [ "$(shasum -a 256 $DD | cut -d' ' -f1)" = "$(git show HEAD:$DD | shasum -a 256 | cut -d' ' -f1)" ] \
      && [ "$(shasum -a 256 $BU | cut -d' ' -f1)" = "$(git show HEAD:$BU | shasum -a 256 | cut -d' ' -f1)" ]; do
  [ "$(date -u +%Y%m%d%H%M)" -ge 202610021800 ] && { log "prerequisites not met by 18:00Z 10-02 -- stopping"; log "END"; exit 1; }; sleep 60; done
log "released: def_day $(shasum -a 256 $DD | cut -c1-16), builder $(shasum -a 256 $BU | cut -c1-16)"
rm -rf $W/efr510
DEF_DAY_FORCE="D=K562_HAEMATOPOIETIC_AND_LYMPHOID_TISSUE,E=HCT116_LARGE_INTESTINE,F=A549_LUNG" \
DEF_DAY_LODO=1 DEF_DAY_LODO_HI_PCT=96.64697522769474 DEF_DAY_SHARPEN_PT=stercile:2.5,2.0,1.75 \
DEF_DAY_PROGRAMME_SOURCE=h1_2025_full DEF_DAY_PRIOR_GATE=le1 \
  bash $DD $EFR $W/efr510 510efr 1.75 0.5 "" 4 > $W/efr510.nohup 2>&1; RC=$?
log "def_day rc=$RC: $(/usr/bin/grep -E 'REHEARSAL|same-line mode|LODO key ON|per-target sharpening|prior gate|CONVERGED at|NOT converged|shipping the nearest|final:|postbuild|STOP|DONE' $W/efr510/def_day.log | tr '\n' ' ' | cut -c1-1000)"
if /usr/bin/grep -q "DONE. Next" $W/efr510/def_day.log; then
  FIN="$(/usr/bin/grep -o 'final: [^ ]* factors [^ ]*' $W/efr510/def_day.log | tail -1)"
  log "final rho: $(uv run python scripts/def_day_helpers.py rho "$W/efr510/$(echo "$FIN" | awk '{print $2}').ntclevel.json" 2>&1 | tail -1)"
  log "same-line artifact: $(python3 -c 'import numpy as np,sys; z=np.load(sys.argv[1]); print({k: (z[k].tolist() if z[k].size < 4 else z[k].shape) for k in z.files if k.startswith(("w_","contexts"))})' $W/efr510/sameline.npz 2>&1 | cut -c1-400)"
  log "meta: $(python3 -c 'import json,sys; c=json.load(open(sys.argv[1]))["config"]; print("sameline", json.dumps(c.get("sameline"))[:300], "| gate", json.dumps(c.get("de_prior_gate"))[:300], "| pt", c.get("bulk_sharpen_pt"))' $W/efr510/510efr_final.h5ad.meta.json 2>&1 | cut -c1-900)"
  log "PASS: the 09-27 October recipe completes on efr with the forced in-mix same-line call"
else
  log "FAIL: the day did not complete (see $W/efr510/def_day.log)"
fi
rm -f $W/efr510/*.h5ad
[ "$(shasum -a 256 $BU | cut -d' ' -f1)" = "$(git show HEAD:$BU | shasum -a 256 | cut -d' ' -f1)" ] && log "builder restored" || log "WARNING builder differs from HEAD after the day"
log "END"
