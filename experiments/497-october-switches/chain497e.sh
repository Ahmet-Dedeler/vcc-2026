#!/usr/bin/env bash
# 497e: every candidate October switch at once on an October-like panel (efr = 238's ef-bundle-r, the panel 465 used):
#   DEF_DAY_LODO=1 + DEF_DAY_LODO_HI_PCT=96.64697522769474 (on efr the percentile lands at 8-16, 543) +
#   DEF_DAY_SHARPEN_PT=stercile:2.5,2.0,1.75 + DEF_DAY_PROGRAMME_SOURCE=h1_2025_full (flags pt) + DEF_DAY_PRIOR_GATE=le1.
# No known answer exists for this configuration: the checks are that the day completes (DONE), postbuild PASSes with the
# allowed diffs, every final rho <= 1, the builder is restored, and the cap / tercile / gate audits are recorded.
# Serial after chain496p (process exit); needs 497's def_day.sh at HEAD. Never submits.
set -uo pipefail
REPO=~/Code/virtual-cell-challenge-2026; cd "$REPO" || exit 1
X=experiments/497-october-switches; W=$X/work; LOG=$W/chain497e.log; DD=scripts/def_day.sh; BU=scripts/make_submission_pm_cd4.py
EFR=experiments/238-essential-free-panel/work/ef-bundle-r
log(){ printf '%s %s\n' "$(date -u +%FT%TZ)" "$*" >> "$LOG"; }
echo $$ > $W/chain497e.pid
log "armed (pid $$): waiting for chain496p's process to exit, 497's def_day.sh at HEAD, builder == HEAD"
P496="$(cat experiments/496-uncalled-sign-hook/work/chain496p.pid 2>/dev/null)"
until ! ps -p "${P496:-1}" -o args= 2>/dev/null | /usr/bin/grep -q chain496p \
      && ! ps -axo comm,args | awk '$1 ~ /^bash$/ && /chain4(7[69]b?|9[5-8])\.sh|chain496p\.sh|run_ka_4|run_arm_generic\.sh|def_day\.sh|run_dump4/' | /usr/bin/grep -q . \
      && [ "$(shasum -a 256 $DD | cut -c1-16)" = e29e793a7964df03 ] \
      && [ "$(shasum -a 256 $DD | cut -d' ' -f1)" = "$(git show HEAD:$DD | shasum -a 256 | cut -d' ' -f1)" ] \
      && [ "$(shasum -a 256 $BU | cut -d' ' -f1)" = "$(git show HEAD:$BU | shasum -a 256 | cut -d' ' -f1)" ]; do
  [ "$(date -u +%Y%m%d%H%M)" -ge 202610021400 ] && { log "prerequisites not met by 14:00Z 10-02 -- stopping"; log "END"; exit 1; }; sleep 60; done
FREE="$(vm_stat | awk '/Pages free|Pages inactive/ {gsub("\\.","",$NF); s+=$NF} END {printf "%d", s*16384/2^30}')"
log "released (free+inactive ${FREE} GiB): def_day $(shasum -a 256 $DD | cut -c1-16), builder $(shasum -a 256 $BU | cut -c1-16)"
rm -rf $W/efr497
DEF_DAY_LODO=1 DEF_DAY_LODO_HI_PCT=96.64697522769474 DEF_DAY_SHARPEN_PT=stercile:2.5,2.0,1.75 DEF_DAY_PROGRAMME_SOURCE=h1_2025_full \
  DEF_DAY_PRIOR_GATE=le1 bash $DD $EFR $W/efr497 497efr 1.75 auto pt 4 > $W/efr497.nohup 2>&1; RC=$?
log "def_day rc=$RC: $(/usr/bin/grep -E 'LODO key ON|per-target sharpening|prior gate|same-line mode|CONVERGED at|NOT converged|shipping the nearest|final:|postbuild|STOP|DONE' $W/efr497/def_day.log | tr '\n' ' ' | cut -c1-900)"
if /usr/bin/grep -q "DONE. Next" $W/efr497/def_day.log; then
  FIN="$(/usr/bin/grep -o 'final: [^ ]* factors [^ ]*' $W/efr497/def_day.log | tail -1)"
  log "final rho: $(uv run python scripts/def_day_helpers.py rho "$W/efr497/$(echo "$FIN" | awk '{print $2}').ntclevel.json" 2>&1 | tail -1)"
  log "LODO cap on efr: $(python3 -c 'import json,sys; d=json.load(open(sys.argv[1])); print("hi", d["hi"], "hi_pct", d["hi_pct"], "frac_at_cap", d["frac_at_cap"])' $W/efr497/pert_scale_lodo_def_lo0.5_pct.json 2>&1)"
  log "meta: $(python3 -c 'import json,sys; c=json.load(open(sys.argv[1]))["config"]; print({k: c.get(k) for k in ("bulk_sharpen_q", "bulk_sharpen_pt", "bulk_ctx_scale")}, "gate", json.dumps(c.get("de_prior_gate"))[:200])' $W/efr497/497efr_final.h5ad.meta.json 2>&1 | cut -c1-700)"
  log "PASS: every October switch together completes on efr"
else
  log "FAIL: the all-switches day did not complete on efr (see $W/efr497/def_day.log)"
fi
rm -f $W/efr497/*.h5ad
[ "$(shasum -a 256 $BU | cut -d' ' -f1)" = "$(git show HEAD:$BU | shasum -a 256 | cut -d' ' -f1)" ] && log "builder restored" || log "WARNING builder differs from HEAD after the day"
log "END"
