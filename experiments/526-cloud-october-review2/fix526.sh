#!/usr/bin/env bash
# 526 applied (2026-09-29): the inert known answer for the three production edits taken from 526's review --
# def_day.sh refuses a fresh run whose builder is not HEAD (finding 3), def_pin_panel_artifacts.py writes its record and
# before-image BEFORE replacing the builder (finding 3), def_preflight.py fails closed on a control file without ntc_id
# when the manifest declares NTC guides (finding 2; DEF_DAY_ALLOW_NO_NTC_ID=1 overrides). Waits for build527's process to
# exit and the builder to equal HEAD (its committer), then runs efp relaxed through def_day.sh; the final artifact must be
# 474/476/497/502's ac49d319 byte for byte (this also proves 523's and 527's builder knobs inert on the October path).
# commit_after_chain.sh commits the three scripts (commit526.spec). Never submits.
set -uo pipefail
REPO=~/Code/virtual-cell-challenge-2026; cd "$REPO" || exit 1
X=experiments/526-cloud-october-review2; W=$X/work; LOG=$W/fix526.log
DD=scripts/def_day.sh; BU=scripts/make_submission_pm_cd4.py; EFP=experiments/238-essential-free-panel/work/ef-bundle
mkdir -p $W
log(){ printf '%s %s\n' "$(date -u +%FT%TZ)" "$*" >> "$LOG"; }
eq_head(){ [ "$(shasum -a 256 "$1" | cut -d' ' -f1)" = "$(git show HEAD:"$1" | shasum -a 256 | cut -d' ' -f1)" ]; }
fi_gib(){ vm_stat | awk '/Pages free/{f=$3}/Pages inactive/{i=$3}END{gsub(/\./,"",f);gsub(/\./,"",i);printf "%d",(f+i)*16384/2^30}'; }
busy(){ ps -axo comm,args | awk '$1 ~ /^bash$/ && /build527\.sh|run_ka_[0-9]|run_arm_generic\.sh|def_day\.sh|run_dump4/' | /usr/bin/grep -q . && return 0
        [ "$(uv run python -c 'import sys; sys.path.insert(0,"scripts"); import def_pipeline as P; print(len(P.other_builders()))' 2>/dev/null)" != 0 ]; }
echo $$ > $W/fix526.pid
log "armed (pid $$): waiting for build527 to exit, the builder == HEAD, no other build, free+inactive >= 30 GiB"
until ! busy && eq_head $BU && [ "$(fi_gib)" -ge 30 ]; do
  [ "$(date -u +%Y%m%d%H%M)" -ge 202610021200 ] && { log "prerequisites not met by 12:00Z 10-02 -- stopping"; log "END"; exit 1; }
  sleep 60; done
log "released: builder $(shasum -a 256 $BU | cut -c1-16) (HEAD), def_day $(shasum -a 256 $DD | cut -c1-16) (526-patched, uncommitted)"
bash -n $DD || { log "def_day.sh syntax"; log "END"; exit 1; }
rm -rf $W/efp_inert
DEF_DAY_CD4_THRESHOLDS=relaxed bash $DD $EFP $W/efp_inert 474efp 1.75 auto pt:noprog 4 > $W/efp_inert.nohup 2>&1
/usr/bin/grep -q "DONE. Next" $W/efp_inert/def_day.log || { log "KA526 FAIL: def_day did not finish: $(tail -2 $W/efp_inert/def_day.log | tr '\n' ' ' | cut -c1-300)"; log "END"; exit 1; }
GOT="$(uv run python -c 'import json,sys; print(json.load(open(sys.argv[1]))["output"]["sha256"])' $W/efp_inert/474efp_final.h5ad.meta.json)"
rm -f $W/efp_inert/*.h5ad
eq_head $BU || { log "KA526 FAIL: builder not HEAD after def_day (trap)"; log "END"; exit 1; }
[ "$GOT" = ac49d319f5982507003b1f2070660d2d03ee351fd13c9fdd30df8f5fb08ec0b3 ] || { log "KA526 FAIL: final ${GOT:0:16} != ac49d319"; log "END"; exit 1; }
log "KA526 PASS: efp relaxed through the 526-patched def_day.sh on builder $(shasum -a 256 $BU | cut -c1-16) retraces ac49d319 byte for byte"
log "END"
