#!/usr/bin/env bash
# 502 applied (2026-09-28, after 446bj2z read "zs TRANSFERS"): DEF_DAY_DE_MAG / DEF_DAY_MDE_SIGN pass-throughs in
# def_day.sh + helpers (patch_defday502.py on 497's e29e793a / 60b2a777 -> 610db782 / 20142142 == repro/*.new), then:
#   tests: bash -n; test_riders502.py 9/9; test_helpers497.py 9/9 on the 502 helpers; 485's repro_fixes on the 502 copies;
#   (a) inert KA: efp relaxed with the 502 script and neither rider == ac49d319 byte for byte (474/476/497's retrace);
#   (b) rider on: efp relaxed + DEF_DAY_DE_MAG=dn_k250_zs -> DONE, postbuild PASS with de_mag allowed, riders_match, the
#       final meta's de_mag.variant == dn_k250_zs, and the output differs from (a) (the rider acted).
# Serial after the whole queue (chain511's process exit, no other build) and after the 496/511 committers (builder,
# emitter, def_day and helpers all == HEAD). Any failure in apply/tests/(a) restores the pre-502 files. Never submits.
# Committed by commit_after_chain.sh (spec commit502.spec).
set -uo pipefail
REPO=~/Code/virtual-cell-challenge-2026; cd "$REPO" || exit 1
X=experiments/502-october-de-riders; W=$X/work; LOG=$W/chain502.log
DD=scripts/def_day.sh; HP=scripts/def_day_helpers.py; BU=scripts/make_submission_pm_cd4.py; EMP=src/vcc2026/rc_emit_k.py
EFP=experiments/238-essential-free-panel/work/ef-bundle
log(){ printf '%s %s\n' "$(date -u +%FT%TZ)" "$*" >> "$LOG"; }
eq_head(){ [ "$(shasum -a 256 "$1" | cut -d' ' -f1)" = "$(git show HEAD:"$1" | shasum -a 256 | cut -d' ' -f1)" ]; }
echo $$ > $W/chain502.pid
log "armed (pid $$): waiting for chain511's process to exit, no other build, and builder/emitter/def_day/helpers == HEAD"
# 519 finding 2: a kill between the patch's two file writes (or during a run) can leave one file new and one old; the
# release gate below requires the pre-502 sha, so such a tree would wait until 10-03 and never repair itself.
if [ -f $W/pre502.def_day.sh ] && [ -f $W/pre502.def_day_helpers.py ]; then
  DDSHA="$(shasum -a 256 $DD | cut -c1-16)"; HPSHA="$(shasum -a 256 $HP | cut -c1-16)"
  case "$DDSHA $HPSHA" in
    "e29e793a7964df03 60b2a7771ef857d4"|"610db78277239f50 20142142e4574b7f") ;;
    *) cp $W/pre502.def_day.sh $DD; cp $W/pre502.def_day_helpers.py $HP
       log "reverted a partial 502 pair to the pre-502 copies (was $DDSHA $HPSHA)";;
  esac
fi
P511="$(cat experiments/511-tercile-ties/work/chain511.pid 2>/dev/null)"
until ! ps -p "${P511:-1}" -o args= 2>/dev/null | /usr/bin/grep -q chain511 \
      && ! ps -axo comm,args | awk '$1 ~ /^bash$/ && /chain4(7[69]b?|9[5-8])\.sh|chain496p\.sh|chain497e\.sh|chain51[01]\.sh|run_ka_4|run_arm_generic\.sh|def_day\.sh|run_dump4/' | /usr/bin/grep -q . \
      && [ "$(shasum -a 256 $DD | cut -c1-16)" = e29e793a7964df03 ] && eq_head $DD && eq_head $HP && eq_head $BU && eq_head $EMP; do
  [ "$(date -u +%Y%m%d%H%M)" -ge 202610031200 ] && { log "prerequisites not met by 12:00Z 10-03 -- stopping"; log "END"; exit 1; }; sleep 60; done
log "released: builder $(shasum -a 256 $BU | cut -c1-16), emitter $(shasum -a 256 $EMP | cut -c1-16), def_day $(shasum -a 256 $DD | cut -c1-16)"
cp $DD $W/pre502.def_day.sh; cp $HP $W/pre502.def_day_helpers.py
back(){ cp $W/pre502.def_day.sh $DD; cp $W/pre502.def_day_helpers.py $HP; rm -f $W/efp_inert/*.h5ad $W/efp_zs/*.h5ad; log "RESTORED the two pre-502 files ($1)"; log "END"; exit 1; }
uv run python $X/patch_defday502.py $DD $HP >> $LOG 2>&1 || back "patch failed"
cmp -s $DD $X/repro/def_day.sh.new && cmp -s $HP $X/repro/def_day_helpers.py.new || back "applied files differ from repro/*.new"
bash -n $DD || back "def_day.sh syntax"
uv run python $X/test_riders502.py $HP > $W/test_riders502.log 2>&1 || back "test_riders502 failed"
uv run python experiments/497-october-switches/test_helpers497.py $HP > $W/test_helpers497.log 2>&1 || back "test_helpers497 on the 502 helpers failed"
(cd $X/repro && uv run python repro_fixes.py > last-fixes-chain.log 2>&1) || back "repro_fixes on the 502 files failed"
log "applied: def_day $(shasum -a 256 $DD | cut -c1-16), helpers $(shasum -a 256 $HP | cut -c1-16); test_riders502 9/9; test_helpers497 9/9; repro_fixes PASS"
memwait(){ local n=0 fi; while :; do
  fi="$(vm_stat | awk '/Pages free/{f=$3}/Pages inactive/{i=$3}END{gsub(/\./,"",f);gsub(/\./,"",i);printf "%.1f",(f+i)*16384/2^30}')"
  awk -v v="$fi" 'BEGIN{exit !(v>=30)}' && return 0
  n=$((n+1)); [ $n -ge 120 ] && return 1; log "waiting for memory: free+inactive $fi GiB < 30"; sleep 60; done; }
# --- (a) inert KA
memwait || back "memory never reached 30 GiB before (a)"
rm -rf $W/efp_inert
DEF_DAY_CD4_THRESHOLDS=relaxed bash $DD $EFP $W/efp_inert 474efp 1.75 auto pt:noprog 4 > $W/efp_inert.nohup 2>&1
/usr/bin/grep -q "DONE. Next" $W/efp_inert/def_day.log || back "(a) did not finish: $(tail -2 $W/efp_inert/def_day.log | tr '\n' ' ')"
GOTA="$(uv run python -c 'import json,sys; print(json.load(open(sys.argv[1]))["output"]["sha256"])' $W/efp_inert/474efp_final.h5ad.meta.json)"
[ "$GOTA" = ac49d319f5982507003b1f2070660d2d03ee351fd13c9fdd30df8f5fb08ec0b3 ] || back "(a) final ${GOTA:0:16} != ac49d319"
eq_head $BU || back "builder not restored after (a)"
log "(a) PASS: the 502 script with neither rider retraces efp relaxed byte for byte (${GOTA:0:16})"
# --- (b) rider on. A failure restores the pre-502 files (519 finding 2): an applied-but-unproven pass-through must not
#     stay in the tree uncommitted; the October resolver then refuses a recipe with de_mag (the knob is missing).
memwait || back "memory never reached 30 GiB before (b)"
rm -rf $W/efp_zs
DEF_DAY_DE_MAG=dn_k250_zs DEF_DAY_CD4_THRESHOLDS=relaxed bash $DD $EFP $W/efp_zs 502zs 1.75 auto pt:noprog 4 > $W/efp_zs.nohup 2>&1
if /usr/bin/grep -q "DONE. Next" $W/efp_zs/def_day.log \
   && /usr/bin/grep -q "postbuild PASS" $W/efp_zs/def_day.log \
   && /usr/bin/grep -q "DE rider for the level builds: VCC_DE_MAG=dn_k250_zs" $W/efp_zs/def_day.log; then
  V="$(uv run python -c 'import json,sys; c=json.load(open(sys.argv[1]))["config"]; print((c.get("de_mag") or {}).get("variant"))' $W/efp_zs/502zs_final.h5ad.meta.json)"
  GOTB="$(uv run python -c 'import json,sys; print(json.load(open(sys.argv[1]))["output"]["sha256"])' $W/efp_zs/502zs_final.h5ad.meta.json)"
  if [ "$V" = dn_k250_zs ] && [ "$GOTB" != "$GOTA" ]; then
    log "(b) PASS: rider on -> DONE, postbuild PASS, riders_match, meta de_mag $V, output ${GOTB:0:16} != (a)"
  else log "(b) FAIL: meta de_mag '$V', output ${GOTB:0:16} (a ${GOTA:0:16})"; back "(b) failed"; fi
else log "(b) FAIL: $(/usr/bin/grep -E 'STOP|postbuild|DE rider|DONE' $W/efp_zs/def_day.log | tail -3 | tr '\n' ' ')"; back "(b) failed"; fi
eq_head $BU || log "WARNING builder differs from HEAD after (b)"
rm -f $W/efp_inert/*.h5ad $W/efp_zs/*.h5ad
log "END def_day $(shasum -a 256 $DD | cut -c1-16), helpers $(shasum -a 256 $HP | cut -c1-16) (committed by commit_after_chain.sh)"
