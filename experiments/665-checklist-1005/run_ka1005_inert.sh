#!/usr/bin/env bash
# 665 (coordinator, 2026-10-05): 10-05 CHECKLIST v2 step 7, the inert known-answer build on the post-bundle tree (538 + 554 +
# 590 + 641 applied). 550 RESULTS "10-05 command sequence" step 4, plus the five knobs 550's list misses. The output must be
# ac49d319... (526/efp_inert/474efp_final.h5ad.meta.json); postbuild.json must carry dose_guard per context; def_day.log must
# have "DONE. Next". NEVER chain502 / fix526 / chain476.
# Memory: two iOS Simulators keep free+inactive at 21-28 GiB and def_pipeline refuses below 23 at each build stage, so this
# waits for >= 27 GiB on two readings 60 s apart (as 664 did for 524dr) and retries only on that refusal. Cap 10-05 16:00Z.
# A finished run is judged by the checks below and ends with "KA1005 PASS" / "KA1005 FAIL <why>".
set -uo pipefail
REPO=~/Code/virtual-cell-challenge-2026; cd "$REPO" || exit 1
X=experiments/665-checklist-1005; mkdir -p $X/work; LOG=$X/work/ka1005.log
KA=experiments/502-october-de-riders/work/ka1005-inert; WANT=ac49d319f5982507003b1f2070660d2d03ee351fd13c9fdd30df8f5fb08ec0b3
CAP=202610051600; FLOOR="${KA_FLOOR:-27}"   # release floor (GiB); 10-05 09:4xZ re-armed at 25.5 once the simulators were gone
log(){ printf '%s %s\n' "$(date -u +%FT%TZ)" "$*" >> "$LOG"; }
fi_gib(){ vm_stat | awk '/Pages free/{f=$3}/Pages inactive/{i=$3}END{gsub(/\./,"",f);gsub(/\./,"",i);printf "%.1f",(f+i)*16384/2^30}'; }
busy(){ uv run python -c 'import sys; sys.path.insert(0,"scripts"); import def_pipeline as P; raise SystemExit(0 if P.other_builders() else 1)' \
        || ps -axo args= | grep -vE 'claude|grep' | grep -qE '(^|/| )vcc prep( |$)'; }
echo $$ > $X/work/ka1005.pid
log "armed (pid $$): release at free+inactive >= $FLOOR GiB twice 60 s apart; cap $CAP"
N=0
while [ "$(date -u +%Y%m%d%H%M)" -lt $CAP ]; do
  a="$(fi_gib)"
  if awk -v v="$a" -v f="$FLOOR" 'BEGIN{exit !(v>=f)}' && ! busy; then
    sleep 60; b="$(fi_gib)"
    if awk -v v="$b" -v f="$FLOOR" 'BEGIN{exit !(v>=f)}' && ! busy; then
      N=$((N + 1))
      [ -e "$KA" ] && mv "$KA" "$KA.aborted-$(date -u +%Y%m%dT%H%M%SZ)"
      log "attempt $N: free+inactive $a then $b GiB -> def_day (inert, efp relaxed) into $KA"
      env -u DEF_DAY_LODO -u DEF_DAY_LODO_HI -u DEF_DAY_LODO_HI_PCT -u DEF_DAY_LODO_STRENGTH \
          -u DEF_DAY_DE_MAG -u DEF_DAY_MDE_SIGN -u DEF_DAY_CALLSET_BLOCK -u DEF_DAY_CENSUS_MAX_ZERO \
          -u DEF_DAY_RHO_SCALE -u DEF_DAY_RC_CALLSET_CPM -u DEF_DAY_CD4_STACK -u DEF_DAY_HOOK_ESTIMATOR \
          -u DEF_DAY_PROGRAMME_SOURCE -u DEF_DAY_PRIOR_GATE -u DEF_DAY_SHARPEN_PT -u DEF_DAY_FORCE \
          -u DEF_DAY_KEEP_PINS -u DEF_DAY_ALLOW_NO_NTC_ID \
          -u DEF_DAY_ALLOW_CELLS_PER_PERT -u DEF_DAY_ALLOW_DEPTH -u DEF_DAY_ALLOW_LOW_DOSE -u DEF_DAY_ALLOW_WIDE_LIBS \
          -u DEF_DAY_INCUMBENT_META \
          DEF_DAY_CD4_THRESHOLDS=relaxed \
          bash scripts/def_day.sh experiments/238-essential-free-panel/work/ef-bundle "$KA" 474efp 1.75 auto pt:noprog 4 \
          > $X/work/attempt$N.nohup 2>&1
      RC=$?
      log "attempt $N: def_day rc=$RC; def_day.log tail: $(tail -2 $KA/def_day.log 2>/dev/null | tr '\n' ' ' | cut -c1-300)"
      if /usr/bin/grep -q "DONE. Next" $KA/def_day.log 2>/dev/null; then
        GOT="$(python3 -c 'import json,sys; print(json.load(open(sys.argv[1]))["output"]["sha256"])' $KA/474efp_final.h5ad.meta.json 2>/dev/null)"
        # def_postbuild.dose_guard: {"floors", "source", "per_context": {ctx: {...}, "pooled": {...}}} (rows dropped)
        DG="$(python3 -c 'import json,sys; g=json.load(open(sys.argv[1])).get("dose_guard") or {}; pc=[k for k in (g.get("per_context") or {}) if k != "pooled"]; print("per-context" if len(pc) >= 3 else "MISSING per_context " + ",".join(pc))' $KA/postbuild.json 2>/dev/null)"
        log "attempt $N: output.sha256 ${GOT:-<none>} (want ${WANT:0:16}); postbuild dose_guard: ${DG:-<unreadable>}"
        if [ "$GOT" = "$WANT" ] && [ "$DG" = "per-context" ]; then log "KA1005 PASS"; log "END"; exit 0; fi
        log "KA1005 FAIL output or dose_guard mismatch"; log "END"; exit 1
      fi
      if /usr/bin/grep -qh "refusing" $KA/driver-*.log 2>/dev/null; then
        log "attempt $N: def_pipeline's memory refusal -- retrying after 5 min"; sleep 300; continue
      fi
      log "KA1005 FAIL def_day did not finish for another reason (see $KA)"; log "END"; exit 1
    fi
  fi
  sleep 60
done
log "KA1005 FAIL cap $CAP reached ($N attempts)"; log "END"; exit 1
