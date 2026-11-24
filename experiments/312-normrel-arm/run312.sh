#!/usr/bin/env bash
# 312 -- BUILD the normrel arm on the REAL A/B/C challenge panel, for a board slot.
# Modelled on rehearsal311.sh: gated, two invocations, builder restored, nothing submitted here.
set -uo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"; REPO="$(cd "$HERE/../.." && pwd)"; cd "$REPO"
W="$HERE/work/build"; mkdir -p "$W"; LOG="$W/run312.log"
BUNDLE="$REPO/data/challenge"
NAME="the-cat-sat-on-the-modem"
log(){ printf '%s %s\n' "$(date -u +%FT%TZ)" "$*" | tee -a "$LOG"; }
gate_ok(){
  local pct fi others preps
  pct="$(memory_pressure | awk '/free percentage/{gsub(/%/,"",$5); print $5}')"
  fi="$(vm_stat | awk '/Pages free/{f=$3}/Pages inactive/{i=$3}END{gsub(/\./,"",f);gsub(/\./,"",i);printf "%.1f",(f+i)*16384/2^30}')"
  others="$(uv run python -c 'import sys; sys.path.insert(0,"scripts"); import def_pipeline as P; print(len(P.other_builders()))' 2>/dev/null || echo 99)"
  preps="$(ps -axo args= | grep -vE 'claude|grep' | grep -cE '(^|/| )vcc prep( |$)' | tr -d ' ')"
  log "gate: free ${pct}%  free+inactive ${fi} GiB  other builders ${others}  vcc preps ${preps}"
  [ "${pct:-0}" -ge 40 ] && awk -v v="$fi" 'BEGIN{exit !(v>=23)}' && [ "$others" = "0" ] && [ "$preps" = "0" ]
}
wait_gate(){ local n=0; until gate_ok; do n=$((n+1)); [ $n -ge 240 ] && { log "GATE never opened"; exit 3; }; sleep 30; done; }

export VCC_NORMREL_EXP=-1.5          # THE ARM (097 a6 / 308 normrel1.5), shape-only
log "ARMED pid=$$  VCC_NORMREL_EXP=$VCC_NORMREL_EXP  builder: $(shasum -a 256 scripts/make_submission_pm_cd4.py | cut -d' ' -f1)  git: $(git rev-parse --short HEAD)"
wait_gate
log "DRIVER invocation 1: preflight ... postbuild  (panel A/B/C, real challenge bundle)"
/usr/bin/time -l uv run python scripts/def_pipeline.py \
    --bundle "$BUNDLE" --work "$W" --name "$NAME" \
    --whitener-id union_pruned_def --emitter rc_callset \
    --frozen-inputs "$W/frozen_inputs.json" --min-free-gib 23 \
    --stages preflight,identity,cd4-extract,census,whitener,repin-frozen,pin-artifacts,repin-panel-support,build,postbuild \
    > "$W/driver.log" 2>&1
RC=$?; log "invocation 1 rc=$RC  builder now: $(shasum -a 256 scripts/make_submission_pm_cd4.py | cut -d' ' -f1)"
if [ $RC -ne 0 ]; then
  log "FAILED; restoring builder"; [ -f "$W/pin-artifacts.json" ] && uv run python scripts/def_pin_panel_artifacts.py --restore --record "$W/pin-artifacts.json" | tee -a "$LOG"; exit 1
fi
wait_gate
log "DRIVER invocation 2: prep + restore-builder"
/usr/bin/time -l uv run python scripts/def_pipeline.py \
    --bundle "$BUNDLE" --work "$W" --name "$NAME" \
    --whitener-id union_pruned_def --emitter rc_callset \
    --frozen-inputs "$W/frozen_inputs.json" --min-free-gib 23 \
    --stages prep --restore-builder-after \
    > "$W/driver-prep.log" 2>&1
RC=$?; log "invocation 2 rc=$RC  builder now: $(shasum -a 256 scripts/make_submission_pm_cd4.py | cut -d' ' -f1)"
[ $RC -ne 0 ] && { log "prep FAILED; restoring"; uv run python scripts/def_pin_panel_artifacts.py --restore --record "$W/pin-artifacts.json" | tee -a "$LOG"; exit 1; }
META="$W/$NAME.h5ad.meta.json"
log "meta normrel_exp=$(uv run python -c 'import json,sys; print(json.load(open(sys.argv[1]))["config"].get("normrel_exp","<ABSENT -- ARM DID NOT LAND>"))' "$META")"
log "postbuild: $(grep -hE 'checks pass' "$W"/logs/postbuild.log 2>/dev/null | tail -1)  FAIL: $(grep -c '^  FAIL' "$W"/logs/postbuild.log 2>/dev/null)"
log "artifact: $(ls -la "$W/$NAME.h5ad" 2>/dev/null | awk '{print $5, $9}')"
log "git status scripts/: $(git status --short scripts/ | tr '\n' ' ')"
log "DONE"
