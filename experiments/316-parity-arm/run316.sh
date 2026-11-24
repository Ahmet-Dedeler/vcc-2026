#!/usr/bin/env bash
# 316 -- BUILD the PARITY-CORRECT normrel arm on the REAL A/B/C challenge panel.
# Modelled byte-for-byte on run312.sh's gated two-invocation pattern, with one difference: the
# postbuild config-drift FAIL on the deliberate arm key is EXPECTED (312 saw exactly one, 33/34),
# so invocation 1's rc=1 is inspected rather than treated as a build failure. Nothing is submitted
# here and `vcc submit` is never run.
set -uo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"; REPO="$(cd "$HERE/../.." && pwd)"; cd "$REPO"
W="$HERE/work/build"; mkdir -p "$W"; LOG="$W/run316.log"
BUNDLE="$REPO/data/challenge"
NAME="a-duck-reviews-the-diff"
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
restore(){ [ -f "$W/pin-artifacts.json" ] && uv run python scripts/def_pin_panel_artifacts.py --restore --record "$W/pin-artifacts.json" | tee -a "$LOG"; }

export VCC_NORMREL_EXP=-1.5          # 312's dose, unchanged -- the arm here is the RENDERING
export VCC_NORMREL_PM=blend          # THE ARM (313 section 8): PerturbMatch inside the norm match
log "ARMED pid=$$  VCC_NORMREL_EXP=$VCC_NORMREL_EXP  VCC_NORMREL_PM=$VCC_NORMREL_PM"
log "builder: $(shasum -a 256 scripts/make_submission_pm_cd4.py | cut -d' ' -f1)  driver: $(shasum -a 256 scripts/def_pipeline.py | cut -d' ' -f1)  git: $(git rev-parse --short HEAD)"
wait_gate
log "DRIVER invocation 1: preflight ... postbuild  (panel A/B/C, real challenge bundle)"
/usr/bin/time -l uv run python scripts/def_pipeline.py \
    --bundle "$BUNDLE" --work "$W" --name "$NAME" \
    --whitener-id union_pruned_def --emitter rc_callset \
    --frozen-inputs "$W/frozen_inputs.json" --min-free-gib 23 \
    --stages preflight,identity,cd4-extract,census,whitener,repin-frozen,pin-artifacts,repin-panel-support,build,postbuild \
    > "$W/driver.log" 2>&1
RC=$?; log "invocation 1 rc=$RC  builder now: $(shasum -a 256 scripts/make_submission_pm_cd4.py | cut -d' ' -f1)"
log 'restoring builder (the pin is finished; the prep does not need it)'
restore
PB="$W/logs/postbuild.log"
NFAIL="$(grep -c '^  FAIL' "$PB" 2>/dev/null || echo 99)"
log "postbuild: $(grep -hE 'checks pass' "$PB" 2>/dev/null | tail -1)  FAIL lines: $NFAIL"
grep -h '^  FAIL' "$PB" 2>/dev/null | tee -a "$LOG"
if [ ! -s "$W/$NAME.h5ad" ]; then log "NO ARTIFACT -- stopping"; exit 1; fi
if [ "$NFAIL" != "1" ]; then log "UNEXPECTED postbuild failure count ($NFAIL != 1) -- stopping before prep"; exit 1; fi
if ! grep -q '^  FAIL  config differences' "$PB"; then log "the single FAIL is NOT the config-drift line -- stopping"; exit 1; fi
META="$W/$NAME.h5ad.meta.json"
log "meta normrel_exp=$(uv run python -c 'import json,sys; print(json.load(open(sys.argv[1]))["config"].get("normrel_exp","<ABSENT -- ARM DID NOT LAND>"))' "$META")"
log "meta normrel_pm=$(uv run python -c 'import json,sys; print(json.load(open(sys.argv[1]))["config"].get("normrel_pm","<ABSENT -- ARM DID NOT LAND>"))' "$META")"
log "artifact: $(ls -la "$W/$NAME.h5ad" | awk '{print $5, $9}')  sha256 $(uv run python -c 'import json,sys; print(json.load(open(sys.argv[1]))["output"]["sha256"])' "$META")"

wait_gate
log "INVOCATION 2: vcc prep -> submissions/316-normrel15-pm.vcc  (the driver's prep stage is --dry-run,"
log "  i.e. validation only, so the package is written by this real prep; ~22-24 GiB, gated above)"
/usr/bin/time -l vcc prep "$W/$NAME.h5ad" -g "$BUNDLE/gene_names.csv" \
    --perts "$BUNDLE/pert_counts.csv" --contexts A,B,C \
    -o "$REPO/submissions/316-normrel15-pm.vcc" > "$W/prep.log" 2>&1
RC=$?; log "prep rc=$RC"
tail -20 "$W/prep.log" | tee -a "$LOG"
log "vcc: $(ls -la "$REPO/submissions/316-normrel15-pm.vcc" 2>/dev/null | awk '{print $5, $9}')"
log "git status scripts/: $(git status --short scripts/ | tr '\n' ' ')"
log "DONE"
