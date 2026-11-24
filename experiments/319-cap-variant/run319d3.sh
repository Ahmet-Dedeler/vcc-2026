#!/usr/bin/env bash
# 319 -- BUILD the SUPPRESSION-ONLY normrel arm (`VCC_NORMREL_MODE=cap`) on the REAL A/B/C panel.
# Modelled on run316.sh, which is run312.sh's gated two-invocation pattern.  ONE variable against
# entry Q519I2pLWRCW7PAEx5ul (`312-normrel15`): the same exponent -1.5, `VCC_NORMREL_PM` UNSET, only
# the ratio clamp added.  Nothing is submitted here and `vcc submit` is never run.
set -uo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"; REPO="$(cd "$HERE/../.." && pwd)"; cd "$REPO"
W="$HERE/work/build_d3"; mkdir -p "$W"; LOG="$W/run319d3.log"
BUNDLE="$REPO/data/challenge"
NAME="a-moth-in-the-lint-trap"
PKG="$REPO/submissions/319-normrel15-cap.vcc"
log(){ printf '%s %s\n' "$(date -u +%FT%TZ)" "$*" | tee -a "$LOG"; }
free_gib(){ df -k "$REPO" | awk 'NR==2{printf "%.1f", $4/1048576}'; }
gate_ok(){
  local pct fi others preps disk
  pct="$(memory_pressure | awk '/free percentage/{gsub(/%/,"",$5); print $5}')"
  fi="$(vm_stat | awk '/Pages free/{f=$3}/Pages inactive/{i=$3}END{gsub(/\./,"",f);gsub(/\./,"",i);printf "%.1f",(f+i)*16384/2^30}')"
  others="$(uv run python -c 'import sys; sys.path.insert(0,"scripts"); import def_pipeline as P; print(len(P.other_builders()))' 2>/dev/null || echo 99)"
  preps="$(ps -axo args= | grep -vE 'claude|grep' | grep -cE '(^|/| )vcc prep( |$)' | tr -d ' ')"
  disk="$(free_gib)"
  log "gate: free ${pct}%  free+inactive ${fi} GiB  other builders ${others}  vcc preps ${preps}  disk ${disk} GiB"
  # the build writes 13.2 GiB (h5ad 11.8 + whitener 1.2 + cd4_compact 0.2) and the prep 2.6 GiB
  # beside it, measured on 316; 20 GiB free is that with ~4 GiB of slack.
  [ "${pct:-0}" -ge 40 ] && awk -v v="$fi" 'BEGIN{exit !(v>=23)}' && [ "$others" = "0" ] \
    && [ "$preps" = "0" ] && awk -v v="$disk" 'BEGIN{exit !(v>=20)}'
}
wait_gate(){ local n=0; until gate_ok; do n=$((n+1)); [ $n -ge 240 ] && { log "GATE never opened"; exit 3; }; sleep 30; done; }
restore(){ [ -f "$W/pin-artifacts.json" ] && uv run python scripts/def_pin_panel_artifacts.py --restore --record "$W/pin-artifacts.json" | tee -a "$LOG"; }

export VCC_NORMREL_EXP=-3.0          # DOSE 2 on the cap form: where `plain` turned over
export VCC_NORMREL_MODE=cap          # THE ARM: ratio clamped at 1, so only OVERSIZED rows are suppressed
unset VCC_NORMREL_PM                 # the shipped weight-first rendering, i.e. 312's, not 316's
log "ARMED pid=$$  VCC_NORMREL_EXP=$VCC_NORMREL_EXP  VCC_NORMREL_MODE=$VCC_NORMREL_MODE  VCC_NORMREL_PM=<unset>"
log "builder: $(shasum -a 256 scripts/make_submission_pm_cd4.py | cut -d' ' -f1)  driver: $(shasum -a 256 scripts/def_pipeline.py | cut -d' ' -f1)  git: $(git rev-parse --short HEAD)"
wait_gate
log "DRIVER invocation 1: preflight ... postbuild  (panel A/B/C, real challenge bundle)"
/usr/bin/time -l uv run python scripts/def_pipeline.py \
    --bundle "$BUNDLE" --work "$W" --name "$NAME" \
    --whitener-id union_pruned_def --emitter rc_callset \
    --frozen-inputs "$W/frozen_inputs.json" --min-free-gib 23 \
    --stages preflight,identity,cd4-extract,census,whitener,repin-frozen,pin-artifacts,repin-panel-support,build,postbuild \
    > "$W/driver.log" 2>&1
RC=$?; log "invocation 1 rc=$RC  builder now: $(shasum -a 256 scripts/make_submission_pm_cd4.py | cut -d' ' -f1)  disk $(free_gib) GiB"
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
log "meta normrel_mode=$(uv run python -c 'import json,sys; print(json.load(open(sys.argv[1]))["config"].get("normrel_mode","<ABSENT -- ARM DID NOT LAND>"))' "$META")"
log "meta normrel_pm=$(uv run python -c 'import json,sys; print(json.load(open(sys.argv[1]))["config"].get("normrel_pm","<absent, as intended>"))' "$META")"
log "artifact: $(ls -la "$W/$NAME.h5ad" | awk '{print $5, $9}')  sha256 $(uv run python -c 'import json,sys; print(json.load(open(sys.argv[1]))["output"]["sha256"])' "$META")"

# the one-key config diff, BEFORE the prep: if it is not exactly one key this is not the arm.
uv run python "$HERE/configdiff319.py" "$META" | tee -a "$LOG"
RC=$?; [ $RC -ne 0 ] && { log "CONFIG DIFF is not exactly one key -- stopping before prep"; exit 1; }

wait_gate
log "INVOCATION 2: vcc prep -> $PKG  (the driver's prep stage is --dry-run, i.e. validation only,"
log "  so the package is written by this real prep; ~22-24 GiB RSS, gated above)"
/usr/bin/time -l vcc prep "$W/$NAME.h5ad" -g "$BUNDLE/gene_names.csv" \
    --perts "$BUNDLE/pert_counts.csv" --contexts A,B,C \
    -o "$PKG" > "$W/prep.log" 2>&1
RC=$?; log "prep rc=$RC  disk $(free_gib) GiB"
tail -20 "$W/prep.log" | tee -a "$LOG"
log "vcc: $(ls -la "$PKG" 2>/dev/null | awk '{print $5, $9}')"
cp "$META" "$REPO/submissions/319-normrel15-cap.h5ad.meta.json"
log "meta copied beside the package: $(shasum -a 256 "$REPO/submissions/319-normrel15-cap.h5ad.meta.json" | cut -d' ' -f1)"
uv run python scripts/check_submission_novelty.py check "$PKG" 2>&1 | tee -a "$LOG"
log "novelty rc=${PIPESTATUS[0]}"
log "git status scripts/: $(git status --short scripts/ | tr '\n' ' ')"
log "DONE"
