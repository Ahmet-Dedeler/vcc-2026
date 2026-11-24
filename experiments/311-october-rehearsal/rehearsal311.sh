#!/usr/bin/env bash
# 256 -- D/E/F REHEARSAL WITH NOTHING EXPORTED: the Oct-22 driver on 238's ef-bundle on its own
# DEFAULTS (256 wired the pinned basis into DEF_ENV), through `postbuild` and `vcc prep --dry-run`.
# Adapted from 251's rehearsal251.sh: the ONLY intended differences are (1) no export at all,
# (2) the work dir, (3) the gate at 23 GiB, (4) a final block that reads the pinned block off the meta.
# 244's invocation (work/run1) with: 244's fitted operator handed to --whitener-out (the driver's
# skip branch must ACCEPT it -- 247 B#18's axis proof; 244 section 5 observed rc 0 for exactly this
# file), a throwaway frozen table, --restore-builder-after (the pin is a rehearsal edit), and the
# mandate's 30 GiB gate.  Two invocations: everything through `postbuild`, then `prep` (+ restore)
# under this script's own gate, so a refused prep cannot leave the builder pinned (244 section 7.2).
set -uo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"; REPO="$(cd "$HERE/../.." && pwd)"; cd "$REPO"
W="$HERE/work/rehearsal"; mkdir -p "$W"; LOG="$W/rehearsal311.log"
BUNDLE="$REPO/experiments/238-essential-free-panel/work/ef-bundle"
OP244="$REPO/experiments/244-october-rehearsal/work/run1/whitener_union_pruned_def.npz"
NAME="lamp-on-the-landing"
log(){ printf '%s %s\n' "$(date -u +%FT%TZ)" "$*" | tee -a "$LOG"; }
gate_ok(){
  local pct fi others
  pct="$(memory_pressure | awk '/free percentage/{gsub(/%/,"",$5); print $5}')"
  fi="$(vm_stat | awk '/Pages free/{f=$3}/Pages inactive/{i=$3}END{gsub(/\./,"",f);gsub(/\./,"",i);printf "%.1f",(f+i)*16384/2^30}')"
  others="$(uv run python -c 'import sys; sys.path.insert(0,"scripts"); import def_pipeline as P; print(len(P.other_builders()))' 2>/dev/null || echo 99)"
  preps="$(ps -axo args= | grep -vE 'claude|grep' | grep -cE '(^|/| )vcc prep( |$)' | tr -d ' ')"
  log "gate: memory_pressure free ${pct}%  free+inactive ${fi} GiB  other builders ${others}  vcc preps ${preps}"
  [ "${pct:-0}" -ge 40 ] && awk -v v="$fi" 'BEGIN{exit !(v>=23)}' && [ "$others" = "0" ] && [ "$preps" = "0" ]
}
wait_gate(){ local n=0; until gate_ok; do n=$((n+1)); [ $n -ge 240 ] && { log "GATE never opened; giving up"; exit 3; }; sleep 30; done; }

log "ARMED pid=$$  builder before: $(shasum -a 256 scripts/make_submission_pm_cd4.py | cut -d' ' -f1)  git: $(git rev-parse --short HEAD) $(git status --short scripts/ | tr '\n' ' ')"
unset VCC_CD4_DEFLATE_BASIS VCC_CD4_DEFLATE_K VCC_UNION_SOURCES VCC_UNION2_SOURCES VCC_UNION2_G VCC_CD4_SIZING VCC_CD4_TARGET_G VCC_ALPHA VCC_CD4_TARGET_F VCC_H1_UNION_G
log "env: NOTHING exported (every VCC_* unset); the driver's defaults must supply every arm INCLUDING the pinned basis"
wait_gate
log "DRIVER invocation 1: preflight ... postbuild"
/usr/bin/time -l uv run python scripts/def_pipeline.py \
    --bundle "$BUNDLE" --work "$W" --name "$NAME" \
    --whitener-id union_pruned_def --whitener-out "$OP244" --emitter rc_callset \
    --frozen-inputs "$W/frozen_inputs.json" --cd4-thresholds relaxed --min-free-gib 23 \
    --stages preflight,identity,cd4-extract,census,whitener,repin-frozen,pin-artifacts,repin-panel-support,build,postbuild \
    > "$W/driver.log" 2>&1
RC=$?; log "DRIVER invocation 1 rc=$RC  peak RSS $(awk '/maximum resident set size/{printf "%.2f", $1/2^30}' "$W/driver.log" | tail -1) GiB  builder now: $(shasum -a 256 scripts/make_submission_pm_cd4.py | cut -d' ' -f1)"
if [ $RC -ne 0 ]; then
  log "invocation 1 FAILED; restoring the builder from the pin record if it was pinned"
  [ -f "$W/pin-artifacts.json" ] && uv run python scripts/def_pin_panel_artifacts.py --restore --record "$W/pin-artifacts.json" | tee -a "$LOG"
  exit 1
fi
wait_gate
log "DRIVER invocation 2: prep (--dry-run) + restore-builder"
/usr/bin/time -l uv run python scripts/def_pipeline.py \
    --bundle "$BUNDLE" --work "$W" --name "$NAME" \
    --whitener-id union_pruned_def --whitener-out "$OP244" --emitter rc_callset \
    --frozen-inputs "$W/frozen_inputs.json" --cd4-thresholds relaxed --min-free-gib 23 \
    --stages prep --restore-builder-after \
    > "$W/driver-prep.log" 2>&1
RC=$?; log "DRIVER invocation 2 rc=$RC  peak RSS $(awk '/maximum resident set size/{printf "%.2f", $1/2^30}' "$W/driver-prep.log" | tail -1) GiB  builder now: $(shasum -a 256 scripts/make_submission_pm_cd4.py | cut -d' ' -f1)"
if [ $RC -ne 0 ]; then
  log "invocation 2 FAILED; restoring the builder from the pin record"
  uv run python scripts/def_pin_panel_artifacts.py --restore --record "$W/pin-artifacts.json" | tee -a "$LOG"
  exit 1
fi
log "builder after restore: $(shasum -a 256 scripts/make_submission_pm_cd4.py | cut -d' ' -f1)  (intended $(cat "$HERE/work/builder.intended.sha256"))  git: $(git status --short scripts/ | tr '\n' ' ')"
# --- 256: read the pinned block off the finished meta, and the driver's own [arms] verdict
META="$W/$NAME.h5ad.meta.json"
log "meta config.cd4_deflate_basis=$(uv run python -c 'import json,sys; c=json.load(open(sys.argv[1]))["config"]; print(c.get("cd4_deflate_basis"), "k=", c.get("cd4_deflate_k","<absent>"))' "$META")  deflation.mode=$(uv run python -c 'import json,sys; d=json.load(open(sys.argv[1]))["effects"]["cd4"].get("deflation") or {}; print(d.get("mode"), d.get("basis"), "captured", d.get("captured_frac_of_centred"))' "$META")"
grep -E '^  (OK|NO ) config\.cd4_deflate' "$W/driver.log" | tee -a "$LOG" >/dev/null
log "postbuild: $(grep -hE 'checks pass' "$W"/logs/postbuild.log 2>/dev/null | tail -1)  FAIL lines: $(grep -c '^  FAIL' "$W"/logs/postbuild.log 2>/dev/null)"
log "DONE"
