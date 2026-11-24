#!/usr/bin/env bash
# 417 bulk DUMP runner: run the driver to the build stage with VCC_BULK_DUMP=<out.npz>; the builder writes the finished bulk
# matrix (after VCC_PERT_SCALE / VCC_BULK_GENE_GAIN) and exits before emission, so the driver's post-build assert fails BY
# DESIGN. The arm env (VCC_* knobs) is inherited from the caller. usage: run_dump417.sh <name> <out.npz> <builder_sha>
set -uo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"; REPO="$(cd "$HERE/../.." && pwd)"; cd "$REPO"
NAME="${1:?name}"; OUT="${2:?out.npz}"; WANT_B="${3:?builder sha}"
W="$HERE/work/dump_$NAME"; mkdir -p "$W"; LOG="$W/run_dump.log"
BUNDLE="$REPO/data/challenge"
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
  [ "${pct:-0}" -ge 40 ] && awk -v v="$fi" 'BEGIN{exit !(v>=23)}' && [ "$others" = "0" ] \
    && [ "$preps" = "0" ] && awk -v v="$disk" 'BEGIN{exit !(v>=20)}'
}
wait_gate(){ local n=0; until gate_ok; do n=$((n+1)); [ $n -ge 240 ] && { log "GATE never opened"; exit 3; }; sleep 30; done; }
restore(){ [ -f "$W/pin-artifacts.json" ] && uv run python scripts/def_pin_panel_artifacts.py --restore --record "$W/pin-artifacts.json" | tee -a "$LOG"; }

export VCC_BULK_DUMP="$OUT"
log "ARMED dump pid=$$ out=$OUT env: $(env | grep '^VCC_' | sort | tr '\n' ' ')"
B="$(shasum -a 256 scripts/make_submission_pm_cd4.py | cut -d' ' -f1)"
[ "$B" = "$WANT_B" ] || { log "BUILDER hash $B != $WANT_B -- stopping"; exit 4; }
[ -e "$OUT" ] && { log "$OUT exists -- stopping"; exit 4; }
wait_gate
/usr/bin/time -l uv run python scripts/def_pipeline.py \
    --bundle "$BUNDLE" --work "$W" --name "dump-$NAME" \
    --whitener-id union_pruned_def --emitter rc_callset \
    --frozen-inputs "$W/frozen_inputs.json" --min-free-gib 23 \
    --stages preflight,identity,cd4-extract,census,whitener,repin-frozen,pin-artifacts,repin-panel-support,build \
    > "$W/driver.log" 2>&1
RC=$?; log "driver rc=$RC (non-zero is expected: no h5ad after a dump)"
restore
log "builder after restore: $(shasum -a 256 scripts/make_submission_pm_cd4.py | cut -d' ' -f1)"
L="$(grep -h '\[bulk-dump\] wrote' "$W/logs/build.log" 2>/dev/null | tail -1)"
if [ -s "$OUT" ] && [ -n "$L" ]; then log "DUMP OK: $L"; log "dump sha256 $(shasum -a 256 "$OUT" | cut -d' ' -f1)"; else log "DUMP FAILED (no file or no log line) -- stopping"; exit 1; fi
log "DONE"
