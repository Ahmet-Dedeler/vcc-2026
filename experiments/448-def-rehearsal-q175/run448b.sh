#!/usr/bin/env bash
# 448b — third re-level pass on ef-bundle-r (the secant step from relevel_iter.py): F only moves, D/E stay at their
# converged lvl2 factors. Tests whether one bracketing secant step lands F inside 1 % of rho*. Waits for chain449 (446br2)
# to end so heavy builds stay serial; re-runs the driver's pin stages (the whitener/key fits are reused), builds, reads
# ntclevel, restores the builder.
set -uo pipefail
REPO=~/Code/virtual-cell-challenge-2026; cd "$REPO" || exit 1
X=experiments/448-def-rehearsal-q175; W=$X/work/efr; LOG=$X/work/run448b.log
B=experiments/238-essential-free-panel/work/ef-bundle-r
log(){ printf '%s %s\n' "$(date -u +%FT%TZ)" "$*" | tee -a "$LOG"; }
until grep -q " END$" experiments/449-repair-floor/work/chain449.log 2>/dev/null && ! ps -axo comm,args | awk '$1 ~ /^bash$/ && /chain449.sh/' | grep -q .; do
  [ "$(date -u +%Y%m%d%H%M)" -ge 202609251400 ] && { log "chain449 never ended by 14:00Z"; exit 1; }; sleep 30; done
for v in $(env | grep -o '^VCC_[A-Z0-9_]*'); do unset "$v"; done
BSHA0="$(shasum -a 256 scripts/make_submission_pm_cd4.py | cut -d' ' -f1)"
log "start: builder ${BSHA0:0:16}"
[ "$BSHA0" = "b6c3f13141e707012588cc354f62758889db9c38c3757f9765ebea37a8a733f8" ] || { log "builder is not b6c3f131 -- stopping"; exit 1; }
DRV=(uv run python scripts/def_pipeline.py --bundle "$B" --work "$W" --whitener-id union_pruned_def --emitter rc_callset
     --frozen-inputs "$W/frozen_inputs.json" --min-free-gib 23)
restore(){ [ -f "$W/pin-artifacts.json" ] && uv run python scripts/def_pin_panel_artifacts.py --restore --record "$W/pin-artifacts.json" >> "$LOG" 2>&1
           log "builder after restore: $(shasum -a 256 scripts/make_submission_pm_cd4.py | cut -c1-16) (want b6c3f1314)"; }
trap restore EXIT
"${DRV[@]}" --name 448b-setup --stages preflight,identity,cd4-extract,census,whitener,key,repin-frozen,pin-artifacts,repin-panel-support > "$W/driver-setup-b.log" 2>&1
log "setup rc=$?: $(grep -E '^--- \[' "$W/driver-setup-b.log" | tr '\n' ' ' | cut -c1-600)"
grep -q "^--- \[repin-panel-support\].*rc=0" "$W/driver-setup-b.log" || { log "setup did not complete -- stopping"; exit 1; }
DUMP="$W/dump_h1x2_cd4x1.npz"; DSHA="$(shasum -a 256 "$DUMP" | cut -d' ' -f1)"
[ "${DSHA:0:16}" = "16cdf3d3e3d2e8f8" ] || { log "dump sha changed (${DSHA:0:16}) -- stopping"; exit 1; }
F3="$(uv run python $X/relevel_iter.py $W/448-lvl0.ntclevel.json:D=1,E=1,F=1 $W/448-lvl1.ntclevel.json:D=1.2084,E=1.2562,F=1.2871 \
      $W/448-lvl2.ntclevel.json:D=1.2218,E=1.2721,F=1.2356 | tee -a "$LOG" | grep -o 'VCC_BULK_CTX_SCALE=[A-Z0-9=.,]*' | cut -d= -f2-)"
log "lvl3 factors $F3"
export VCC_BULK_LOAD="$DUMP:$DSHA" VCC_BULK_SHARPEN_Q=1.75 VCC_BULK_CTX_SCALE="$F3"
( "${DRV[@]}" --name 448-lvl3 --stages build > "$W/driver-lvl3.log" 2>&1 )
[ -s "$W/448-lvl3.h5ad" ] || { log "lvl3 build failed: $(tail -3 $W/driver-lvl3.log)"; exit 1; }
LV3="$(uv run python experiments/417-bulk-sharpen/ntclevel417.py "$W/448-lvl3.h5ad" --bundle "$B" --contexts D,E,F --sd 0.01 2>&1)"
echo "$LV3" >> "$LOG"
log "lvl3 (factors $F3) rho: $(echo "$LV3" | grep -E '^  [DEF]:' | awk '{print $1, $13}' | tr '\n' ' ')"
uv run python $X/relevel_iter.py $W/448-lvl0.ntclevel.json:D=1,E=1,F=1 $W/448-lvl1.ntclevel.json:D=1.2084,E=1.2562,F=1.2871 \
  $W/448-lvl2.ntclevel.json:D=1.2218,E=1.2721,F=1.2356 "$W/448-lvl3.ntclevel.json:$F3" 2>&1 | tail -4 | tee -a "$LOG"
uv run python experiments/428-zero-gene-audit/audit428.py "$W/448-lvl3.h5ad" --bundle "$B" 2>&1 | head -1 | tee -a "$LOG"
rm -f "$W/448-lvl3.h5ad"
log "END"
