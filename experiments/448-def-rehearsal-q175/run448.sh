#!/usr/bin/env bash
# 448 — October rehearsal of the CURRENT recipe on 238's essential-free realistic-coverage synthetic panel (ef-bundle-r),
# the panel whose re-level needs factors > 1.2 at q 1.75 (438's level file): driver setup stages (day's operator, day's key,
# bare-env pins) -> stack dump (H1 x2, CD4 x1 = 438's October default) -> un-levelled LOAD build at q 1.75 (+ MDE coarse,
# zero repair from SHIPPED_ENV) -> ntclevel417 (sd 0.01) -> levelled build -> ntclevel verify (one correction pass if
# any context is > 0.5 % off in scale) -> audits -> builder restore. Exercises 447's widened bound end to end.
set -uo pipefail
REPO=~/Code/virtual-cell-challenge-2026; cd "$REPO" || exit 1
X=experiments/448-def-rehearsal-q175; W=$X/work/efr; LOG=$X/work/run448.log; mkdir -p "$W"
B=experiments/238-essential-free-panel/work/ef-bundle-r
log(){ printf '%s %s\n' "$(date -u +%FT%TZ)" "$*" | tee -a "$LOG"; }
for v in $(env | grep -o '^VCC_[A-Z0-9_]*'); do unset "$v"; done
BSHA0="$(shasum -a 256 scripts/make_submission_pm_cd4.py | cut -d' ' -f1)"
log "start: builder ${BSHA0:0:16} driver $(shasum -a 256 scripts/def_pipeline.py | cut -c1-16); exported VCC_*: $(env | grep -c '^VCC_')"
[ "$BSHA0" = "b6c3f13141e707012588cc354f62758889db9c38c3757f9765ebea37a8a733f8" ] || { log "builder is not b6c3f131 -- stopping"; exit 1; }
DRV=(uv run python scripts/def_pipeline.py --bundle "$B" --work "$W" --whitener-id union_pruned_def --emitter rc_callset
     --frozen-inputs "$W/frozen_inputs.json" --min-free-gib 23)
restore(){ [ -f "$W/pin-artifacts.json" ] && uv run python scripts/def_pin_panel_artifacts.py --restore --record "$W/pin-artifacts.json" >> "$LOG" 2>&1
           log "builder after restore: $(shasum -a 256 scripts/make_submission_pm_cd4.py | cut -c1-16) (want b6c3f1314)"; }
trap restore EXIT
# 1. setup: light stages, the day's operator and key, bare-env pins (no build, no restore yet)
/usr/bin/time -l "${DRV[@]}" --name 448-setup \
  --stages preflight,identity,cd4-extract,census,whitener,key,repin-frozen,pin-artifacts,repin-panel-support > "$W/driver-setup.log" 2>&1
log "setup rc=$?: $(grep -E '^--- \[' "$W/driver-setup.log" | tr '\n' ' ' | cut -c1-600)"
grep -q "^--- \[repin-panel-support\].*rc=0" "$W/driver-setup.log" || { log "setup did not complete -- stopping"; exit 1; }
# 2. the stack dump: H1 x2 (VCC_H1_UNION_G=0.56), CD4 at the DEF_ENV default (x1)
DUMP="$W/dump_h1x2_cd4x1.npz"
( export VCC_H1_UNION_G=0.56 VCC_BULK_DUMP="$DUMP"; "${DRV[@]}" --name 448-dump --stages build > "$W/driver-dump.log" 2>&1 )
[ -s "$DUMP" ] || { log "no dump written -- stopping ($(grep -h 'bulk-dump\|Error\|REFUS' $W/driver-dump.log $W/logs/build.log 2>/dev/null | tail -2))"; exit 1; }
DSHA="$(shasum -a 256 "$DUMP" | cut -d' ' -f1)"; log "dump ${DSHA:0:16}: $(grep -h '\[h1-union\]\|\[cd4\] sizing\|bulk-dump' $W/logs/build.log 2>/dev/null | tail -3 | tr '\n' ' ' | cut -c1-300)"
# 3. the un-levelled LOAD build at q 1.75
export VCC_BULK_LOAD="$DUMP:$DSHA" VCC_BULK_SHARPEN_Q=1.75
( "${DRV[@]}" --name 448-lvl0 --stages build > "$W/driver-lvl0.log" 2>&1 )
[ -s "$W/448-lvl0.h5ad" ] || { log "lvl0 build failed: $(tail -3 $W/driver-lvl0.log)"; exit 1; }
LV0="$(uv run python experiments/417-bulk-sharpen/ntclevel417.py "$W/448-lvl0.h5ad" --bundle "$B" --contexts D,E,F --sd 0.01 2>&1)"
echo "$LV0" >> "$LOG"; S1="$(echo "$LV0" | grep -o 'VCC_BULK_CTX_SCALE (relative to THIS artifact) = [A-Z0-9=.,]*' | awk '{print $NF}')"
log "lvl0 rho/cliffs: $(echo "$LV0" | grep -E '^  [DEF]:' | awk '{print $1, "rho", $13, "cliff", $15}' | tr '\n' ' ') -> factors $S1"
uv run python experiments/428-zero-gene-audit/audit428.py "$W/448-lvl0.h5ad" --bundle "$B" 2>&1 | head -1 | tee -a "$LOG"
rm -f "$W/448-lvl0.h5ad"
[ -n "$S1" ] || { log "no re-level string -- stopping"; exit 1; }
# 4. the levelled build (the string is relative to lvl0, i.e. the factors themselves)
export VCC_BULK_CTX_SCALE="$S1"
( VCC_RC_AUDIT_DUMP="$W/rc_audit_lvl1.jsonl" "${DRV[@]}" --name 448-lvl1 --stages build > "$W/driver-lvl1.log" 2>&1 )
[ -s "$W/448-lvl1.h5ad" ] || { log "lvl1 build failed: $(tail -3 $W/driver-lvl1.log)"; exit 1; }
LV1="$(uv run python experiments/417-bulk-sharpen/ntclevel417.py "$W/448-lvl1.h5ad" --bundle "$B" --contexts D,E,F --sd 0.01 2>&1)"
echo "$LV1" >> "$LOG"; S2="$(echo "$LV1" | grep -o 'VCC_BULK_CTX_SCALE (relative to THIS artifact) = [A-Z0-9=.,]*' | awk '{print $NF}')"
log "lvl1 (factors $S1) rho: $(echo "$LV1" | grep -E '^  [DEF]:' | awk '{print $1, $13}' | tr '\n' ' ') residual EV-scale $S2"
uv run python experiments/428-zero-gene-audit/audit428.py "$W/448-lvl1.h5ad" --bundle "$B" 2>&1 | head -1 | tee -a "$LOG"
uv run python experiments/431-rc-audit-dump/analyse431.py "$W/rc_audit_lvl1.jsonl" 2>&1 | grep -i "mde_lambda\|median" | head -4 | tee -a "$LOG"
grep -h "\[arms\]\|REFUSING\|bulk-ctx-scale\|bulk-load" "$W/driver-lvl1.log" "$W/logs/build.log" 2>/dev/null | head -8 | cut -c1-200 | tee -a "$LOG"
# 5. one correction pass if any context is > 0.5 % off in scale (rho ~ f^1.86 compression above 1.1)
NEED="$(python3 -c "
s='$S2'; d={k:float(v) for k,v in (p.split('=') for p in s.split(','))}
print('yes' if any(abs(v-1)>0.005 for v in d.values()) else 'no')")"
if [ "$NEED" = "yes" ]; then
  S3="$(python3 -c "
a={k:float(v) for k,v in (p.split('=') for p in '$S1'.split(','))}; b={k:float(v) for k,v in (p.split('=') for p in '$S2'.split(','))}
print(','.join(f'{k}={a[k]*b[k]**(2/1.86):.4f}' for k in a))")"
  log "correction pass: $S1 x residual $S2 (compression-adjusted) -> $S3"
  export VCC_BULK_CTX_SCALE="$S3"
  rm -f "$W/448-lvl1.h5ad"
  ( "${DRV[@]}" --name 448-lvl2 --stages build > "$W/driver-lvl2.log" 2>&1 )
  [ -s "$W/448-lvl2.h5ad" ] && { LV2="$(uv run python experiments/417-bulk-sharpen/ntclevel417.py "$W/448-lvl2.h5ad" --bundle "$B" --contexts D,E,F --sd 0.01 2>&1)"; echo "$LV2" >> "$LOG"
    log "lvl2 (factors $S3) rho: $(echo "$LV2" | grep -E '^  [DEF]:' | awk '{print $1, $13}' | tr '\n' ' ') residual $(echo "$LV2" | grep -o 'VCC_BULK_CTX_SCALE (relative to THIS artifact) = [A-Z0-9=.,]*' | awk '{print $NF}')"; rm -f "$W/448-lvl2.h5ad"; } \
    || log "lvl2 build failed: $(tail -3 $W/driver-lvl2.log)"
else
  log "no correction needed (all contexts within 0.5 % of the EV scale)"; rm -f "$W/448-lvl1.h5ad"
fi
log "END"
