#!/usr/bin/env bash
# 463: adopt the zero-repair floor -2 into the October driver's defaults, IF 446br2 reads floor_ok on 09-27 (449's bins:
# nmae raw delta <= -0.002 vs 446b). Launched by hand after the read, with the board delta text as the argument:
#   bash experiments/463-floor-adopt/chain463.sh "+0.00xx overall, nmae -0.00xx raw (446br2, entry <id>)"
# Then 455's pattern: patch SHIPPED_ENV VCC_DE_ZERO_REPAIR -5 -> -2 with scripts/shipped_env_update.py, and prove row 1
# again: the driver's DEFAULTS plus only 446b's two day-of exports (stack dump, re-level) must rebuild 446br2 byte for byte
# (artifact sha 5581e305...). On FAIL the driver is restored from its pre-patch copy. Waits for any running build.
# Never commits, never submits.
set -uo pipefail
REPO=~/Code/virtual-cell-challenge-2026; cd "$REPO" || exit 1
TXT="${1:?board delta text for ARM_BOARD_DELTA}"
X=experiments/463-floor-adopt; W=$X/work/ka446br2; LOG=$X/work/chain463.log; mkdir -p "$W"
log(){ printf '%s %s\n' "$(date -u +%FT%TZ)" "$*" >> "$LOG"; }
echo $$ > $X/work/chain463.pid
log "armed (pid $$): waiting for no running build"
until ! ps -axo comm,args | awk '$1 ~ /^bash$/ && /def_day\.sh|run_arm_generic\.sh|chain46[01]\.sh|chain458d\.sh/' | /usr/bin/grep -q .; do
  [ "$(date -u +%Y%m%d%H%M)" -ge 202609272000 ] && { log "a build still running at 20:00Z 09-27 -- stopping"; exit 1; }; sleep 60; done
BF="$(shasum -a 256 scripts/make_submission_pm_cd4.py | cut -d' ' -f1)"
[ "$BF" = "$(git show HEAD:scripts/make_submission_pm_cd4.py | shasum -a 256 | cut -d' ' -f1)" ] || { log "builder != HEAD (pinned?) -- stopping"; exit 1; }
cp scripts/def_pipeline.py "$X/work/def_pipeline.pre463.py"
for v in $(env | grep -o '^VCC_[A-Z0-9_]*'); do unset "$v"; done
uv run python scripts/shipped_env_update.py VCC_DE_ZERO_REPAIR=-2 "zero-repair floor -2: $TXT" >> "$LOG" 2>&1 || { log "patch failed"; exit 1; }
log "driver patched: $(shasum -a 256 scripts/def_pipeline.py | cut -c1-16) (was $(shasum -a 256 $X/work/def_pipeline.pre463.py | cut -c1-16))"
D="$REPO/experiments/417-bulk-sharpen/work/dump_h1x2_cd4x15.npz"; DS=$(shasum -a 256 "$D" | cut -d' ' -f1)
export VCC_BULK_LOAD="$D:$DS" VCC_BULK_CTX_SCALE="A=1.1004,B=1.1517,C=1.1720"
log "builder $(shasum -a 256 scripts/make_submission_pm_cd4.py | cut -c1-16); exported: $(env | grep '^VCC_' | tr '\n' ' ')"
/usr/bin/time -l uv run python scripts/def_pipeline.py --bundle data/challenge --work "$W" --name 463-ka-446br2 \
  --whitener-id union_pruned_def --emitter rc_callset --frozen-inputs "$W/frozen_inputs.json" --min-free-gib 23 \
  --incumbent-meta "$REPO/submissions/446br2-cap15-partial50-clean4q80-pertscale3-ballast-sharpen175-stackload-ctxlevel-mde010coarse-zr2.h5ad.meta.json" \
  --stages preflight,identity,cd4-extract,census,whitener,repin-frozen,pin-artifacts,repin-panel-support,build,postbuild \
  --restore-builder-after > "$W/driver.log" 2>&1
RC=$?; log "driver rc=$RC; builder after: $(shasum -a 256 scripts/make_submission_pm_cd4.py | cut -c1-16)"
[ -f "$W/pin-artifacts.json" ] && [ "$(shasum -a 256 scripts/make_submission_pm_cd4.py | cut -d' ' -f1)" != "$BF" ] && \
  { uv run python scripts/def_pin_panel_artifacts.py --restore --record "$W/pin-artifacts.json" >> "$LOG" 2>&1; log "builder after explicit restore: $(shasum -a 256 scripts/make_submission_pm_cd4.py | cut -c1-16)"; }
A="$W/463-ka-446br2.h5ad"; WANT=5581e3057b523ad3c3da4de2bdcd7cc44e8f31d0291784d012b9c7145e24b49c
if [ -s "$A" ]; then
  GOT=$(uv run python -c 'import json,sys; print(json.load(open(sys.argv[1]))["output"]["sha256"])' "$A.meta.json")
  if [ "$GOT" = "$WANT" ]; then log "ROW1 PASS: artifact sha $GOT == 446br2's"; else log "ROW1 FAIL: artifact sha $GOT != 446br2's $WANT"; cp "$X/work/def_pipeline.pre463.py" scripts/def_pipeline.py; log "driver restored to its pre-patch bytes"; fi
  /usr/bin/grep -E "\[arms\]|NO  config|ARM DRIFT|checks pass|FAIL" "$W/driver.log" | head -20 >> "$LOG"
  rm -f "$A"
else log "NO ARTIFACT"; tail -20 "$W/driver.log" >> "$LOG"; cp "$X/work/def_pipeline.pre463.py" scripts/def_pipeline.py; log "driver restored to its pre-patch bytes"; fi
log "driver now $(shasum -a 256 scripts/def_pipeline.py | cut -c1-16)"
log "END"
