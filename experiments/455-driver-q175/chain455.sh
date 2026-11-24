#!/usr/bin/env bash
# 455: after 453's chain ends (it runs the driver), patch the driver's SHIPPED q 1.25 -> 1.75 and prove row 1 again:
# the driver's DEFAULTS plus only the two day-of exports (446b's stack dump and re-level) must rebuild 446b byte for byte
# (artifact sha 9559923a...). On FAIL the driver is restored from its pre-patch copy. Never commits, never submits.
set -uo pipefail
REPO=~/Code/virtual-cell-challenge-2026; cd "$REPO" || exit 1
X=experiments/455-driver-q175; W=$X/work/ka446b; LOG=$X/work/chain455.log; mkdir -p "$W"
log(){ printf '%s %s\n' "$(date -u +%FT%TZ)" "$*" >> "$LOG"; }
echo $$ > $X/work/chain455.pid
log "armed (pid $$): waiting for 453's chain to end"
until /usr/bin/grep -qE " END$|stopping|did not finish|failed" experiments/453-abc-dress-rehearsal/work/chain453.log 2>/dev/null \
      && ! ps -axo comm,args | awk '$1 ~ /^(bash)$/ && /chain453\.sh|def_day\.sh/' | /usr/bin/grep -q .; do
  [ "$(date -u +%Y%m%d%H%M)" -ge 202609261800 ] && { log "453 not done by 18:00Z -- stopping"; exit 1; }
  sleep 60
done
log "453 chain ended: $(tail -1 experiments/453-abc-dress-rehearsal/work/chain453.log | cut -c1-120)"
cp scripts/def_pipeline.py "$X/work/def_pipeline.pre455.py"
for v in $(env | grep -o '^VCC_[A-Z0-9_]*'); do unset "$v"; done
uv run python $X/patch455.py >> "$LOG" 2>&1 || { log "patch failed"; exit 1; }
D="$REPO/experiments/417-bulk-sharpen/work/dump_h1x2_cd4x15.npz"; DS=$(shasum -a 256 "$D" | cut -d' ' -f1)
export VCC_BULK_LOAD="$D:$DS" VCC_BULK_CTX_SCALE="A=1.1004,B=1.1517,C=1.1720"
log "builder $(shasum -a 256 scripts/make_submission_pm_cd4.py | cut -c1-16) driver $(shasum -a 256 scripts/def_pipeline.py | cut -c1-16); exported: $(env | grep '^VCC_' | tr '\n' ' ')"
/usr/bin/time -l uv run python scripts/def_pipeline.py --bundle data/challenge --work "$W" --name 455-ka-446b \
  --whitener-id union_pruned_def --emitter rc_callset --frozen-inputs "$W/frozen_inputs.json" --min-free-gib 23 \
  --incumbent-meta "$REPO/submissions/446b-cap15-partial50-clean4q80-pertscale3-ballast-sharpen175-stackload-ctxlevel-mde010coarse-zr5.h5ad.meta.json" \
  --stages preflight,identity,cd4-extract,census,whitener,repin-frozen,pin-artifacts,repin-panel-support,build,postbuild \
  --restore-builder-after > "$W/driver.log" 2>&1
RC=$?; log "driver rc=$RC; builder after: $(shasum -a 256 scripts/make_submission_pm_cd4.py | cut -c1-16)"
[ -f "$W/pin-artifacts.json" ] && [ "$(shasum -a 256 scripts/make_submission_pm_cd4.py | cut -c1-16)" != "4feae27d37ffa436" ] && \
  { uv run python scripts/def_pin_panel_artifacts.py --restore --record "$W/pin-artifacts.json" >> "$LOG" 2>&1; log "builder after explicit restore: $(shasum -a 256 scripts/make_submission_pm_cd4.py | cut -c1-16)"; }
A="$W/455-ka-446b.h5ad"; WANT=9559923a774eb2489e08a7df38eb6b9b0dbb52241166e3a4c8988793d47457d6
if [ -s "$A" ]; then
  GOT=$(uv run python -c 'import json,sys; print(json.load(open(sys.argv[1]))["output"]["sha256"])' "$A.meta.json")
  if [ "$GOT" = "$WANT" ]; then log "ROW1 PASS: artifact sha $GOT == 446b's"; else log "ROW1 FAIL: artifact sha $GOT != 446b's $WANT"; cp "$X/work/def_pipeline.pre455.py" scripts/def_pipeline.py; log "driver restored to its pre-patch bytes"; fi
  /usr/bin/grep -E "\[arms\]|NO  config|ARM DRIFT|checks pass|FAIL" "$W/driver.log" | head -20 >> "$LOG"
  rm -f "$A"
else log "NO ARTIFACT"; tail -20 "$W/driver.log" >> "$LOG"; cp "$X/work/def_pipeline.pre455.py" scripts/def_pipeline.py; log "driver restored to its pre-patch bytes"; fi
log "driver now $(shasum -a 256 scripts/def_pipeline.py | cut -c1-16)"
log "END"
