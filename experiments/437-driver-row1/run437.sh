#!/usr/bin/env bash
# 437 — row 1 for the NEW driver (432's proposal, installed 2026-09-23; scripts/def_pipeline.py 882a2bf5):
# do the driver's DEFAULTS (SHIPPED_ENV + ABC_ENV), with only the two day-of values exported (VCC_BULK_LOAD,
# VCC_BULK_CTX_SCALE), rebuild the 09-25 candidate 433h byte for byte (artifact sha 5d9e2e38...)?
# Same stage list and flags as run_arm_generic.sh's invocation 1 (no `key` stage on A/B/C), builder restored after.
REPO=~/Code/virtual-cell-challenge-2026; cd "$REPO" || exit 1
X=experiments/437-driver-row1; W=$X/work/ka433h; LOG=$X/work/run437.log; mkdir -p "$W"
log(){ printf '%s %s\n' "$(date -u +%FT%TZ)" "$*" | tee -a "$LOG"; }
for v in $(env | grep -o '^VCC_[A-Z0-9_]*'); do unset "$v"; done
D="$REPO/experiments/417-bulk-sharpen/work/dump_h1x2_cd4x15.npz"; DS=$(shasum -a 256 "$D" | cut -d' ' -f1)
export VCC_BULK_LOAD="$D:$DS" VCC_BULK_CTX_SCALE="A=1.0088,B=1.05,C=1.0736"
log "builder $(shasum -a 256 scripts/make_submission_pm_cd4.py | cut -c1-16) driver $(shasum -a 256 scripts/def_pipeline.py | cut -c1-16) emitter $(shasum -a 256 src/vcc2026/rc_emit_k.py | cut -c1-16); exported: $(env | grep '^VCC_' | tr '\n' ' ')"
/usr/bin/time -l uv run python scripts/def_pipeline.py --bundle data/challenge --work "$W" --name 437-ka-433h \
  --whitener-id union_pruned_def --emitter rc_callset --frozen-inputs "$W/frozen_inputs.json" --min-free-gib 23 \
  --incumbent-meta "$REPO/submissions/433h-cap15-partial50-clean4q80-pertscale3-ballast-sharpen125-stackload-mde010coarse-zr5-ctxlevel.h5ad.meta.json" \
  --stages preflight,identity,cd4-extract,census,whitener,repin-frozen,pin-artifacts,repin-panel-support,build,postbuild \
  --restore-builder-after > "$W/driver.log" 2>&1
RC=$?; log "driver rc=$RC; builder after: $(shasum -a 256 scripts/make_submission_pm_cd4.py | cut -c1-16)"
[ -f "$W/pin-artifacts.json" ] && [ "$(shasum -a 256 scripts/make_submission_pm_cd4.py | cut -c1-8)" != "2c90d9e1" ] && \
  { uv run python scripts/def_pin_panel_artifacts.py --restore --record "$W/pin-artifacts.json" | tee -a "$LOG"; log "builder after explicit restore: $(shasum -a 256 scripts/make_submission_pm_cd4.py | cut -c1-16)"; }
A="$W/437-ka-433h.h5ad"
if [ -s "$A" ]; then
  GOT=$(uv run python -c 'import json,sys; print(json.load(open(sys.argv[1]))["output"]["sha256"])' "$A.meta.json")
  WANT=5d9e2e380bcff3dc21ee1730ed06fe6c2417ba8c1b6859f35da1b2bffc4e188d
  [ "$GOT" = "$WANT" ] && log "ROW1 PASS: artifact sha $GOT == 433h's" || log "ROW1 FAIL: artifact sha $GOT != 433h's $WANT"
  grep -E "\[arms\]|NO  config|ARM DRIFT|checks pass|FAIL" "$W/driver.log" | head -20 | tee -a "$LOG"
  rm -f "$A"
else log "NO ARTIFACT"; tail -20 "$W/driver.log" | tee -a "$LOG"; fi
log "END"
