#!/usr/bin/env bash
# 439 — the INSTALLED October driver (5463fe28) end to end on an October-like synthetic panel (238's essential-free,
# realistic-coverage `ef-bundle-r`), with nothing exported: every September knob from the driver's SHIPPED_ENV + DEF_ENV,
# the day's operator from the `whitener` stage, the day's key from the `key` stage (scripts/def_pert_scale_key.py),
# its registration by `pin-artifacts --pert-scale`, bare-env pins, and the guarded restore in `finally`.
# 432 rehearsed the same stack with its own wrappers; this is the driver's own path. No LOAD/CTX (the b0 analogue).
REPO=~/Code/virtual-cell-challenge-2026; cd "$REPO" || exit 1
X=experiments/439-def-driver-e2e; W=$X/work/efr; LOG=$X/work/run439.log; mkdir -p "$W"
B=experiments/238-essential-free-panel/work/ef-bundle-r
log(){ printf '%s %s\n' "$(date -u +%FT%TZ)" "$*" | tee -a "$LOG"; }
for v in $(env | grep -o '^VCC_[A-Z0-9_]*'); do unset "$v"; done
log "start: builder $(shasum -a 256 scripts/make_submission_pm_cd4.py | cut -c1-16) driver $(shasum -a 256 scripts/def_pipeline.py | cut -c1-16) pin $(shasum -a 256 scripts/def_pin_panel_artifacts.py | cut -c1-16) key $(shasum -a 256 scripts/def_pert_scale_key.py | cut -c1-16); exported VCC_*: $(env | grep -c '^VCC_')"
/usr/bin/time -l uv run python scripts/def_pipeline.py --bundle "$B" --work "$W" --name 439-efr \
  --whitener-id union_pruned_def --emitter rc_callset --frozen-inputs "$W/frozen_inputs.json" --min-free-gib 23 \
  --stages preflight,identity,cd4-extract,census,whitener,key,repin-frozen,pin-artifacts,repin-panel-support,build,postbuild \
  --restore-builder-after > "$W/driver.log" 2>&1
RC=$?; log "driver rc=$RC; builder after: $(shasum -a 256 scripts/make_submission_pm_cd4.py | cut -c1-16) (want 63dfd7d0e4027388)"
grep -E "^--- \[|stage .* failed|\[arms\]|NO  config|REFUSING|restore|not restoring|checks pass|FAIL" "$W/driver.log" | cut -c1-220 | tee -a "$LOG"
A="$W/439-efr.h5ad"; [ -s "$A" ] && { log "artifact $(ls -la "$A" | awk '{print $5}') bytes"; \
  uv run python experiments/428-zero-gene-audit/audit428.py "$A" --bundle "$B" 2>&1 | head -1 | tee -a "$LOG"; rm -f "$A"; log "artifact deleted (the meta stays)"; }
log "END"
