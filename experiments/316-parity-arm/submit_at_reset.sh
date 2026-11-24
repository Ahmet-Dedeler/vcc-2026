#!/usr/bin/env bash
# 316 -- submit the parity-correct arm as the FIRST slot of 2026-09-08, then capture the record.
# Preregistration f15ac1d5..., frozen and committed before this runs. Written once; never edited while live.
set -uo pipefail
cd ~/Code/virtual-cell-challenge-2026
LOG=experiments/316-parity-arm/work/submit316.log
log(){ printf '%s %s\n' "$(date -u +%FT%TZ)" "$*" | tee -a "$LOG"; }
log "armed; waiting for the 2026-09-08 slot allowance (00:00Z)"
while [ "$(date -u +%Y%m%d)" -lt 20260908 ]; do sleep 120; done
log "date rolled to $(date -u +%F); submitting"
vcc submit submissions/316-normrel15-pm.vcc -m "a duck reviews the diff" --wait >> "$LOG" 2>&1
RC=$?; log "vcc submit rc=$RC"
ENTRY=$(grep -oE 'entry: [A-Za-z0-9]+' "$LOG" | tail -1 | awk '{print $2}')
log "entry=$ENTRY"
if [ -n "$ENTRY" ]; then
  uv run python scripts/capture_board_entry.py "$ENTRY" --package submissions/316-normrel15-pm.vcc \
     --note "316 parity-correct normrel15-pm; prereg f15ac1d5" >> "$LOG" 2>&1
  log "captured rc=$?"
fi
log "DONE"
