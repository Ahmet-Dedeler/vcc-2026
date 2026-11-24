#!/usr/bin/env bash
# 670 (coordinator [2d6f46], 2026-10-16): delete a prediction h5ad once its package is READY (AGENTS §10: only after the .vcc exists and
# its sha is logged). It waits for "<ts> READY <TAG> <sha64> <pkg>" in the READY log, requires sha256(<pkg>) == <sha64>, requires that no
# process holds the h5ad open (lsof), then removes it and logs to work/SWEEP.log. Any mismatch → it keeps the file and says why.
#   nohup bash scripts/detach.sh bash experiments/670-coord-1015/sweep_after_ready.sh <TAG> <h5ad> <CAP YYYYmmddHHMM> > …/sweep-<TAG>.nohup 2>&1 &
set -uo pipefail
REPO=~/Code/virtual-cell-challenge-2026; cd "$REPO" || exit 1
TAG="${1:?tag}"; H5="${2:?h5ad}"; CAP="${3:?cap}"
RDY="${SW_READY:-experiments/701-kup-callset/work/READY-science.log}"; W="${SW_W:-experiments/670-coord-1015/work}"
LOG="$W/sweep-$TAG.log"
log(){ printf '%s %s\n' "$(date -u +%FT%TZ)" "$*" >> "$LOG"; }
echo $$ > "$W/sweep-$TAG.pid"; log "armed (pid $$): sweep $H5 after READY $TAG in $RDY (cap $CAP)"
until L="$(/usr/bin/grep -E " READY $TAG [0-9a-f]{64} " "$RDY" 2>/dev/null | tail -1)"; [ -n "$L" ]; do
  [ "$(date -u +%Y%m%d%H%M)" -ge "$CAP" ] && { log "no READY $TAG by cap: kept $H5"; log "END"; exit 0; }
  sleep 60; done
SHA="$(echo "$L" | awk '{print $4}')"; PKG="$(echo "$L" | awk '{print $5}')"
[ -s "$PKG" ] || { log "package $PKG missing: kept $H5"; log "END"; exit 1; }
GOT="$(shasum -a 256 "$PKG" | cut -d' ' -f1)"
[ "$GOT" = "$SHA" ] || { log "package sha $GOT != READY $SHA: kept $H5"; log "END"; exit 1; }
[ -f "$H5" ] || { log "$H5 already gone"; log "END"; exit 0; }
for i in 1 2 3 4 5 6 7 8 9 10; do lsof "$H5" >/dev/null 2>&1 || break; sleep 60; done
lsof "$H5" >/dev/null 2>&1 && { log "$H5 still open after 10 min: kept"; log "END"; exit 1; }
SZ="$(stat -f %z "$H5")"; rm -f "$H5" && printf '%s swept %s (%s bytes): .vcc %s sha %s == READY %s line; restorable via scripts/restore_prediction_h5ad.py (sweep_after_ready)\n' \
  "$(date -u +%FT%TZ)" "$H5" "$SZ" "$PKG" "$SHA" "$TAG" >> "$W/SWEEP.log"
log "swept $H5"; log "END"
