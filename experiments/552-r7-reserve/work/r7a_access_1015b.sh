#!/usr/bin/env bash
# R7a access, v2 (coordinator [400b6e]'s guard, 10-14 23:5xZ): start r7a551.py access only when
#   (1) no bash process is running experiments/702-pie-bulk-arm/build702.sh (its r1 -> ntclevel -> r2 has builder-free
#       moments; R7a's 14 GiB fits must not starve r2's 23 GiB gate), matched as interpreter + path (AGENTS 6.2);
#   (2) no bash process is running experiments/671-dnmag-arm/dn671.sh (DN build + its KAu);
#   (3) def_pipeline.other_builders() is empty (no def_pipeline build in its build stage).
# Then 551's own gate admits each heavy step (22 GiB free+inactive, pressure 55 %, swap 400 MiB, the recorded window).
# Resumable. Start cap: the ledger's window extension minus 30 min (10-16 07:3xZ, [d80a10]; was a fixed 12:00Z). Terminal line: "<ts> END rc=<rc>".
set -uo pipefail
REPO=~/Code/virtual-cell-challenge-2026; cd "$REPO" || exit 1
LOG=experiments/551-r6-c10b/work/r7a/access_wrapper_1015.log
log(){ printf '%s %s\n' "$(date -u +%FT%TZ)" "$*" >> "$LOG"; }
alive(){ ps -axo comm,args | awk -v p="$1" '$1 ~ /^bash$/ && index($0, p) {n++} END {exit n ? 0 : 1}'; }
# start cap = the ledger R7 section's window extension minus 30 min (10-16, [d80a10]; fallback 10-16 17:30Z)
r7cap(){ local w c
  w="$(awk '/^## Reserve R7/{s=1;next} s&&/^## /{s=0} s' docs/CONFIRMATION-LEDGER.md \
       | /usr/bin/grep -oE 'R7a WINDOW END EXTENDED TO [0-9]{4}-[0-9]{2}-[0-9]{2}T[0-9]{2}:[0-9]{2}Z' | tail -1 | awk '{print $NF}')"
  c="$(date -j -u -v-30M -f '%Y-%m-%dT%H:%MZ' "$w" +%Y%m%d%H%M 2>/dev/null)"; echo "${c:-202610161730}"; }
log "v2 armed (pid $$): waiting for build702.sh and dn671.sh to exit and no def_pipeline build"
until ! alive experiments/702-pie-bulk-arm/build702.sh && ! alive experiments/671-dnmag-arm/dn671.sh \
      && ! uv run python -c 'import sys; sys.path.insert(0,"scripts"); import def_pipeline as P; raise SystemExit(0 if P.other_builders() else 1)' 2>/dev/null; do
  CAP="$(r7cap)"; [ "$(date -u +%Y%m%d%H%M)" -ge "$CAP" ] && { log "start cap $CAP passed before the conditions held (ledger window minus 30 min)"; log "END rc=9"; exit 9; }
  sleep 60
done
log "conditions hold: starting r7a551.py access"
uv run python experiments/551-r6-c10b/r7a551.py access > experiments/551-r6-c10b/work/r7a/access.log 2>&1
rc=$?
log "access rc=$rc: $(tail -3 experiments/551-r6-c10b/work/r7a/access.log | tr '\n' ' ' | cut -c1-400)"
log "END rc=$rc"
