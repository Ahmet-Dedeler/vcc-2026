#!/usr/bin/env bash
# 670 (coordinator [2d6f46], 2026-10-16) v2 of release_hold690: lift the 690 HOLD only after R7a's wrapper (r7a_after705.sh, armed 04:36:32Z,
# pid 75920) has EXITED and written its terminal " END" line, per R7 owner [d80a10]: order build705 -> R7a -> run690. r7a_after705 always ends
# with "<ts> END" (also when it skips because build705 ran late or the release landed). Cap 10-18 12:00Z (HOLD then stays; coordinator decides).
set -uo pipefail
REPO=~/Code/virtual-cell-challenge-2026; cd "$REPO" || exit 1
LOG=experiments/670-coord-1015/work/release_hold690.log; HOLD=experiments/690-def-readiness-1015/work/HOLD
RL=experiments/552-r7-reserve/work/r7a_after705.log
log(){ printf '%s %s\n' "$(date -u +%FT%TZ)" "$*" >> "$LOG"; }
alive(){ ps -axo comm,args | awk '$1 ~ /^bash$/ && /experiments\/552-r7-reserve\/work\/r7a_after705\.sh/' | /usr/bin/grep -q .; }
cap(){ [ "$(date -u +%Y%m%d%H%M)" -ge 202610181200 ]; }
log "v2 armed (pid $$): HOLD stays until r7a_after705.sh exits with a terminal END after its 04:36:32Z arm"
until ! alive && awk '$1 > "2026-10-16T04:36:32Z" && / END$/' "$RL" | /usr/bin/grep -q .; do
  cap && { log "cap reached; HOLD left in place"; log "END"; exit 1; }; sleep 60; done
log "R7a wrapper ended: $(tail -2 "$RL" | tr '\n' ' ' | cut -c1-200)"
rm -f "$HOLD" && log "HOLD removed: run690 (R-NEW + R-KA) may start"
log "END"
