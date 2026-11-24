#!/usr/bin/env bash
# 670 (coordinator [400b6e], 2026-10-14): 10-15 slot 2 chooser (docs/PLAN-2026-10-15.md, "10-15 s2 rule").
#   nohup bash experiments/670-coord-1015/choose1015s2.sh "<name for 702>" "<name for ku30>" > …/choose1015s2.nohup 2>&1 &
# Reads experiments/701-kup-callset/work/READY-science.log (lines "<ts> READY <tag> <sha256> <package path>"):
#   * the PIE bulk arm (tag pie702*) READY before 18:00Z 10-15 -> it takes slot 2;
#   * else at 18:00Z, K-UP (tag ku30) if READY -> slot 2; if neither, it waits for whichever is READY first until 21:00Z.
# Then it execs arm1015.sh 1015s2 with prior = slot 1's log (submit-1015s1.log). work/NOGO-1015s2 stops it. Never submits itself.
set -uo pipefail
REPO=~/Code/virtual-cell-challenge-2026; cd "$REPO" || exit 1
NAME_PIE="${1:?name for the PIE arm}"; NAME_KU="${2:?name for ku30}"
W=experiments/670-coord-1015/work; LOG=$W/choose1015s2.log; RDY=experiments/701-kup-callset/work/READY-science.log
log(){ printf '%s %s\n' "$(date -u +%FT%TZ)" "$*" >> "$LOG"; }
now(){ date -u +%Y%m%d%H%M; }
echo $$ > $W/choose1015s2.pid
log "armed (pid $$): PIE if READY before 202610151800, else ku30; prior $W/submit-1015s1.log"
pick(){ /usr/bin/grep -E " READY $1[A-Za-z0-9_]* [0-9a-f]{64} " "$RDY" 2>/dev/null | tail -1; }
CH=""
until [ -n "$CH" ]; do
  [ -f "$W/NOGO-1015s2" ] && { log "NOGO-1015s2 present: stopping"; log "END"; exit 1; }
  P="$(pick pie702)"; K="$(pick ku30)"
  if [ -n "$P" ]; then CH="$P"; NAME="$NAME_PIE"; NOTE="pie702 slot 2 on 10-15"
  elif [ "$(now)" -ge 202610151800 ] && [ -n "$K" ]; then CH="$K"; NAME="$NAME_KU"; NOTE="ku30 slot 2 on 10-15"
  elif [ "$(now)" -ge 202610152100 ]; then log "nothing READY by 21:00Z: slot 2 unused"; log "END"; exit 1
  fi
  [ -z "$CH" ] && sleep 60
done
TAG="$(echo "$CH" | awk '{print $3}')"; PKG="$(echo "$CH" | awk '{print $5}')"
RL="experiments/428-zero-gene-audit/work/build_$TAG/run$TAG.log"
log "chosen: $TAG pkg=$PKG run log $RL name='$NAME'"
exec bash experiments/670-coord-1015/arm1015.sh 1015s2 "$TAG" "$RDY" "$RL" "$PKG" "$NAME" "$NOTE" "$W/submit-1015s1.log" 202610150001 202610152300
