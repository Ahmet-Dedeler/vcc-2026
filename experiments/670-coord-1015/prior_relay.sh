#!/usr/bin/env bash
# 670 (coordinator [3c7a7e], 2026-10-20): prior relay for a chained slot.
# Why: arm1015.sh waits for its prior log's terminal line, then `kill -0`s the FIRST pid that log names. Process ids here wrapped past
# 99999 on 10-20 (new processes got 83xxx while 98102 was still alive), so after the prior arm1015 exits, a new process can take its
# number and hold the next slot until its cap. arm1015.sh is executed by live waiters and cannot be edited. Instead, the next slot's
# prior is a STUB that this relay writes only once the source slot is really over: its submit log has a terminal line AND no bash
# arm1015.sh <SRC-SL> process is alive (checked by command name and arguments, not by pid). The stub carries the source's entry id
# (so arm1015 can capture or skip it) and an END line, and never contains a process id.
#   nohup bash scripts/detach.sh bash experiments/670-coord-1015/prior_relay.sh <SRC-SL> <STUB> <CAP YYYYmmddHHMM> > …/prior_relay-<SRC-SL>.nohup 2>&1 &
set -uo pipefail
REPO=~/Code/virtual-cell-challenge-2026; cd "$REPO" || exit 1
SRC="${1:?source slot label}"; STUB="${2:?stub path}"; CAP="${3:?cap}"
W="${PR_W:-experiments/670-coord-1015/work}"; S="$W/submit-$SRC.log"; LOG="$W/prior_relay-$SRC.log"
log(){ printf '%s %s\n' "$(date -u +%FT%TZ)" "$*" >> "$LOG"; }
terminal(){ /usr/bin/grep -qE " END$|SUBMIT rc=|NOT SUBMITTING|ERROR " "$S" 2>/dev/null; }
live(){ ps -axo comm,args | awk -v s="$SRC" '$1 ~ /^bash$/ { for (i = 2; i < NF; i++) if ($i ~ /arm1015\.sh$/ && $(i+1) == s) n++ } END {exit n ? 0 : 1}'; }
[ -e "$STUB" ] && { echo "stub $STUB exists already: refusing to overwrite"; exit 2; }
log "armed: relays $S into $STUB once it is terminal and no arm1015 $SRC is alive; cap $CAP"
while :; do
  if terminal && ! live; then
    sleep 5; terminal && ! live || continue
    ID="$(/usr/bin/grep -oE 'entry[: ]+[A-Za-z0-9]{20}|entry id: [A-Za-z0-9]{20}' "$S" | tail -1 | /usr/bin/grep -oE '[A-Za-z0-9]{20}$')"
    T="$(date -u +%FT%TZ)"
    if [ -n "$ID" ]; then printf '%s relay of slot %s: terminal, entry %s\n%s END\n' "$T" "$SRC" "$ID" "$T" > "$STUB"
    else printf '%s relay of slot %s: terminal, it did not submit\n%s END\n' "$T" "$SRC" "$T" > "$STUB"; fi
    log "wrote $STUB (${ID:-no entry}): $(tail -1 "$S" | cut -c1-160)"; log "END"; exit 0
  fi
  [ "$(date -u +%Y%m%d%H%M)" -ge "$CAP" ] && { log "cap $CAP reached; no stub written (the next slot waits to its own cap)"; log "END"; exit 1; }
  sleep 30
done
