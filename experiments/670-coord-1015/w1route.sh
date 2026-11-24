#!/usr/bin/env bash
# 670 (coordinator [3c7a7e], 2026-10-20): the preregistered W1 BELOW route, made automatic (SIXA-DECISION Amendment 2; 673 PREREGISTRATION).
# Why: if W1 (10-21 s1, H1 x2.5 on R0's recipe) reads BELOW, wait708 must be re-armed with WAIT708_R710=HARM BEFORE the release, so that
# decide708 itself blocks any H1 > x2 move. Cross-session messaging was paused on 10-20 (too many sends without a user message), and the
# coordinator's session may be asleep at the read (~10-21 01:30Z), so this does it without a message.
# It waits for submit-1021s1.log's "CAPTURED own entry <id>", reads it with read695.py --w1 (against R-ABC2's entry), and:
#   - COMPOSES / ABOVE: logs the bin and exits. Nothing else is touched.
#   - BELOW, and no release yet: stops wait708 (matched by command name + script path, never by a pid from a log), relaunches it from
#     source with WAIT708_R710=HARM (OCT22-CHECKLIST 4b's command), and checks that its new armed line says R710=HARM.
#   - BELOW after the release, or if wait708 already started day708: does NOT touch wait708; logs "APPLY BY HAND" (SIXA backstop).
#   - no capture by the cap / the slot did not submit: logs it and exits.
#   DRY=1: prints what it would do instead of stopping or launching anything (for tests).
#   nohup bash scripts/detach.sh bash experiments/670-coord-1015/w1route.sh <CAP YYYYmmddHHMM> > experiments/670-coord-1015/work/w1route.nohup 2>&1 &
set -uo pipefail
REPO=~/Code/virtual-cell-challenge-2026; cd "$REPO" || exit 1
CAP="${1:?cap}"; DRY="${DRY:-0}"
W="${WR_W:-experiments/670-coord-1015/work}"; S="${WR_SUBMIT:-$W/submit-1021s1.log}"; LOG="${WR_LOG:-$W/w1route.log}"
R2="${WR_RABC2:-data/board-entries/nnHQfcZHsT3as9EafA3i.json}"; MAN="${WR_MAN:-data/challenge-def/manifest.json}"
W708=experiments/690-def-readiness-1015/wait708.sh; W708LOG=experiments/690-def-readiness-1015/work/wait708.log
DAYLOG=experiments/708-def-panel-repricing/work/day708-d708def.log
log(){ printf '%s %s\n' "$(date -u +%FT%TZ)" "$*" >> "$LOG"; }
w708pids(){ ps -axo pid,comm,args | awk '$2 ~ /^bash$/ && /wait708\.sh/ && !/w1route/ {print $1}'; }
log "armed (pid $$, DRY=$DRY): waits for W1's capture in $S; BELOW re-arms wait708 with WAIT708_R710=HARM before the release; cap $CAP"
until ID="$(/usr/bin/grep -oE 'CAPTURED own entry [A-Za-z0-9]{20}' "$S" 2>/dev/null | tail -1 | awk '{print $4}')"; [ -n "$ID" ]; do
  if /usr/bin/grep -qE "NOT SUBMITTING|ERROR submit|WARN capture" "$S" 2>/dev/null && ! /usr/bin/grep -q "CAPTURED own entry" "$S"; then
    log "W1's slot ended without a capture: $(tail -2 "$S" | tr '\n' ' ' | cut -c1-200)"; log "END"; exit 1; fi
  [ "$(date -u +%Y%m%d%H%M)" -ge "$CAP" ] && { log "cap $CAP reached with no W1 capture"; log "END"; exit 1; }
  sleep 60; done
OUT="$(uv run python experiments/670-coord-1015/read695.py --w1 "data/board-entries/$ID.json" "$R2" 2>&1)"; RC=$?
log "read695 --w1 $ID rc $RC: $(printf '%s' "$OUT" | tr '\n' ' ' | cut -c1-400)"
[ $RC -eq 0 ] || { log "reader failed: nothing changed; APPLY BY HAND"; log "END"; exit 1; }
V="$(printf '%s' "$OUT" | python3 -c 'import json,sys; print(json.load(sys.stdin)["verdict"])')"
case "$V" in BELOW*) ;; *) log "W1 bin: $V -> no change to wait708"; log "END"; exit 0;; esac
log "W1 bin: BELOW -> the H1 > x2 block applies (SIXA Amendment 2)"
if [ -s "$MAN" ] || [ -e "$DAYLOG" ]; then log "release present or day708 started: NOT touching wait708; APPLY BY HAND at SIXA time (backstop)"; log "END"; exit 3; fi
P="$(w708pids)"; log "wait708 now: ${P:-none alive}"
# readiness may re-arm wait708 by hand on the same read: if a live wait708's latest armed line already says R710=HARM, do nothing
if [ -n "$P" ] && /usr/bin/grep 'wait708 armed' "$W708LOG" | tail -1 | /usr/bin/grep -q 'R710=HARM'; then
  log "wait708's latest armed line already says R710=HARM (re-armed by hand?): nothing to do"; log "END"; exit 0; fi
if [ "$DRY" = 1 ]; then log "DRY: would stop ${P:-nothing} and launch: WAIT708_R710=HARM nohup bash scripts/detach.sh bash $W708 >> experiments/690-def-readiness-1015/work/wait708.nohup"; log "END"; exit 0; fi
for p in $P; do kill "$p" 2>/dev/null && log "stopped wait708 $p"; done; sleep 3
[ -z "$(w708pids)" ] || { log "wait708 still alive after the stop: APPLY BY HAND"; log "END"; exit 1; }
T0="$(date -u +%FT%TZ)"
WAIT708_R710=HARM nohup bash scripts/detach.sh bash "$W708" >> experiments/690-def-readiness-1015/work/wait708.nohup 2>&1 &
armed_harm(){ awk -v t0="$T0" 'substr($0, 1, 20) >= t0 && /wait708 armed/ && /R710=HARM/ {f = 1} END {exit f ? 0 : 1}' "$W708LOG"; }
sleep 15
if ! armed_harm && [ -z "$(w708pids)" ]; then   # it never came up: one more launch, so wait708 is never left dead
  log "no wait708 alive and no HARM armed line 15 s after the launch: launching once more"
  WAIT708_R710=HARM nohup bash scripts/detach.sh bash "$W708" >> experiments/690-def-readiness-1015/work/wait708.nohup 2>&1 &
  sleep 30; fi
N="$(w708pids | wc -l | tr -d ' ')"
if [ "${N:-0}" -gt 1 ]; then   # a hand re-arm raced this one: keep only the newest (highest start time), all of them read R710 from the same env rule
  KEEP="$(ps -axo pid=,lstart=,comm=,args= | awk '$7 ~ /^bash$/ && /wait708\.sh/ && !/w1route/' | python3 -c '
import sys, datetime as dt
rows = [l.split(None, 6) for l in sys.stdin if l.strip()]
print(max(rows, key=lambda r: dt.datetime.strptime(" ".join(r[1:6]), "%a %b %d %H:%M:%S %Y"))[0] if rows else "")')"
  for p in $(w708pids); do [ "$p" = "$KEEP" ] || { kill "$p" 2>/dev/null && log "stopped a duplicate wait708 $p (kept $KEEP)"; }; done; fi
if armed_harm; then
  log "wait708 re-armed with R710=HARM: $(/usr/bin/grep 'wait708 armed' "$W708LOG" | tail -1 | cut -c1-200)"
else log "re-armed wait708 shows no R710=HARM armed line after $T0 (alive: $(w708pids | tr '\n' ' ')): CHECK BY HAND"; fi
log "END"
