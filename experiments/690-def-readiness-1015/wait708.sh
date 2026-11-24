#!/usr/bin/env bash
# 690 wait708 (readiness [1efb7d], 2026-10-19; the coordinator [223d39] asked for a detached waiter, since sessions can be gone on the day).
# Runs 708's day-of re-pricing on D/E/F's real panel ONCE, after BOTH:
#   (a) a REAL release: data/challenge-def/manifest.json exists. T_REL = the real watcher's "RELEASE DETECTED" line in
#       experiments/691-release-watcher/work/watch691.log (else the manifest's mtime). Rehearsal watchers write elsewhere.
#   (b) R0 is done, judged on experiments/670-coord-1015/work/submit-R0.log lines LATER than T_REL:
#       - "SUBMITTED ok" (R0 uploaded), or
#       - the last line ends " END" (NOT SUBMITTING / ERROR / a NOT READY build) AND no bash `arm1015.sh R0` process is alive, on two
#         checks 120 s apart. name_guard kills and relaunches arm1015 after a content-filter refusal; a killed arm leaves no END line.
#   So day708 never competes with R0's build or upload. day708 itself waits for a free machine at each step (gate602).
# Then (same command as OCT22-CHECKLIST 4b): R705=CLEARED R707=WEAK R710=POSITIVE bash experiments/708-def-panel-repricing/day708.sh
#   data/challenge-def experiments/692-def-r0/work/r0 d708def
#   -> experiments/708-def-panel-repricing/work/decision708-d708def.txt (day708 log: work/day708-d708def.log).
# It writes decisions only: no build, no submit; docs/october-decisions.json is untouched. The coordinator decides any 6a
# replacement after reading the decision file.
# Refuses to run twice (an existing decision708-d708def.txt or day708-d708def.log = stop). Cap: CAP708 (default 202610270000).
# Terminal line, always: "<ts> END".
#   nohup bash scripts/detach.sh bash experiments/690-def-readiness-1015/wait708.sh > experiments/690-def-readiness-1015/work/wait708.nohup 2>&1 &
set -uo pipefail
REPO=~/Code/virtual-cell-challenge-2026; cd "$REPO" || exit 1
E=experiments/690-def-readiness-1015; W=$E/work; LOG="${LOG708:-$W/wait708.log}"
if [ -z "${WAIT708_FROZEN:-}" ]; then
  FZ="$W/.frozen-$$-wait708.sh"; cp "$REPO/$E/wait708.sh" "$FZ" || exit 1
  WAIT708_FROZEN=1 exec bash "$FZ" "$@"
fi
CAP="${CAP708:-202610270000}"
# 10-20 (SIXA Amendment 2): 710's board read is POSITIVE; if W1 (673, H1 x2.5 on R0's recipe) reads BELOW, re-arm with
# WAIT708_R710=HARM so decide708 itself blocks any H1 move above x2. Only HARM blocks; anything else does not.
R710V="${WAIT708_R710:-POSITIVE}"
# test knobs (production sets none): paths, poll / re-check seconds, the arm1015 R0 process pattern, the day708 script
MAN="${MAN708:-data/challenge-def/manifest.json}"; WLOG="${WLOG708:-experiments/691-release-watcher/work/watch691.log}"
SUB="${SUB708:-experiments/670-coord-1015/work/submit-R0.log}"; POLL="${POLL708:-60}"; RECHECK="${RECHECK708:-120}"
ARMPAT="${ARMPAT708:-arm1015\.sh R0 }"
X=experiments/708-def-panel-repricing; XW="${XW708:-$X/work}"; DAY="${DAY708:-$X/day708.sh}"
DEC=$XW/decision708-d708def.txt; DLOG=$XW/day708-d708def.log
log(){ printf '%s %s\n' "$(date -u +%FT%TZ)" "$*" >> "$LOG"; }
trap 'log "END"; rm -f "$W/.frozen-$$-wait708.sh"' EXIT
case "$R710V" in POSITIVE|HARM|NOT_HARM|NONE) ;; *) log "STOP: WAIT708_R710 must be POSITIVE|HARM|NOT_HARM|NONE (got $R710V)"; exit 2;; esac
capped(){ [ "$(date -u +%Y%m%d%H%M)" -ge "$CAP" ]; }
t_rel(){
  local t; t="$(/usr/bin/grep -E "RELEASE DETECTED" "$WLOG" 2>/dev/null | head -1 | cut -c1-20)"
  [ -n "$t" ] && { echo "$t"; return; }
  [ -f "$MAN" ] && date -u -r "$(stat -f %m "$MAN")" +%FT%TZ
}
arm_alive(){ ps -axo comm=,args= | awk -v p="$ARMPAT" '$1 ~ /bash$/ && $0 ~ p' | /usr/bin/grep -q .; }
# lines of submit-R0.log stamped after $1 (the release)
after(){ awk -v t="$1" 'substr($0, 1, 20) > t' "$SUB" 2>/dev/null; }
r0_done(){   # echoes the reason and returns 0 when R0 is done
  local t="$1" a; a="$(after "$t")"
  echo "$a" | /usr/bin/grep -q "SUBMITTED ok" && { echo "R0 SUBMITTED ok ($(echo "$a" | /usr/bin/grep -E 'entry id:' | tail -1 | cut -c1-80))"; return 0; }
  echo "$a" | tail -1 | /usr/bin/grep -qE " END$" && ! arm_alive && { echo "R0 terminal without a submit: $(echo "$a" | /usr/bin/grep -E 'NOT SUBMITTING|ERROR' | tail -1 | cut -c1-200)"; return 0; }
  return 1
}
echo $$ > "${LOG%.log}.pid"
log "wait708 armed (pid $$): cap $CAP; waits for a real release ($MAN) and R0 done ($SUB), then day708 d708def with R705=CLEARED R707=WEAK R710=$R710V"
[ -e "$DEC" ] || [ -e "$DLOG" ] && { log "STOP: $DEC or $DLOG already exists (day708 d708def ran); refusing to run twice"; exit 1; }
last=""
until [ -f "$MAN" ]; do
  [ "$last" != r ] && log "waiting: no release ($MAN absent)"; last=r
  capped && { log "STOP: no release by cap $CAP"; exit 0; }
  sleep "$POLL"
done
T="$(t_rel)"; log "release: $MAN present; T_REL $T"
last=""
while :; do
  capped && { log "STOP: R0 not done by cap $CAP (day708 not run)"; exit 0; }
  if why="$(r0_done "$T")"; then
    sleep "$RECHECK"
    if why2="$(r0_done "$T")"; then log "R0 done (checked twice, $RECHECK s apart): $why2"; break; fi
    log "R0 looked done but not on the re-check (a relaunch?): $why"
  else
    [ "$last" != w ] && log "waiting: R0 not done yet (last R0 line after T_REL: $(after "$T" | tail -1 | cut -c1-120))"; last=w
  fi
  sleep "$POLL"
done
[ -e "$DEC" ] || [ -e "$DLOG" ] && { log "STOP: $DEC or $DLOG appeared meanwhile (someone ran day708 d708def); not running it again"; exit 1; }
log "day708: R705=CLEARED R707=WEAK R710=$R710V bash $DAY data/challenge-def experiments/692-def-r0/work/r0 d708def"
t0=$(date +%s)
R705=CLEARED R707=WEAK R710="$R710V" bash "$DAY" data/challenge-def experiments/692-def-r0/work/r0 d708def > "$XW/day708-d708def.nohup" 2>&1
rc=$?
log "day708 rc $rc in $(( ($(date +%s) - t0) / 60 )) min; its tail: $(tail -2 "$DLOG" 2>/dev/null | tr '\n' ' ' | cut -c1-300)"
if [ -s "$DEC" ]; then
  log "decision file $DEC:"
  sed 's/^/  | /' "$DEC" >> "$LOG"
else
  log "no decision file (day708 did not reach its decisions; see $DLOG)"
fi
