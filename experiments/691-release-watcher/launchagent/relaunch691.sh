#!/bin/bash
# 691 LaunchAgent payload (DRAFT; installed only after Ahmet's explicit yes, relayed in his words by the coordinator).
# launchd runs this at login (RunAtLoad) and every 30 min (StartInterval). Between 10-18 18:00Z and 10-27 00:00Z it
# re-launches the D/E/F release watcher and the R0 adapter when a reboot or a crash killed them, with the commands recorded
# in docs/PLAN-2026-10-15.md. Outside that window it exits at once. It re-launches a process only if
#   - no copy of it is running, AND
#   - its log exists and its last line is NOT "END" (it was killed mid-run; a run that ended on purpose (a STOP, the cap,
#     an ambiguous release) waits for a person and is never restarted here).
# Both scripts resume after a restart (691 RESULTS v2). Every decision is appended to launchagent.log.
# DRY691=1 prints the decisions without launching anything (tests).
REPO=~/Code/virtual-cell-challenge-2026
cd "$REPO" || exit 0
LA="${LA691:-experiments/691-release-watcher/work/launchagent.log}"   # *691 overrides: tests only
NOW="${NOW691:-$(date -u +%Y%m%d%H%M)}"
FROM="${FROM691:-202610181800}"; UNTIL="${UNTIL691:-202610270000}"
[ "$NOW" -ge "$FROM" ] && [ "$NOW" -le "$UNTIL" ] || exit 0
say(){ printf '%s %s\n' "$(date -u +%FT%TZ)" "$*" >> "$LA"; [ "${DRY691:-0}" = 1 ] && echo "$*"; }
running(){ ps -axo comm=,args= | awk -v s="$1" '$1 ~ /bash$/ && index($0, s) > 0' | grep -q .; }
ended(){ [ -s "$1" ] && tail -1 "$1" | grep -qE " END$"; }
maybe(){  # $1 label, $2 script path, $3 log, $4 command (running = any bash whose argv holds the script's basename:
          # the scripts exec a frozen copy named .frozen-<pid>-<basename>)
  if running "${2##*/}"; then return 0; fi
  if [ ! -s "$3" ]; then say "$1: no log at $3 (never armed): not launching"; return 0; fi
  if ended "$3"; then return 0; fi
  if [ "${DRY691:-0}" = 1 ]; then say "$1: WOULD re-launch: $4"; return 0; fi
  say "$1: not running and its log does not end in END: re-launching"
  /bin/bash -c "$4" >> "$LA" 2>&1
}
maybe watch691 "${S691:-experiments/691-release-watcher/watch691.sh}" "${L691:-experiments/691-release-watcher/work/watch691.log}" \
  "AUTO_R0=1 AUTO_ALLOW_DEPTH691=1 AUTO_ALLOW_WIDE691=1 nohup bash scripts/detach.sh bash experiments/691-release-watcher/watch691.sh >> experiments/691-release-watcher/work/watch691.nohup 2>&1 &"
maybe r0ready692 "${S692:-experiments/692-def-r0/r0ready692.sh}" "${L692:-experiments/692-def-r0/work/r0def.run.log}" \
  "nohup bash scripts/detach.sh bash experiments/692-def-r0/r0ready692.sh data/challenge-def experiments/692-def-r0/work/r0 r0def 202610270000 >> experiments/692-def-r0/work/r0ready692.nohup 2>&1 &"
exit 0
