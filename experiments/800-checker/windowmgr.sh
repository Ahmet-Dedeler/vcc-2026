#!/bin/bash
# windowmgr.sh — keep the lab inside charter Amendment 1's reserved windows without manual pausing.
# Every 60 s: if window800.sh prints WINDOW and more than ONE lab worker is alive, stop the extras (lowest priority first;
# the lead session is the second agent) and record them in resume/paused.txt; when the window clears, resume each paused
# worker by session with resume/after-window.md. Capped: exits at the CAP time (default 2026-11-06T00:00Z) or on a STOP file.
#   nohup bash scripts/detach.sh bash experiments/800-checker/windowmgr.sh > experiments/800-checker/resume/windowmgr.log 2>&1 &
R=~/Code/virtual-cell-challenge-2026; cd "$R" || exit 1
E=experiments/800-checker; P=$E/resume/paused.txt; STOP=$E/resume/STOP-windowmgr
PRIORITY="803-de-engine 810-strong-direction 812-de-signs 811-lolo-model 813-public-data"   # kept first -> last
CAP=$(date -j -u -f "%Y-%m-%dT%H:%M:%SZ" "${WINDOWMGR_CAP:-2026-11-06T00:00:00Z}" +%s)
log(){ echo "$(date -u +%FT%TZ) $*"; }
# a worker is alive only if its recorded pid is a live process NAMED claude: pids are recycled (10-17: 805's stale pid 9077 had
# become an unrelated `sleep 30`), and this script TERMs whole process groups, so a bare kill -0 could hit a stranger's group
alive(){ local p="experiments/$1/claude-run.pid"; [ -f "$p" ] && [ "$(ps -o comm= -p "$(cat "$p")" 2>/dev/null | tr -d ' ')" = "claude" ]; }
touch "$P"
log "start (priority: $PRIORITY)"
while [ ! -f "$STOP" ] && [ "$(date -u +%s)" -lt "$CAP" ]; do
  if ! W="$(bash $E/window800.sh)"; then
    kept=0
    for n in $PRIORITY $(ls -d experiments/8[0-9][0-9]-* 2>/dev/null | xargs -n1 basename | grep -v -w -e 803-de-engine -e 810-strong-direction -e 812-de-signs -e 811-lolo-model -e 813-public-data); do
      alive "$n" || continue
      if [ "$kept" -lt 1 ]; then kept=1; continue; fi
      pid=$(cat "experiments/$n/claude-run.pid"); kill -TERM -- "-$pid" 2>/dev/null || kill -TERM "$pid" 2>/dev/null
      grep -qx "$n" "$P" || echo "$n" >> "$P"
      log "WINDOW ($W): paused $n (pid $pid)"
    done
  elif [ -s "$P" ]; then
    while read -r n; do
      [ -n "$n" ] || continue
      if alive "$n"; then log "$n already alive"; continue; fi
      bash $E/lab_subagent.sh --resume "$n" $E/resume/after-window.md && log "resumed $n"
    done < "$P"
    : > "$P"
  fi
  sleep 60
done
log "END"
