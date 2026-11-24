#!/usr/bin/env bash
# Lead-session monitor: one stdout line per event. Bash (not zsh) so globs iterate.
# Portable to macOS bash 3.2: no associative arrays (a name containing a variable-like word such
# as "free" was being evaluated arithmetically as an array index -- two false FINISHED events).
#   FINISHED <exp>   a supervised experiment wrote RESULTS.md
#   DEAD <exp>       its process is gone with no RESULTS.md (agent_supervisor may relaunch)
#   ARMED-SUBMIT     new lines in the armed-submission log
#   SLOT             30 min before / at the 00:00Z reset
#   CLOUD            the 221 bundle instance changed state
#   MEM              system free memory under 25%
#   HEARTBEAT        every 30 min: what is running, free memory
cd ~/Code/virtual-cell-challenge-2026
SEEN=work/lead_monitor.seen; : > "$SEEN"; first=1; lasthb=0; lastcloud=""; cloudt=0; sublines=0; slot30=0; slot0=0
mark(){ grep -qxF -- "$1" "$SEEN" && return 1; echo "$1" >> "$SEEN"; return 0; }
SUBLOGS="experiments/237-gate-validation/submit-armed.log experiments/249-cd4-deflate/submit-armed.log"
while :; do
  now=$(date -u +%s); hm=$(date -u +%H%M)
  running=""
  for d in experiments/2[0-9][0-9]-*/; do
    n=$(basename "$d"); [ -f "$d/claude-run.pid" ] || continue
    pid=$(cat "$d/claude-run.pid" 2>/dev/null)
    if [ -f "$d/RESULTS.md" ]; then
      if mark "done:$n"; then [ $first -eq 0 ] && echo "FINISHED $n (RESULTS.md written)"; fi
    elif kill -0 "$pid" 2>/dev/null; then running="$running $n"
    elif mark "dead:$n"; then [ $first -eq 0 ] && echo "DEAD $n pid=$pid gone, no RESULTS.md (supervisor loop may relaunch)"; fi
  done
  for SUBLOG in $SUBLOGS; do
    [ -f "$SUBLOG" ] || continue
    c=$(wc -l < "$SUBLOG"); key="$(echo "$SUBLOG" | tr '/' '_')"; prev=$(grep "^lines:$key=" "$SEEN" | tail -1 | cut -d= -f2); prev=${prev:-0}
    if [ "$c" -gt "$prev" ]; then [ $first -eq 0 ] && { echo "ARMED-SUBMIT ($SUBLOG):"; tail -n $((c-prev)) "$SUBLOG" | cut -c1-300; }; echo "lines:$key=$c" >> "$SEEN"; fi
  done
  if [ "$hm" -ge 2330 ] && [ $slot30 -eq 0 ]; then slot30=1; echo "SLOT 30 minutes to the 00:00Z reset; armed pids: $(cat experiments/237-gate-validation/submit-armed.log.pid 2>/dev/null) / $(cat experiments/249-cd4-deflate/submit-armed.log.pid 2>/dev/null)"; fi
  if [ "$hm" -lt 0100 ] && [ "$hm" -ge 0000 ] && [ $slot0 -eq 0 ] && [ $slot30 -eq 1 ]; then slot0=1; echo "SLOT 00:00Z reset passed"; fi
  free=$(memory_pressure 2>/dev/null | grep -o 'free percentage: [0-9]*' | grep -o '[0-9]*$'); free=${free:-0}
  if [ "$free" -lt 25 ]; then echo "MEM free=${free}% running=[${running# }]"; fi
  if [ $((now - cloudt)) -ge 900 ]; then cloudt=$now
    st=$(aws ec2 describe-instances --region us-east-1 --instance-ids ${CLOUD_INSTANCE_ID:-} --query 'Reservations[].Instances[].State.Name' --output text 2>/dev/null || echo "query-failed")
    if [ "$st" != "$lastcloud" ]; then [ $first -eq 0 ] && echo "CLOUD bundle instance ${CLOUD_INSTANCE_ID:-} state: ${lastcloud:-?} -> $st"; lastcloud=$st; fi
  fi
  if [ $((now - lasthb)) -ge 1800 ]; then lasthb=$now; echo "HEARTBEAT $(date -u +%H:%M)Z running=[${running# }] free=${free}% cloud=$lastcloud"; fi
  first=0; sleep 60
done
