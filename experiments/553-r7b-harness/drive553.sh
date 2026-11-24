#!/bin/bash
# 553: serial, gated driver for this directory's jobs. Queue: work/queue553.txt, re-read before every job:
#   <label> <min_free_gib> <rss_cap_gb> <abort|cont> <script> [args...]      (# and blank lines skipped)
# <script> is a file in experiments/553-r7b-harness; it runs as `/usr/bin/time -l uv run python <that file> args...`.
# A label is skipped once work/done553/<label>.ok or .fail exists. min_free_gib 0 = a light job (< 2 GB, no gate).
# Heavy-job gate (mandate): no heavy job 2026-10-04 05:00Z -> 2026-10-05 18:00Z; vm_stat free+inactive >= min_free_gib
# (24, or 27 for the operator fit); memory_pressure free >= 40 %; def_pipeline.other_builders() == 0; no bash running
# run_arm_generic.sh / def_day.sh / build5NN*.sh / rebuild5NN*.sh; no run537.py python alive (worker 537; never signalled).
# Watchdog: the job's process tree RSS is sampled every 5 s; above <rss_cap_gb> it is killed (marked .fail).
# A failed `abort` job (a known answer) writes work/KA-FAIL-553 and stops the driver.
# Stop: touch work/STOP553. Cap: 2026-10-06 12:00Z.
REPO=~/Code/virtual-cell-challenge-2026
HERE="$REPO/experiments/553-r7b-harness"
cd "$REPO" || exit 1
mkdir -p "$HERE/work/logs" "$HERE/work/done553"
LOG="$HERE/work/drive553.log"; Q="$HERE/work/queue553.txt"
echo $$ > "$HERE/work/drive553.pid"
log(){ printf '%s %s\n' "$(date -u +%FT%TZ)" "$*" >> "$LOG"; }
log "armed (pid $$)"
chains_alive(){
  # interpreter + script path from argv (AGENTS 6.2): argv[0] is a shell and its first non-flag argument is the chain
  # script -- never a substring anywhere in the args (headless workers carry these names in their prompts)
  ps -axo args= | awk '{n0=$1; sub(/.*\//, "", n0); if (n0 !~ /^(bash|sh|zsh)$/) next;
    s=""; for (i=2; i<=NF; i++) if ($i !~ /^-/) {s=$i; break}; sub(/.*\//, "", s);
    if (s ~ /^(run_arm_generic|def_day)[.]sh$/ || s ~ /^(re)?build5[0-9][0-9][^\/]*[.]sh$/) c++} END {exit c ? 0 : 1}'
}
w537_alive(){ ps -axo args= | awk '$1 ~ /python[0-9.]*$/ && $2 ~ /537-arm-pipeline\/run537[.]py$/ {c++} END {exit c ? 0 : 1}'; }
last_why=""; last_t=0
override_vals(){
  # the owner's answer to BLOCKED.md: GATE-OVERRIDE.json's operator_fit floors, printed "fi pct swap" (empty if absent)
  [ -f "$HERE/GATE-OVERRIDE.json" ] || return 0
  uv run python -c 'import json,sys; o=json.load(open(sys.argv[1]))["operator_fit"]; print(int(o["min_free_inactive_gib"]), int(o["min_pressure_pct"]), int(o["min_swap_free_mib"]))' "$HERE/GATE-OVERRIDE.json" 2>/dev/null
}
gate(){
  local need=$1 label=$2 minp=40 minsw=0
  if [ "$label" = "r7b_operator" ]; then
    ov=$(override_vals)
    if [ -n "$ov" ]; then set -- $ov; need=$1; minp=$2; minsw=$3; log "operator-fit gate from GATE-OVERRIDE.json: fi>=${need} GiB pressure>=${minp}% swapfree>=${minsw} MiB"; fi
  fi
  while true; do
    [ -e "$HERE/work/STOP553" ] && { log "STOP553 present"; log "END"; exit 0; }
    now=$(date -u +%Y%m%d%H%M)
    [ "$now" -ge 202610061200 ] && { log "past the 10-06 12:00Z cap"; log "END"; exit 0; }
    why=""
    if [ "$now" -ge 202610040500 ] && [ "$now" -lt 202610051800 ]; then why="$why blackout-10-04/05"; fi
    chains_alive && why="$why chain-running"
    w537_alive && why="$why run537-alive"
    nb=$(uv run python -c 'import sys; sys.path.insert(0,"scripts"); import def_pipeline as P; print(len(P.other_builders()))' 2>/dev/null)
    [ "${nb:-9}" = "0" ] || why="$why builders${nb}"
    fp=$(memory_pressure | awk -F': ' '/free percentage/ {gsub("%","",$2); print $2}')
    fi=$(vm_stat | awk '/Pages free|Pages inactive/ {s+=$NF} END {printf "%d", s*16384/2^30}')
    sw=$(sysctl -n vm.swapusage | awk '{for(i=1;i<=NF;i++) if($i=="free") {v=$(i+2); u=substr(v,length(v)); x=substr(v,1,length(v)-1); if(u=="G") x*=1024; printf "%d", x}}')
    [ "${fp:-0}" -ge "$minp" ] || why="$why free${fp}%<${minp}%"
    [ "${fi:-0}" -ge "$need" ] || why="$why fi${fi}<${need}GiB"
    [ "${sw:-0}" -ge "$minsw" ] || why="$why swapfree${sw}<${minsw}MiB"
    if [ -z "$why" ]; then log "gate ok: free ${fp}% free+inactive ${fi} GiB swap free ${sw} MiB builders 0 (need ${need} GiB, ${minp}%, ${minsw} MiB)"; return 0; fi
    t=$(date +%s)
    if [ "$why" != "$last_why" ] || [ $((t - last_t)) -ge 1800 ]; then log "gate closed:$why"; last_why="$why"; last_t=$t; fi
    sleep 120
  done
}
tree_rss_kb(){
  # sum of RSS (KiB) over the process tree rooted at $1
  ps -axo pid=,ppid=,rss= | awk -v root="$1" '{pp[$1]=$2; rs[$1]=$3} END {s=0; for (p in pp) {q=p; while (q != "" && q != 0 && q != 1) {if (q == root) {s+=rs[p]; break}; q=pp[q]}} print s}'
}
next_line(){
  [ -f "$Q" ] || return 1
  while read -r label gib cap mode script rest; do
    case "$label" in ""|\#*) continue;; esac
    [ -e "$HERE/work/done553/$label.ok" ] && continue
    [ -e "$HERE/work/done553/$label.fail" ] && continue
    echo "$label $gib $cap $mode $script $rest"; return 0
  done < "$Q"
  return 1
}
while true; do
  [ -e "$HERE/work/STOP553" ] && { log "STOP553 present"; log "END"; exit 0; }
  [ "$(date -u +%Y%m%d%H%M)" -ge 202610061200 ] && { log "past the cap"; log "END"; exit 0; }
  line=$(next_line) || { sleep 300; continue; }
  set -- $line; label=$1; gib=$2; cap=$3; mode=$4; script=$5; shift 5; args="$*"
  [ "$gib" -gt 0 ] && gate "$gib" "$label"
  line2=$(next_line) || continue
  [ "$line2" = "$line" ] || continue
  log "START $label ($script $args)"
  capkb=$(awk -v g="$cap" 'BEGIN {printf "%d", g*1e9/1024}')
  /usr/bin/time -l uv run python "$HERE/$script" $args > "$HERE/work/logs/$label.log" 2>&1 &
  JP=$!
  peak=0; killed=0
  while kill -0 "$JP" 2>/dev/null; do
    r=$(tree_rss_kb "$JP"); [ "${r:-0}" -gt "$peak" ] && peak=$r
    if [ "${r:-0}" -gt "$capkb" ]; then log "WATCHDOG: $label tree rss ${r} KiB > cap ${cap} GB: killing"; pkill -P "$JP"; kill "$JP"; killed=1; fi
    sleep 5
  done
  wait "$JP"; rc=$?
  [ $killed -eq 1 ] && rc=137
  tpk=$(/usr/bin/grep 'maximum resident' "$HERE/work/logs/$label.log" | tail -1 | awk '{printf "%.2f GB", $1/1e9}')
  real=$(/usr/bin/grep ' real ' "$HERE/work/logs/$label.log" | tail -1 | awk '{print $1 " s"}')
  log "END $label rc=$rc peakRSS(time -l)=$tpk tree-peak=$(awk -v k="$peak" 'BEGIN{printf "%.2f GB", k*1024/1e9}') real=$real"
  if [ $rc -eq 0 ]; then touch "$HERE/work/done553/$label.ok"; else
    touch "$HERE/work/done553/$label.fail"
    if [ "$mode" = "abort" ]; then log "KNOWN ANSWER / REQUIRED STEP FAILED: $label; stopping"; touch "$HERE/work/KA-FAIL-553"; log "END"; exit 1; fi
  fi
done
