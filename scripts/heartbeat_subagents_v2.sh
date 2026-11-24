#!/usr/bin/env bash
# Heartbeat + completion watcher for headless subagents: one line when any recent claude-run exits, one status line every 30 min
# with the board top-3 from the newest data/leaderboard snapshot. Arm it from a Claude session as a persistent Monitor:
#   Monitor(command="bash scripts/heartbeat_subagents.sh", persistent=true)
# Why (2026-09-06): a session that only waits is not re-invoked; Monitor events are. Written by the coordinator session.
cd ~/Code/virtual-cell-challenge-2026 || exit 1
declare -A seen
last_beat=0
while true; do
  now=$(date +%s)
  while IFS= read -r pidf; do
    [ -n "$pidf" ] || continue
    d=$(dirname "$pidf"); name=$(basename "$d"); pid=$(cat "$pidf" 2>/dev/null)
    [ -n "$pid" ] || continue
    if ! kill -0 "$pid" 2>/dev/null; then
      if [ -z "${seen[$name]:-}" ]; then
        seen[$name]=1
        res=$(grep -c '"type":"result"' "$d/claude-run.jsonl" 2>/dev/null || echo 0)
        hr=no; [ -f "$d/RESULTS.md" ] && hr=yes
        err=$(tail -c 240 "$d/claude-run.err" 2>/dev/null | tr '\n' ' ')
        echo "DONE $name pid=$pid result_events=$res RESULTS.md=$hr err=[$err]"
      fi
    fi
  done < <(find experiments -maxdepth 2 -name claude-run.pid -mmin -1440 2>/dev/null)
  if [ $((now - last_beat)) -ge 1800 ]; then
    last_beat=$now
    running=""
    while IFS= read -r pidf; do
      [ -n "$pidf" ] || continue
      pid=$(cat "$pidf" 2>/dev/null); [ -n "$pid" ] && kill -0 "$pid" 2>/dev/null && running="$running $(basename "$(dirname "$pidf")")"
    done < <(find experiments -maxdepth 2 -name claude-run.pid -mmin -1440 2>/dev/null)
    lb=$(ls -t data/leaderboard/*.json 2>/dev/null | head -1)
    top=$(python3 - "$lb" <<'PYX' 2>/dev/null
import json,sys
try:
    d=json.load(open(sys.argv[1]))
    def lists(o):
        if isinstance(o,list) and o and isinstance(o[0],dict): yield o
        if isinstance(o,dict):
            for v in o.values(): yield from lists(v)
    best=None
    for L in lists(d):
        keys=[k for k in ('score_avg','score','overall','overall_score','scaled') if k in L[0]]
        if keys:
            k=keys[0]; S=sorted(L,key=lambda x:-(float(x.get(k) or 0))); best=S[:3]
            nk=next((n for n in ('team_name','team','name','username') if n in best[0]),None)
            ours=next((i+1 for i,x in enumerate(S) if str(x.get(nk))=='Ahmet Dedeler'),None)
            print(' | '.join(f"{(x.get(nk) if nk else '?')}:{float(x.get(k) or 0):.5f}" for x in best)+f" | our rank {ours}/{len(S)}"); break
    else: print('no score list')
except Exception as e: print('unreadable',e)
PYX
)
    echo "HEARTBEAT $(date -u +%H:%MZ) running:[${running:-none}] board_top3: ${top} ($(basename "${lb:-none}"))"
  fi
  sleep 60
done
