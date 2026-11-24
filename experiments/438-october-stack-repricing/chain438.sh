#!/usr/bin/env bash
# 438: run jobs strictly serially, each behind the memory gate. Each arg is one job: "comp:<panel>:<fold>[:<fold_ops>[:<level>]]"
# or "full:<panel>[:<tag>[:<op>]]". Logs to work/logs/, a line per job to work/chain438.log.
set -uo pipefail
cd "$(dirname "$0")"
export PYTHONDONTWRITEBYTECODE=1
for job in "$@"; do
  n=0; until ./gate438.sh >> work/chain438.log; do n=$((n+1)); [ $n -ge 240 ] && { echo "gate never opened for $job" >> work/chain438.log; exit 3; }; sleep 30; done
  IFS=: read -r kind a b c d <<< "$job"
  if [ "$kind" = comp ]; then
    ops="${c:-held}"; lv="${d:-}"
    log="work/logs/comp438-$a-$b-$ops${lv:+-L}.log"
    echo "$(date -u +%FT%TZ) start $job" >> work/chain438.log
    /usr/bin/time -l uv run python comp438.py "$a" "$b" --fold-ops "$ops" ${lv:+--level "$lv"} > "$log" 2>&1
    rc=$?
  else
    tag="${b:-$a}"; log="work/logs/full438-$tag.log"
    echo "$(date -u +%FT%TZ) start $job" >> work/chain438.log
    /usr/bin/time -l uv run python full438.py "$a" --tag "$tag" ${c:+--op "$c"} > "$log" 2>&1
    rc=$?
  fi
  echo "$(date -u +%FT%TZ) end $job rc=$rc $(grep 'maximum resident' "$log" | awk '{printf "maxrss %.2f GB", $1/1e9}')" >> work/chain438.log
done
