#!/usr/bin/env bash
# 438 serial chain, v2 (chain438.sh is live in another waiter; never edit a running bash script). Jobs:
#   comp:<panel>:<fold>[:<fold_ops>[:<level>]]   compr:... the same with --resume
#   full:<panel>[:<tag>[:<op>]]                  fullr:... the same with --resume
set -uo pipefail
cd "$(dirname "$0")"
export PYTHONDONTWRITEBYTECODE=1
for job in "$@"; do
  n=0; until ./gate438.sh >> work/chain438.log; do n=$((n+1)); [ $n -ge 240 ] && { echo "gate never opened for $job" >> work/chain438.log; exit 3; }; sleep 30; done
  IFS=: read -r kind a b c d <<< "$job"
  extra=""; case "$kind" in compr|fullr) extra="--resume";; esac
  if [ "$kind" = comp ] || [ "$kind" = compr ]; then
    ops="${c:-held}"; lv="${d:-}"
    log="work/logs/comp438-$a-$b-$ops${lv:+-L}${extra:+-r}.log"
    echo "$(date -u +%FT%TZ) start $job" >> work/chain438.log
    /usr/bin/time -l uv run python comp438.py "$a" "$b" --fold-ops "$ops" ${lv:+--level "$lv"} $extra > "$log" 2>&1
    rc=$?
  else
    tag="${b:-$a}"; log="work/logs/full438-$tag${extra:+-r}.log"
    echo "$(date -u +%FT%TZ) start $job" >> work/chain438.log
    /usr/bin/time -l uv run python full438.py "$a" --tag "$tag" ${c:+--op "$c"} $extra > "$log" 2>&1
    rc=$?
  fi
  echo "$(date -u +%FT%TZ) end $job rc=$rc $(grep 'maximum resident' "$log" | awk '{printf "maxrss %.2f GB", $1/1e9}')" >> work/chain438.log
done
