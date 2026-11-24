#!/usr/bin/env bash
# 420: run comp420.py on the given folds, one fold per process, each behind the memory gate (>= 40 % free, >= 20 GiB
# free+inactive). Usage: bash chain420.sh <fold> [<fold> ...] [-- extra args for comp420.py]
set -uo pipefail
cd "$(dirname "$0")"
export PYTHONDONTWRITEBYTECODE=1
gate(){ local pct fi; pct="$(memory_pressure | awk '/free percentage/{gsub(/%/,"",$5); print $5}')"
  fi="$(vm_stat | awk '/Pages free/{f=$3}/Pages inactive/{i=$3}END{gsub(/\./,"",f);gsub(/\./,"",i);printf "%.1f",(f+i)*16384/2^30}')"
  echo "$(date -u +%FT%TZ) gate: free ${pct}% free+inactive ${fi} GiB" >> work/chain420.log
  [ "${pct:-0}" -ge 40 ] && awk -v v="$fi" 'BEGIN{exit !(v>=20)}'; }
for f in "$@"; do
  n=0; until gate; do n=$((n+1)); [ $n -ge 240 ] && { echo "gate never opened for $f" >> work/chain420.log; exit 3; }; sleep 30; done
  echo "$(date -u +%FT%TZ) start $f" >> work/chain420.log
  /usr/bin/time -l uv run python comp420.py "$f" > "work/comp420-$f.log" 2>&1
  echo "$(date -u +%FT%TZ) end $f rc=$?" >> work/chain420.log
  touch "work/comp420-$f.done"
done
