#!/usr/bin/env bash
# 420: per-context re-levelling variants on PLOLO (comp420.py --level), one fold per process, behind the memory gate.
# Usage: bash chain420L.sh <fold> [<fold> ...]   (specs read from work/level_specs.txt, written by analyse420's factors)
set -uo pipefail
cd "$(dirname "$0")"
export PYTHONDONTWRITEBYTECODE=1
eval "SPECS=($(cat work/level_specs.txt))"
gate(){ local pct fi; pct="$(memory_pressure | awk '/free percentage/{gsub(/%/,"",$5); print $5}')"
  fi="$(vm_stat | awk '/Pages free/{f=$3}/Pages inactive/{i=$3}END{gsub(/\./,"",f);gsub(/\./,"",i);printf "%.1f",(f+i)*16384/2^30}')"
  echo "$(date -u +%FT%TZ) gate: free ${pct}% free+inactive ${fi} GiB" >> work/chain420L.log
  [ "${pct:-0}" -ge 40 ] && awk -v v="$fi" 'BEGIN{exit !(v>=20)}'; }
for f in "$@"; do
  n=0; until gate; do n=$((n+1)); [ $n -ge 240 ] && { echo "gate never opened for $f" >> work/chain420L.log; exit 3; }; sleep 30; done
  echo "$(date -u +%FT%TZ) start $f (${#SPECS[@]} specs)" >> work/chain420L.log
  /usr/bin/time -l uv run python comp420.py "$f" --level "${SPECS[@]}" > "work/comp420L-$f.log" 2>&1
  echo "$(date -u +%FT%TZ) end $f rc=$?" >> work/chain420L.log
done
