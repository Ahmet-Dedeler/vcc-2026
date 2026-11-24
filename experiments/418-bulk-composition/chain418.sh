#!/bin/bash
# 418: run the remaining folds strictly one at a time, each behind the memory gate
# (memory_pressure >= 40 % free AND vm_stat free+inactive >= 20 GiB). Writes only into ./work.
cd "$(dirname "$0")" || exit 1
export NUMBA_CACHE_DIR="$PWD/work/numba-cache"
export PYTHONDONTWRITEBYTECODE=1   # imports of 414/415/218 modules must not drop .pyc files outside this dir
gate() {
  local pct fi
  pct=$(memory_pressure | awk '/free percentage/ {gsub("%","",$5); print $5}')
  fi=$(vm_stat | awk '/Pages free|Pages inactive/ {gsub("\\.","",$3); s+=$3} END {printf "%d", s*16384/2^30}')
  [ "${pct:-0}" -ge 40 ] && [ "${fi:-0}" -ge 20 ]
}
for f in "$@"; do
  [ -f "work/comp418-$f.done" ] && { echo "[chain] $f done, skipping"; continue; }
  until gate; do echo "[chain] $(date -u +%H:%M:%SZ) memory gate closed before $f; waiting 60 s"; sleep 60; done
  echo "[chain] $(date -u +%H:%M:%SZ) start $f"
  /usr/bin/time -l uv run python comp418.py "$f" > "work/comp418-$f.log" 2>&1
  rc=$?
  echo "[chain] $(date -u +%H:%M:%SZ) $f exit $rc; $(grep 'maximum resident' "work/comp418-$f.log")"
  [ $rc -eq 0 ] && touch "work/comp418-$f.done" || exit $rc
done
