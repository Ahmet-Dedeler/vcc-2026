#!/bin/bash
# 418 second pass, strictly serial behind the memory gate: the A/B/C panel render (abc418.py), then the dose
# continuation (comp418.py <fold> --extra) on the four folds. Writes only into ./work.
cd "$(dirname "$0")" || exit 1
export NUMBA_CACHE_DIR="$PWD/work/numba-cache"
export PYTHONDONTWRITEBYTECODE=1   # imports of 414/415/218 modules must not drop .pyc files outside this dir
gate() {
  local pct fi
  pct=$(memory_pressure | awk '/free percentage/ {gsub("%","",$5); print $5}')
  fi=$(vm_stat | awk '/Pages free|Pages inactive/ {gsub("\\.","",$3); s+=$3} END {printf "%d", s*16384/2^30}')
  [ "${pct:-0}" -ge 40 ] && [ "${fi:-0}" -ge 20 ]
}
step() {  # name, then the python args
  local name=$1; shift
  [ -f "work/$name.done" ] && { echo "[chain] $name done, skipping"; return 0; }
  until gate; do echo "[chain] $(date -u +%H:%M:%SZ) memory gate closed before $name; waiting 60 s"; sleep 60; done
  echo "[chain] $(date -u +%H:%M:%SZ) start $name"
  /usr/bin/time -l uv run python "$@" > "work/$name.log" 2>&1
  local rc=$?
  echo "[chain] $(date -u +%H:%M:%SZ) $name exit $rc; $(grep 'maximum resident' "work/$name.log")"
  [ $rc -eq 0 ] && touch "work/$name.done"
  return $rc
}
step abc418 abc418.py || exit 1
for f in orion_hct116 k562_gwps orion_hek293t kolf2; do
  step "comp418x-$f" comp418.py "$f" --extra || exit 1
done
