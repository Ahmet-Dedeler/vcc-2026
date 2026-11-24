#!/bin/zsh
# 250 — strictly serial driver: memory gate before every fold load, /usr/bin/time -l on every
# heavy script, logs under work/logs/.
cd "$(dirname "$0")"
mkdir -p work/logs
sweep () {
  local panel=$1 held=$2
  uv run python gate250.py --wait
  echo "=== sweep $panel/$held  $(date -u +%H:%M:%SZ) ==="
  /usr/bin/time -l uv run python sweep250.py --panel $panel --fold $held > work/logs/sweep-$panel-$held.log 2>&1
  local rc=$?
  grep -E "KA |arms in|FAIL|Error|maximum resident" work/logs/sweep-$panel-$held.log | tail -8
  [[ $rc -ne 0 ]] && echo "!!! sweep $panel/$held exited $rc" && return $rc
  return 0
}
bridge () {
  local held=$1
  uv run python gate250.py --wait
  echo "=== bridge $held  $(date -u +%H:%M:%SZ) ==="
  /usr/bin/time -l uv run python bridge250.py --fold $held > work/logs/bridge-$held.log 2>&1
  local rc=$?
  grep -E "KA |bridge done|FAIL|Error|maximum resident" work/logs/bridge-$held.log | tail -14
  [[ $rc -ne 0 ]] && echo "!!! bridge $held exited $rc" && return $rc
  return 0
}
set -o pipefail
for held in k562_gwps orion_hct116 orion_hek293t; do sweep abc $held || exit 1; done
for panel in real pess efp efr; do
  for held in k562_gwps orion_hct116 orion_hek293t; do sweep $panel $held || exit 1; done
done
sweep abc kolf2 || echo "kolf2 on abc failed (reported apart)"
for held in k562_gwps orion_hct116 orion_hek293t kolf2; do bridge $held || exit 1; done
uv run python gate250.py --wait
echo "=== basis  $(date -u +%H:%M:%SZ) ==="
/usr/bin/time -l uv run python basis250.py > work/logs/basis.log 2>&1
grep -E "KA-|maximum resident|done in" work/logs/basis.log
echo "=== run250 done $(date -u +%H:%M:%SZ) ==="
