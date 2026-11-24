#!/bin/zsh
cd "$(dirname "$0")"
for panel in real pess efp efr abc; do
  for held in k562_gwps orion_hct116 orion_hek293t; do
    uv run python gate250.py --wait
    echo "=== cross $panel/$held  $(date -u +%H:%M:%SZ) ==="
    /usr/bin/time -l uv run python cross250.py --panel $panel --fold $held > work/logs/cross-$panel-$held.log 2>&1 || { echo "!!! cross $panel/$held failed"; exit 1; }
    grep -E "KA |cross done|maximum resident" work/logs/cross-$panel-$held.log
  done
done
echo "=== cross done $(date -u +%H:%M:%SZ) ==="
