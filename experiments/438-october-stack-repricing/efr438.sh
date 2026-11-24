#!/usr/bin/env bash
# efr (no 3-donor operator): folds unlevelled -> fold-route level factors -> the folds' re-levelled twins (resume).
set -uo pipefail
cd "$(dirname "$0")"
export PYTHONDONTWRITEBYTECODE=1
./chain438b.sh comp:efr:k562_gwps comp:efr:orion_hct116 comp:efr:orion_hek293t
uv run python foldlevel438.py efr efr_s1 >> work/chain438.log 2>&1
echo "$(date -u +%FT%TZ) foldlevel438 efr rc=$?" >> work/chain438.log
./chain438b.sh compr:efr:k562_gwps:held:work/efr/level438.json compr:efr:orion_hct116:held:work/efr/level438.json compr:efr:orion_hek293t:held:work/efr/level438.json
