#!/usr/bin/env bash
# 534's swap guard, run by the owner of 534 (not the worker).
# AGENTS §9 wants swap free >= 400 MiB before any heavy start. gate534.py (Part B's gate) does not read swap, and
# editing a live worker's gate is unsafe, so this guard acts through the owner's override file instead. Every 30 s
# until 2026-10-03T00:00Z it writes GATE-OVERRIDE.json:
#   - floor 27 GiB while swap free >= 400 MiB;
#   - floor 33 GiB (the preregistered value; unreachable with the desktop apps open, so effectively closed) below it.
# Writes are atomic (tmp + mv). Every change is logged to swapguard534.log. It only governs STARTS: a running fit is
# never interrupted.
set -u
F=~/Code/virtual-cell-challenge-2026/experiments/534-generalization-audit/work/partB
END=$(date -j -u -f "%Y-%m-%dT%H:%M:%SZ" "2026-10-03T00:00:00Z" +%s)
WHY27="PREREGISTRATION-B §9's 33 GiB free+inactive was unreachable all afternoon: 25.3-29.3 GiB with memory_pressure 78-82 % free and swap unused, because long-lived desktop apps (Aside renderers, ~10 GB) hold memory the 33 rule did not budget for. 27 GiB is above AGENTS §9's 20 GiB standard. It opens when worker 528 is between jobs (28-31 GiB in its logs) and stays shut while a 528 job is resident (~25-26 GiB), so the two heavy workers take turns (528 itself waits below 23 GiB). Unchanged in gate534.py: memory_pressure >= 40 % free, no builder/chain running, strictly serial, the 23:45Z / 00:30-06:00Z windows. Swap free >= 400 MiB (AGENTS §9) is enforced by the owner's swapguard534.sh through this file. Operational change only: no reading rule, design or seed moves. Record it as the owner's deviation in RESULTS-B.md."
state=""
while [ "$(date -u +%s)" -lt "$END" ]; do
  read -r num unit < <(sysctl -n vm.swapusage | sed -E 's/.*free = ([0-9.]+)([MG]).*/\1 \2/')
  free_m=$(awk -v n="$num" -v u="$unit" 'BEGIN{printf "%.0f", (u=="G") ? n*1024 : n}')
  if [ -n "$free_m" ] && [ "$free_m" -ge 400 ]; then
    x=27; why="$WHY27 [swapguard: swap free ${free_m} MiB >= 400]"
  else
    x=33; why="swapguard534.sh: swap free ${free_m:-unreadable} MiB < 400 (AGENTS §9), so the floor is back at the preregistered 33 GiB (effectively closed) until swap recovers. No new heavy start; a running fit continues."
  fi
  if [ "$x" != "$state" ]; then
    python3 - "$F" "$x" "$why" <<'PY'
import json, os, sys
f, x, why = sys.argv[1], float(sys.argv[2]), sys.argv[3]
tmp = os.path.join(f, ".GATE-OVERRIDE.json.tmp")
with open(tmp, "w") as h:
    json.dump({"min_free_gib": x, "by": "Public data scope [2e83f3] (owner of 534), via swapguard534.sh", "reason": why}, h, indent=2)
    h.write("\n")
os.replace(tmp, os.path.join(f, "GATE-OVERRIDE.json"))
PY
    echo "$(date -u +%FT%TZ) floor ${x} (swap free ${free_m:-?} MiB)" >> "$F/swapguard534.log"
    state=$x
  fi
  sleep 30
done
echo "$(date -u +%FT%TZ) END" >> "$F/swapguard534.log"
