#!/bin/bash
# 693: the full v2 == v1 KA, ONE serial stream (~4 GB: the light budget beside a running build); "<ts> END" at the end
cd ~/Code/virtual-cell-challenge-2026 || exit 1
W=experiments/693-dnmag-speed/work
for job in "0 updn" "1 updn" "2 updn" "0 dn" "1 dn" "2 dn"; do
  set -- $job
  uv run python experiments/693-dnmag-speed/ka693.py "$1" "$2" > "$W/ka693_s$1_$2.log" 2>&1
  echo "$(date -u +%FT%TZ) s$1 $2 rc $? $(tail -1 "$W/ka693_s$1_$2.log" | cut -c1-170)" >> "$W/ka693_all.log"
done
echo "$(date -u +%FT%TZ) END" >> "$W/ka693_all.log"
