#!/usr/bin/env bash
# 438's memory gate (AGENTS 9): memory_pressure >= 40 % free AND vm_stat free+inactive >= 20 GiB. Exit 0 = open.
pct="$(memory_pressure | awk '/free percentage/{gsub(/%/,"",$5); print $5}')"
fi="$(vm_stat | awk '/Pages free/{f=$3}/Pages inactive/{i=$3}END{gsub(/\./,"",f);gsub(/\./,"",i);printf "%.1f",(f+i)*16384/2^30}')"
echo "$(date -u +%FT%TZ) gate: free ${pct}% free+inactive ${fi} GiB"
[ "${pct:-0}" -ge 40 ] && awk -v v="$fi" 'BEGIN{exit !(v>=20)}'
