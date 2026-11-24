#!/bin/bash
# 693 chain: after the H1 KA stream (ka693_all.sh) EXITS, (1) test693p --real-ka on the six real records, (2) abc693.py on
# real A/B/C blocks (<= ~4 GB, needs free+inactive >= 12 GiB), then a summary. Terminal line "<ts> END".
cd ~/Code/virtual-cell-challenge-2026 || exit 1
W=experiments/693-dnmag-speed/work; L=$W/chain693.log
log(){ printf '%s %s\n' "$(date -u +%FT%TZ)" "$*" >> "$L"; }
trap 'log END' EXIT
log "chain693 armed (pid $$)"
until /usr/bin/grep -q " END$" "$W/ka693_all.log" 2>/dev/null && ! ps -axo args= | /usr/bin/grep -v grep | /usr/bin/grep -q "bash experiments/693-dnmag-speed/work/ka693_all.sh"; do sleep 60; done
log "H1 KA done: $(/usr/bin/grep -cE 'rc 0 PASS' "$W/ka693_all.log") of 6 PASS"
/usr/bin/grep -E "rc [0-9]" "$W/ka693_all.log" | while IFS= read -r l; do log "  $l"; done
uv run python experiments/693-dnmag-speed/test693p.py --real-ka > "$W/test693p_real.log" 2>&1
log "test693p --real-ka: $(tail -1 "$W/test693p_real.log")"
# coordinator's guard: start only with >= 70 GiB free disk (two ~12 GB restores) and >= 12 GiB free+inactive; cap 10-15 18:00Z
until [ "$(df -g /System/Volumes/Data | awk 'NR==2{print $4}')" -ge 70 ] && \
      [ "$(uv run python -c "import sys; sys.path.insert(0,'scripts'); import def_pipeline as d; print(int(d.free_gib()))" 2>/dev/null || echo 0)" -ge 12 ]; do
  [ "$(date -u +%Y%m%d%H%M)" -ge 202610151800 ] && { log "abc693 NOT RUN: disk/memory guard not met by 10-15 18:00Z"; exit 3; }
  sleep 60; done
log "abc693 starting: disk $(df -g /System/Volumes/Data | awk 'NR==2{print $4}') GiB free"
/usr/bin/time -l uv run python experiments/693-dnmag-speed/abc693.py 3 > "$W/abc693.log" 2>&1
log "abc693 rc $?: $(/usr/bin/grep -E '^(PASS|FAIL)' "$W/abc693.log" | tail -1 | cut -c1-200); $(/usr/bin/grep 'maximum resident' "$W/abc693.log" | tr -s ' ')"
for f in "$W/pred_630arm.h5ad" "$W/pred_dn631.h5ad"; do
  [ -f "$f" ] && { log "deleting restore $f sha256 $(shasum -a 256 "$f" | cut -d' ' -f1) ($(stat -f %z "$f") bytes; regenerable from its .vcc)"; rm -f "$f"; }
done
log "disk after: $(df -g /System/Volumes/Data | awk 'NR==2{print $4}') GiB free"
