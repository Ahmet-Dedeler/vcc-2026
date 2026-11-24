#!/usr/bin/env bash
# 492g PLOLO runner: waits for 04:30Z (the coordinator's <= 4 GB window ends), then runs arms492g.py per fold, serially, each
# started only when vm_stat free+inactive >= 26 GiB (peak ~5.4 GiB, ~15 s per fold), then the pricer. Capped at 12:00Z.
REPO=~/Code/virtual-cell-challenge-2026; cd "$REPO" || exit 1
X=experiments/492-per-target-sharpening; LOG="$X/work/run492g.log"
log(){ printf '%s %s\n' "$(date -u +%FT%TZ)" "$*" | tee -a "$LOG"; }
free_gib(){ vm_stat | awk '/page size of/ {ps=$8} /Pages free/ {f=$3} /Pages inactive/ {i=$3} END {gsub(/\./,"",f); gsub(/\./,"",i); printf "%d", (f+i)*ps/1073741824}'; }
until [ "$(date -u +%Y%m%d%H%M)" -ge 202609270430 ]; do sleep 30; done
log "start (free+inactive $(free_gib) GiB)"
for f in k562_gwps orion_hct116 orion_hek293t kolf2; do
  until [ "$(free_gib)" -ge 26 ]; do
    [ "$(date -u +%Y%m%d%H%M)" -ge 202609271200 ] && { log "memory gate never opened for $f by 12:00Z; stopping"; exit 1; }
    sleep 30
  done
  log "$f start (free+inactive $(free_gib) GiB)"
  /usr/bin/time -l uv run python $X/arms492g.py "$f" > "$X/work/arms492g-$f.log" 2>&1 || { log "$f FAILED"; exit 1; }
  log "$f done: $(grep -E 'ALL KNOWN|peak RSS' "$X/work/arms492g-$f.log" | tr '\n' ' ')"
done
uv run python $X/analyse492g.py > "$X/work/analyse492g.log" 2>&1 || { log "pricer FAILED"; exit 1; }
log "pricer done END"
