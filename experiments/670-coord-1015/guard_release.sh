#!/usr/bin/env bash
# 670 (coordinator [2d6f46], 2026-10-16): once the D/E/F release is detected, every later submit is a FINAL entry, so no A/B/C package may
# go out. On detection this writes NOGO-<label> for every A/B/C slot label 1016s1 … 1021s2 (arm1015.sh / armready.sh stop on it before
# their submit) and logs it. Detection: data/challenge-def/ has a manifest (691 unpacked it), OR watch691's log says the release was seen,
# OR `vcc --json datasets list` lists anything besides the controls (polled every 10 min). It never submits and never touches R0's waiter.
#   nohup bash scripts/detach.sh bash experiments/670-coord-1015/guard_release.sh > experiments/670-coord-1015/work/guard_release.nohup 2>&1 &
set -uo pipefail
REPO=~/Code/virtual-cell-challenge-2026; cd "$REPO" || exit 1
W=experiments/670-coord-1015/work; LOG=$W/guard_release.log
log(){ printf '%s %s\n' "$(date -u +%FT%TZ)" "$*" >> "$LOG"; }
echo $$ > $W/guard_release.pid; log "armed (pid $$): blocks A/B/C slots 1016s1..1021s2 on release; cap 202610270000"
released(){
  [ -s data/challenge-def/manifest.json ] && { echo "data/challenge-def/manifest.json present"; return 0; }
  /usr/bin/grep -qiE "RELEASE (DETECTED|SEEN)|new dataset id" experiments/691-release-watcher/work/watch691.log 2>/dev/null && { echo "watch691 log: release"; return 0; }
  L="$(vcc --json datasets list 2>/dev/null)"; [ -n "$L" ] || return 1
  N="$(printf '%s' "$L" | python3 -c 'import json,sys
try: d=json.load(sys.stdin)
except Exception: print(0); sys.exit()
items=d if isinstance(d,list) else d.get("datasets", d.get("items", []))
print(sum(1 for x in items if str(x.get("id", "")) != "controls" and x.get("available") is True))' 2>/dev/null)"
  # 10-20 18:58Z false fire: Arc LISTED controls-final with "available": false ("final_week_not_active") two days early. Only an
  # AVAILABLE non-controls dataset is a release (as release691.py's diff requires).
  [ "${N:-0}" -gt 0 ] && { echo "vcc datasets list: $N AVAILABLE dataset(s) besides "controls""; return 0; }
  return 1; }
until WHY="$(released)"; do
  [ "$(date -u +%Y%m%d%H%M)" -ge 202610270000 ] && { log "cap reached without a release"; log "END"; exit 0; }
  sleep 600; done
log "RELEASE detected: $WHY"
for d in 16 17 18 19 20 21; do for s in 1 2; do touch "$W/NOGO-10${d}s$s"; done; done
log "NOGO written for 1016s1..1021s2 (A/B/C slots blocked; every submit is FINAL now). R0's waiter is untouched."
log "END"
