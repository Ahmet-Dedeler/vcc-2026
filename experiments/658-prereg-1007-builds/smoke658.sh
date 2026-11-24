#!/usr/bin/env bash
# 658: bulkid658 smoke test on a REAL package, after C10b. Restores 497dr's .vcc (streamed zstd, low memory) and compares
# it with 497dr's built h5ad: bulk sums and --exact. Expected: BULK IDENTICAL. The --exact verdict tells whether a built
# h5ad and its restored package are byte-equal in X (the drivers use --exact only between two RESTORED packages).
# Waits for the C10b driver's pid to exit (capped by a timestamp), then the memory / disk gate. Last log line is one of
#   <ts> SMOKE DONE bulkid rc <n> ...      <ts> SMOKE NOT RUN <reason>
# usage: nohup bash experiments/658-prereg-1007-builds/smoke658.sh <c10b driver pid> [CAP YYYYmmddHHMM] &
set -uo pipefail
REPO=~/Code/virtual-cell-challenge-2026; cd "$REPO" || exit 1
X=experiments/658-prereg-1007-builds; W=$X/work/smoke; LOG=$W/smoke658.log; mkdir -p "$W"
WAITPID="${1:?c10b driver pid}"; CAP="${2:-202610041200}"
H5=submissions/497dr-october-path-lodopct-spt-abc.h5ad; VCC=submissions/497dr-october-path-lodopct-spt-abc.vcc
log(){ printf '%s %s\n' "$(date -u +%FT%TZ)" "$*" >> "$LOG"; }
TERMINAL=0; finish(){ TERMINAL=1; log "$1"; exit 0; }
trap 'rc=$?; [ "$TERMINAL" = 1 ] || log "SMOKE NOT RUN unexpected exit (rc $rc)"' EXIT
echo $$ > "$W/smoke658.pid"
log "armed (pid $$): waiting for pid $WAITPID to exit, then the gate (pressure >= 40 %, free+inactive >= 20 GiB, disk >= 30 GiB, no builder); cap $CAP"
while ps -p "$WAITPID" > /dev/null 2>&1; do
  [ "$(date -u +%Y%m%d%H%M)" -ge "$CAP" ] && finish "SMOKE NOT RUN pid $WAITPID still alive at the cap $CAP"
  sleep 60; done
log "pid $WAITPID exited"
gate(){
  local pct fi disk
  pct="$(memory_pressure | awk '/free percentage/{gsub(/%/,"",$5); print $5}')"
  fi="$(vm_stat | awk '/Pages free/{f=$3}/Pages inactive/{i=$3}END{gsub(/\./,"",f);gsub(/\./,"",i);printf "%d",(f+i)*16384/2^30}')"
  disk="$(df -k "$REPO" | awk 'NR==2{printf "%d", $4/1048576}')"
  log "gate: pressure free ${pct}% free+inactive ${fi} GiB disk ${disk} GiB"
  [ "${pct:-0}" -ge 40 ] && [ "${fi:-0}" -ge 20 ] && [ "${disk:-0}" -ge 30 ] \
    && ! uv run python -c 'import sys; sys.path.insert(0, "scripts"); import def_pipeline; raise SystemExit(0 if def_pipeline.other_builders() else 1)'
}
until gate; do
  [ "$(date -u +%Y%m%d%H%M)" -ge "$CAP" ] && finish "SMOKE NOT RUN the gate never opened before $CAP"
  sleep 120; done
[ -s "$H5" ] && [ -s "$VCC" ] || finish "SMOKE NOT RUN 497dr's h5ad or package is missing"
RST="$W/r497dr.h5ad"; rm -f "$RST"
/usr/bin/time -l uv run python scripts/restore_prediction_h5ad.py "$VCC" "$RST" > "$W/restore.log" 2>&1 || { rm -f "$RST"; finish "SMOKE NOT RUN restore failed (see $W/restore.log)"; }
log "restore: $(/usr/bin/grep -E 'real|maximum resident' "$W/restore.log" | tr -s ' ' | tr '\n' ' ')"
/usr/bin/time -l uv run python $X/bulkid658.py "$H5" "$RST" --exact --out "$W/bulkid_smoke.json" > "$W/bulkid.log" 2>&1; BID=$?
rm -f "$RST"
finish "SMOKE DONE bulkid rc $BID :: $(/usr/bin/grep -E 'bulk:|exact content|BULK IDENTICAL|DIFFERENT|real|maximum resident' "$W/bulkid.log" | tr -s ' ' | tr '\n' ' ' | cut -c1-600)"
