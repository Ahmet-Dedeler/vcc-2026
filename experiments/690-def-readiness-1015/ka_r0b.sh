#!/usr/bin/env bash
# 690 ka_r0b (readiness [401d9b], 2026-10-22): the known-answer guard for planb_r0.sh's lowdose mode. r0b reruns r0's recipe
# with only DEF_DAY_ALLOW_LOW_DOSE=1 added (a postbuild-row override, not a build input), so its final output must be r0's lvl3
# byte for byte (output.sha256 <KA_SHA>). Waits for r0b's def_day "final:" line and its final meta (written before postbuild,
# well before r0ready692's prep and arm1015's submit), then:
#   - match: logs "KA PASS" and exits;
#   - mismatch, or def_day ends without a final meta: touches experiments/670-coord-1015/work/NOGO-R0 (arm1015 R0 stops on it
#     at every wait and before the submit) with the reason, and logs "KA FAIL".
# Cap: KA_CAP (default 8 h). Terminal line, always: "<ts> END".
#   nohup bash scripts/detach.sh bash experiments/690-def-readiness-1015/ka_r0b.sh <KA_SHA> > …/work/ka_r0b.nohup 2>&1 &
set -uo pipefail
REPO=~/Code/virtual-cell-challenge-2026; cd "$REPO" || exit 1
KA="${1:?expected output sha256}"
W=experiments/692-def-r0/work/r0b; DL="$W/def_day.log"; META="$W/r0def_final.h5ad.meta.json"
NOGO=experiments/670-coord-1015/work/NOGO-R0
LOG="${KA_LOG:-experiments/690-def-readiness-1015/work/ka_r0b.log}"
log(){ printf '%s %s\n' "$(date -u +%FT%TZ)" "$*" >> "$LOG"; }
trap 'log "END"' EXIT
end=$(( $(date +%s) + ${KA_CAP:-28800} ))
alive(){ ps -axo comm=,args= | awk '$1 ~ /bash$/ && /scripts\/def_day\.sh data\/challenge-def experiments\/692-def-r0\/work\/r0b /' | /usr/bin/grep -q .; }
log "ka_r0b armed (pid $$): want $W's final output ${KA:0:16}; NOGO on a mismatch: $NOGO"
until /usr/bin/grep -q " final: " "$DL" 2>/dev/null && [ -s "$META" ]; do
  if [ "$(date +%s)" -ge "$end" ]; then log "KA FAIL: no final meta by the cap"; echo "ka_r0b: no final meta in $W by the cap" > "$NOGO"; exit 1; fi
  if [ -s "$DL" ] && ! alive && ! /usr/bin/grep -q " final: " "$DL"; then
    sleep 30; alive || /usr/bin/grep -q " final: " "$DL" || { log "KA FAIL: def_day ended without a final build: $(tail -1 "$DL" | cut -c1-200)"
      echo "ka_r0b: r0b def_day ended without a final build" > "$NOGO"; exit 1; }
  fi
  sleep 30
done
sleep 5
GOT="$(python3 -c 'import json,sys; print(json.load(open(sys.argv[1]))["output"]["sha256"])' "$META" 2>/dev/null)"
if [ "$GOT" = "$KA" ]; then
  log "KA PASS: r0b's final output $GOT == r0's lvl3 (the rerun reproduced it byte for byte); $(/usr/bin/grep ' final: ' "$DL" | tail -1 | cut -c1-160)"
else
  log "KA FAIL: r0b's final output ${GOT:-<none>} != ${KA}; NOGO-R0 written"
  echo "ka_r0b: r0b final output ${GOT:-none} != r0 lvl3 $KA ($(date -u +%FT%TZ))" > "$NOGO"; exit 1
fi
