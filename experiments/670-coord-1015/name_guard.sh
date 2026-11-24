#!/usr/bin/env bash
# 670 (coordinator [223d39], 2026-10-19): content-filter guard for one slot.
# Why: on 10-19 00:01Z the server refused 1019s1's model name ("…flagged by our automated content filter"). arm1015's 536b
# loop then retries the SAME name until its cap, so the slot would never submit. For R0 that would leave the final entry
# looping for days. This guard watches submit-<SL>.log. When a refusal appears after the guard started, and it is safe (the
# refusal came before any entry: no pending upload for the package in ~/.config/vcc/state.json), it stops that slot's
# arm1015, claims the next queue name (scripts/next_model_name.py, AGENTS §14), and re-launches arm1015 with the slot's
# saved arguments and the new name. At most 3 renames; every step is logged.
#   spec file (one value per line, in order): SL TAG READY-LOG RUN-LOG PKG NOTE PRIOR-LOG FLOOR CAP   (the name is NOT in it)
#   nohup bash scripts/detach.sh bash experiments/670-coord-1015/name_guard.sh <spec> <GUARD-CAP YYYYmmddHHMM> > …/name_guard-<SL>.nohup 2>&1 &
set -uo pipefail
REPO=~/Code/virtual-cell-challenge-2026; cd "$REPO" || exit 1
SPEC="${1:?spec}"; GCAP="${2:?guard cap}"
W="${NG_W:-experiments/670-coord-1015/work}"; ARM="${NG_ARM:-experiments/670-coord-1015/arm1015.sh}"
NAMER="${NG_NAMER:-uv run python scripts/next_model_name.py}"; STATE="${NG_STATE:-$HOME/.config/vcc/state.json}"
SL="$(sed -n 1p "$SPEC")"; TAG="$(sed -n 2p "$SPEC")"; RDY="$(sed -n 3p "$SPEC")"; RL="$(sed -n 4p "$SPEC")"; PKG="$(sed -n 5p "$SPEC")"
NOTE="$(sed -n 6p "$SPEC")"; PRIOR="$(sed -n 7p "$SPEC")"; FLOOR="$(sed -n 8p "$SPEC")"; CAP="$(sed -n 9p "$SPEC")"
[ -n "$CAP" ] || { echo "spec $SPEC incomplete"; exit 2; }
S="$W/submit-$SL.log"; LOG="$W/name_guard-$SL.log"; T0="$(date -u +%FT%TZ)"; N=0
log(){ printf '%s %s\n' "$(date -u +%FT%TZ)" "$*" >> "$LOG"; }
armpid(){ ps -axo pid,comm,args | awk -v s="$SL" '$2 ~ /^bash$/ { for (i = 3; i < NF; i++) if ($i ~ /arm1015\.sh$/ && $(i+1) == s) { print $1; exit } }'; }
refused(){ awk -v t0="$1" 'substr($0, 1, 20) > t0 && /content filter/ {f = 1} END {exit f ? 0 : 1}' "$S" 2>/dev/null; }
pending(){ python3 - "$STATE" "$PKG" <<'PY'
import json, os, sys
try:
    d = json.load(open(sys.argv[1]))
except Exception:
    sys.exit(1)
for r in (d.get("profiles", {}).get("default", {}).get("pending_uploads", {}) or {}).values():
    if os.path.realpath(str(r.get("local_path", ""))) == os.path.realpath(sys.argv[2]):
        sys.exit(0)
sys.exit(1)
PY
}
echo $$ > "$W/name_guard-$SL.pid"; log "armed (pid $$, T0 $T0): $SL $TAG, renames on a content-filter refusal (max 3), guard cap $GCAP"
SINCE="$T0"
while [ "$(date -u +%Y%m%d%H%M)" -lt "$GCAP" ]; do
  if /usr/bin/grep -qE "SUBMITTED ok|CAPTURED own entry" "$S" 2>/dev/null && [ -z "$(armpid)" ]; then log "$SL submitted: nothing to guard"; log "END"; exit 0; fi
  if refused "$SINCE"; then
    P="$(armpid)"
    if pending; then log "refusal seen but an upload is pending for $PKG: not touching it"; sleep 60; continue; fi
    [ "$N" -ge 3 ] && { log "3 renames already: leaving it for a human"; log "END"; exit 1; }
    [ -n "$P" ] && { kill "$P" 2>/dev/null; sleep 2; log "stopped arm1015 $SL (pid $P) after a content-filter refusal"; }
    N=$((N + 1)); NEW="$($NAMER --claim "$TAG-rename$N" 2>/dev/null | tail -1)"
    [ -n "$NEW" ] || { log "could not claim a new name"; log "END"; exit 1; }
    printf '%s renamed by name_guard (rename %s): new name claimed from the queue\n' "$(date -u +%FT%TZ)" "$N" >> "$S"
    nohup bash scripts/detach.sh bash "$ARM" "$SL" "$TAG" "$RDY" "$RL" "$PKG" "$NEW" "$NOTE (renamed $N: content filter)" "$PRIOR" "$FLOOR" "$CAP" > "$W/arm-$SL-rename$N.nohup" 2>&1 &
    log "re-launched arm1015 $SL with name '$NEW' (rename $N)"; SINCE="$(date -u +%FT%TZ)"; sleep 30
  fi
  sleep 60
done
log "guard cap $GCAP reached"; log "END"
