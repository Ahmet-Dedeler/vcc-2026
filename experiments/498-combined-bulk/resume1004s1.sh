#!/usr/bin/env bash
# Resume 10-04 slot 1 (498p2g): its upload reached 2.73 / 2.77 GB, then the network dropped (DNS, 00:59-01:20Z) and
# arm1004.sh's 6 resume attempts all failed. The entry JWSghc0o7MN5D6o8Iadc exists in 10-04's allowance; this finishes the
# upload, then records and captures it as arm1004.sh would.  bash experiments/498-combined-bulk/resume1004s1.sh (FORCE=1)
set -uo pipefail
REPO=~/Code/virtual-cell-challenge-2026; cd "$REPO" || exit 1
W=experiments/498-combined-bulk/work; LOG=$W/submit1004s1.log
PKG=submissions/498p2g-cap15-partial50-clean4q80-lodo6-ballast-sharpenp2g-stacklodo05hi6load-ctxlevelr-mde010coarse-zr5.vcc
RL=experiments/428-zero-gene-audit/work/build_498p2g/run498p2g.log
NAME="Even if I knew the world would end tomorrow, I'd still plant an apple tree."
EID=JWSghc0o7MN5D6o8Iadc
log(){ printf '%s %s\n' "$(date -u +%FT%TZ)" "$*" >> "$LOG"; }
if [ "${FORCE:-0}" != 1 ]; then
  pmset -g batt | grep -q "AC Power" || { echo "not on AC power -- refusing (FORCE=1 overrides)"; exit 1; }
  GW="$(route -n get default 2>/dev/null | awk '/gateway/{print $2}')"
  [ "$GW" = "172.20.10.1" ] && { echo "default gateway $GW is an iPhone hotspot -- refusing (FORCE=1 overrides)"; exit 1; }
fi
SHA="$(grep 'package sha256:' "$RL" | tail -1 | awk '{print $NF}')"; GOT="$(shasum -a 256 "$PKG" | cut -d' ' -f1)"
[ -n "$SHA" ] && [ "$GOT" = "$SHA" ] || { log "RESUME refused: sha mismatch got=$GOT want=$SHA"; exit 1; }
log "RESUME (coordinator): vcc submit --resume $EID"
OUT="$(vcc submit "$PKG" -m "$NAME" --wait --resume "$EID" 2>&1)"; RC=$?
log "SUBMIT rc=$RC"; printf '%s\n' "$OUT" >> "$LOG"
if [ $RC -eq 0 ]; then
  uv run python scripts/check_submission_novelty.py record "$PKG" --entry-id "$EID" --note "498p2g slot 1 on 10-04 (resumed): $NAME" >> "$LOG" 2>&1
  uv run python scripts/capture_board_entry.py "$EID" --package "$PKG" --sha256 "$SHA" --note "498p2g slot 1 on 10-04 (resumed)" \
      --wait --attempts 40 --interval 30 >> "$LOG" 2>&1 \
    && log "CAPTURED own entry $EID -> data/board-entries/$EID.json" || log "WARN capture of $EID failed; re-run scripts/capture_board_entry.py $EID"
else
  log "ERROR resume failed rc=$RC; the pending upload stays in ~/.config/vcc/state.json for another --resume"
fi
log "END"
