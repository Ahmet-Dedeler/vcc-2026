#!/usr/bin/env bash
# Resume 10-03 slot 2 (471c) after the coordinator paused it at 10:53Z (Mac on battery over an iPhone hotspot,
# ~0.28 MB/s; ~1.6 GB left would have outlasted the battery and burned cellular data). The entry yzNoMVih0pKqoDjdqrXm
# exists in 10-03's allowance; this only finishes its upload, then records and captures it as arm1003.sh would.
# Run on AC power and a real network:  bash experiments/495-lodo-cap6/resume1003s2.sh   (FORCE=1 skips both checks)
set -uo pipefail
REPO=~/Code/virtual-cell-challenge-2026; cd "$REPO" || exit 1
W=experiments/495-lodo-cap6/work; LOG=$W/submit1003s2.log
PKG=submissions/471c-cap15-partial50-clean4q80-pertscale3-ballast-sharpen175-stackh1fload-ctxlevel-mde010coarse-zr5.vcc
RL=experiments/428-zero-gene-audit/work/build_471c/run471c.log
NAME="If you don't believe it or don't get it, I don't have the time to try to convince you, sorry. (Satoshi)"
EID=yzNoMVih0pKqoDjdqrXm
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
  uv run python scripts/check_submission_novelty.py record "$PKG" --entry-id "$EID" --note "471c slot 2 on 10-03 (resumed): $NAME" >> "$LOG" 2>&1
  uv run python scripts/capture_board_entry.py "$EID" --package "$PKG" --sha256 "$SHA" --note "471c slot 2 on 10-03 (resumed)" \
      --wait --attempts 40 --interval 30 >> "$LOG" 2>&1 \
    && log "CAPTURED own entry $EID -> data/board-entries/$EID.json" || log "WARN capture of $EID failed; re-run scripts/capture_board_entry.py $EID"
else
  log "ERROR resume failed rc=$RC; the pending upload stays in ~/.config/vcc/state.json for another --resume"
fi
log "END"
