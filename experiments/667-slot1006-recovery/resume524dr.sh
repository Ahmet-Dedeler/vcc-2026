#!/usr/bin/env bash
# 667 (coordinator, 2026-10-06 ~22:15Z): resume 524dr's upload after the 18:53Z reboot. 18:32Z's arm_submit created entry
# KTzIDR0Occ2jYzklkARE in 10-06's allowance (spent 1: K16) and the reboot cut its upload; vcc keeps it as a pending upload.
# This resumes THE SAME entry (never a fresh submit), records novelty and captures. If the entry is still queued at 23:30Z,
# guard_inflight cancels it so CPH's 10-07 00:01Z submit is not blocked. Cap: 23:45Z for the resume attempts.
set -uo pipefail
REPO=~/Code/virtual-cell-challenge-2026; cd "$REPO" || exit 1
X=experiments/667-slot1006-recovery; LOG=$X/work/submit524dr_late.log
PKG=submissions/524dr-october-recipe-abc.vcc; SHA=cbeadd654bc09d53ef5f362c8026b82b89a3d66a05b62811f3366269c5021156
NAME="The most dangerous worldviews are the worldviews of those who have never viewed the world."; EID=KTzIDR0Occ2jYzklkARE
log(){ printf '%s %s\n' "$(date -u +%FT%TZ)" "$*" >> "$LOG"; }
echo $$ > $X/work/resume524dr.pid
[ "$(shasum -a 256 $PKG | cut -d' ' -f1)" = "$SHA" ] || { log "resume524dr: package sha mismatch -- stopping"; log "END"; exit 1; }
log "resume524dr (coordinator): resuming $EID after the 18:53Z reboot"
N=0; RC=1; OUT=""
while [ "$(date -u +%Y%m%d%H%M)" -lt 202610062345 ]; do
  N=$((N + 1)); OUT="$(vcc submit "$PKG" -m "$NAME" --wait --resume "$EID" 2>&1)"; RC=$?
  [ $RC -eq 0 ] && break
  log "resume attempt $N failed rc=$RC: $(echo "$OUT" | tail -c 200 | tr '\n' ' ')"
  echo "$OUT" | /usr/bin/grep -qiE "no pending|not pending|nothing to resume|cancel" && break
  sleep 60
done
log "SUBMIT rc=$RC (resume524dr, $N attempts)"; printf '%s\n' "$OUT" >> "$LOG"
if [ $RC -eq 0 ]; then
  uv run python scripts/check_submission_novelty.py record "$PKG" --entry-id "$EID" --note "524dr slot 2 on 10-06 (resumed after the reboot): $NAME" >> "$LOG" 2>&1
  uv run python scripts/capture_board_entry.py "$EID" --package "$PKG" --sha256 "$SHA" --note "524dr on 10-06 (resumed)" --wait --attempts 40 --interval 30 >> "$LOG" 2>&1 \
    && log "CAPTURED own entry $EID -> data/board-entries/$EID.json" || log "WARN capture of $EID failed; re-run scripts/capture_board_entry.py $EID"
fi
log "END"
