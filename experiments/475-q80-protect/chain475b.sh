#!/usr/bin/env bash
# 475b: the `fix` seeds only (the KA `orig` == 435 seed 0 already PASSED at 02:40Z, work/ka475.txt; run435 refuses to
# overwrite its output). Logs into chain475.log and ends with "<ts> END", which chain471c waits for.
REPO=~/Code/virtual-cell-challenge-2026; cd "$REPO" || exit 1
X=experiments/475-q80-protect; LOGD=$X/work/logs; LOG=$X/work/chain475.log
log(){ printf '%s %s\n' "$(date -u +%FT%TZ)" "$*" >> "$LOG"; }
/usr/bin/grep -q PASS $X/work/ka475.txt || { log "no KA PASS on record -- stopping"; log "END"; exit 1; }
ARGS="--bulk ps --mde-c 0.1 --mde-cap 0.1 --mde-min-count 0.3 --mde-coarse"
for s in 0 1 2; do
  log "START fix seed $s"; /usr/bin/time -l uv run python $X/run475.py fix base "$s" $ARGS > "$LOGD/fix-s$s.log" 2>&1; rc=$?
  log "fix seed $s rc=$rc $(grep '^\[475\]' "$LOGD/fix-s$s.log" | tail -1 | cut -c1-120)"
  [ $rc -eq 0 ] || { log "fix seed $s failed: $(grep -E 'Error|error' "$LOGD/fix-s$s.log" | tail -1 | cut -c1-200)"; log "END"; exit 1; }
done
uv run python $X/summary475.py >> "$LOG" 2>&1
log "END"
