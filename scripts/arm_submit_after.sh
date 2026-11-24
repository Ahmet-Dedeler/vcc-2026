#!/usr/bin/env bash
# Arm a preregistered submission to fire once a prior entry is TERMINAL (scored or failed) and a UTC
# floor has passed -- the "second slot of the day, never two in flight" pattern.
#   scripts/arm_submit_after.sh <prior-entry-id> <YYYYmmddHHMM-UTC floor> <package.vcc> <sha256> "<joke name>" <logfile>
# Same checks as arm_submit.sh: hash re-verified, novelty ledger consulted, account can submit, -d EMPTY.
#
# 2026-09-05 00:00Z change -- CAPTURE BEFORE YOU OVERWRITE. The board API resolves only the team's
# CURRENT entry; every earlier id returns not_found (253 work/vcc-status.log, 6/6). `vcc submit --wait`
# prints the scaled members only. So this script now saves the prior entry's COMPLETE
# `vcc --json status` record to data/board-entries/<id>.json (scripts/capture_board_entry.py, validated,
# atomic) BEFORE the floor wait, and refuses to submit if that capture is not on disk: a second slot
# must never destroy the raw evidence the first slot was spent for. It also captures its own entry
# after scoring. It is exec'd by arm_second_slot.sh with the same six positionals as before.
set -uo pipefail
PRIOR="${1:?prior entry id}"; FLOOR="${2:?UTC floor YYYYmmddHHMM}"; PKG="${3:?package}"; SHA="${4:?sha256}"; NAME="${5:?model name}"; LOG="${6:?log}"
REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"; cd "$REPO"
source .venv/bin/activate 2>/dev/null || true
log(){ printf '%s %s\n' "$(date -u +%FT%TZ)" "$*" >> "$LOG"; }
entry_id_of(){ local s a; s="$(cat)"; a="$(printf '%s\n' "$s" | grep -o 'entry [A-Za-z0-9]\{20\}' | head -1 | awk '{print $2}')"; [ -n "$a" ] && { echo "$a"; return; }; printf '%s\n' "$s" | grep -o '\b[A-Za-z0-9]\{20\}\b' | grep -v -E '^[0-9]+$' | head -1; }   # 2026-09-16: prefer the 'entry <id>' line (the scoring-job token is a 20-char lowercase twin); stdin read once
CAP="data/board-entries/$PRIOR.json"
log "ARMED-AFTER prior=$PRIOR floor=${FLOOR}Z pkg=$PKG name='$NAME' pid=$$"

# 1. Capture the prior entry's complete record FIRST. It is already terminal when we get here
#    (arm_submit.sh's --wait returned), so this is one API call; the retries cover API blips.
#    Fail CLOSED: no capture on disk, no second submission.
if uv run python scripts/capture_board_entry.py "$PRIOR" --wait --attempts 40 --interval 30 >> "$LOG" 2>&1; then
  log "CAPTURED prior $PRIOR full status record -> $CAP"
else
  log "ERROR could not capture prior $PRIOR's full status record -- NOT submitting (the second slot must not erase the first slot's raw members); re-run scripts/capture_board_entry.py $PRIOR, then re-arm"
  exit 1
fi

# 2. Floor.
while [ "$(date -u +%Y%m%d%H%M)" -lt "$FLOOR" ]; do sleep 30; done

# 3. Prior terminal -- read from the validated capture (no API dependency); the capture only
#    validates terminal states, so this is a belt-and-braces re-read of the file we will cite.
ST="$(uv run python -c 'import json,sys; print(json.load(open(sys.argv[1]))["entry"]["status"])' "$CAP" 2>/dev/null)"
case "$ST" in
  published|failed|hidden_admin|superseded) log "prior $PRIOR terminal: $ST" ;;
  *) log "ERROR prior $PRIOR capture reads status '$ST' (not terminal) -- NOT submitting"; exit 1 ;;
esac

# 4. The same pre-checks as arm_submit.sh.
GOT="$(shasum -a 256 "$PKG" | cut -d' ' -f1)"
[ "$GOT" = "$SHA" ] || { log "ERROR sha mismatch got=$GOT want=$SHA -- NOT submitting"; exit 1; }
uv run python scripts/check_submission_novelty.py check "$PKG" >> "$LOG" 2>&1 || { log "ERROR novelty check refused -- NOT submitting"; exit 1; }
WHO="$(vcc --json whoami 2>&1)"; echo "$WHO" | grep -q '"can_submit": true' || { log "ERROR account cannot submit -- NOT submitting"; exit 1; }

# 5. Submit, -d EMPTY (AGENTS.md section 14), --wait so the scaled members land in this log.
log "SUBMITTING $PKG as '$NAME' (description empty)"
OUT="$(vcc submit "$PKG" -m "$NAME" --wait 2>&1)"; RC=$?
log "SUBMIT rc=$RC"; printf '%s\n' "$OUT" >> "$LOG"
[ $RC -eq 0 ] && log "SUBMITTED ok" || log "ERROR submit failed rc=$RC"
ID="$(printf '%s\n' "$OUT" | entry_id_of)"
[ $RC -eq 0 ] && [ -n "$ID" ] && uv run python scripts/check_submission_novelty.py record "$PKG" --entry-id "$ID" --note "$NAME" >> "$LOG" 2>&1

# 6. Capture our own complete record too (it stays resolvable while current, but the raw members
#    belong on disk the moment they exist). Best effort: the submission is already made.
if [ -n "$ID" ]; then
  if uv run python scripts/capture_board_entry.py "$ID" --wait --attempts 40 --interval 30 --package "$PKG" --sha256 "$SHA" --note "$NAME" >> "$LOG" 2>&1; then
    log "CAPTURED own entry $ID -> data/board-entries/$ID.json"
  else
    log "WARN could not capture own entry $ID; re-run: scripts/capture_board_entry.py $ID --wait"
  fi
else
  log "WARN no entry id parsed from the submit output; nothing captured"
fi
