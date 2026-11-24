#!/usr/bin/env bash
# Arm a preregistered submission to fire at a UTC slot reset.
#   scripts/arm_submit.sh <YYYYmmddHHMM-UTC> <package.vcc> <expected-sha256> "<joke model name>" <logfile>
# Waits until the UTC minute, re-verifies the package hash, checks the account can submit and
# nothing is in flight, then runs `vcc submit` with -m only (AGENTS.md §14: -d stays EMPTY) and
# --wait so the score lands in the log. Every step is logged; any failed check aborts without
# submitting. Kill with: kill "$(cat <logfile>.pid)".
set -uo pipefail
WHEN="${1:?UTC minute YYYYmmddHHMM}"; PKG="${2:?package}"; SHA="${3:?sha256}"; NAME="${4:?model name}"; LOG="${5:?log}"
REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"; cd "$REPO"
source .venv/bin/activate 2>/dev/null || true
log(){ printf '%s %s\n' "$(date -u +%FT%TZ)" "$*" >> "$LOG"; }
log "ARMED pkg=$PKG when=${WHEN}Z name='$NAME' pid=$$"
while [ "$(date -u +%Y%m%d%H%M)" -lt "$WHEN" ]; do sleep 20; done
log "SLOT reached, running pre-checks"
GOT="$(shasum -a 256 "$PKG" | cut -d' ' -f1)"
[ "$GOT" = "$SHA" ] || { log "ERROR sha mismatch got=$GOT want=$SHA -- NOT submitting"; exit 1; }
WHO="$(vcc --json whoami 2>&1)"; log "whoami: $(echo "$WHO" | head -c 300)"
echo "$WHO" | grep -q '"can_submit": true' || { log "ERROR account cannot submit -- NOT submitting"; exit 1; }
uv run python scripts/check_submission_novelty.py check "$PKG" >> "$LOG" 2>&1 || { log "ERROR novelty check refused (rc=$?) -- NOT submitting"; exit 1; }
log "SUBMITTING $PKG as '$NAME' (description empty)"
OUT="$(vcc submit "$PKG" -m "$NAME" --wait 2>&1)"; RC=$?
if [ $RC -ne 0 ] && ! echo "$OUT" | grep -qi "entry\|scor"; then log "submit attempt 1 failed rc=$RC (retrying once in 90 s): $(echo "$OUT" | tail -c 200)"; sleep 90; OUT="$(vcc submit "$PKG" -m "$NAME" --wait 2>&1)"; RC=$?; fi
log "SUBMIT rc=$RC"; printf '%s\n' "$OUT" >> "$LOG"
[ $RC -eq 0 ] && log "SUBMITTED ok" || log "ERROR submit failed rc=$RC"
[ $RC -eq 0 ] && uv run python scripts/check_submission_novelty.py record "$PKG" --entry-id "$(echo "$OUT" | grep -o '[A-Za-z0-9]\{20\}' | head -1)" --note "$NAME" >> "$LOG" 2>&1
