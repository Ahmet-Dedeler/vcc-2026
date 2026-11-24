#!/usr/bin/env python3
"""536 (2026-09-30): the slot waiters' submit step retries properly. On 09-30 the 458d upload died when the Mac slept on
battery (connection reset at 02:39Z) and the waiter's single fresh retry sat in `sleep 90` for ~2 h of sleep. This replaces
the two-line submit + one-retry block with a loop of up to 6 attempts before the file's own 12:00Z late-submit cap:
  - an interrupted upload of THIS package (a pending entry in ~/.config/vcc/state.json whose local_path is the package) is
    resumed with `vcc submit --resume <entry>` (no new entry);
  - a failure before any entry existed is retried fresh;
  - a failure after an entry existed with no pending upload stops (a fresh submit could double-submit).
Idempotent (refuses a patched file); the cap is read from the file's own "clock is past 12:00Z" line.
usage: uv run python experiments/536-waiter-retry/patch_retry536.py <arm*.sh> [...]"""
import re
import sys
from pathlib import Path

BLOCK = r'''# 536 (2026-09-30): retry loop -- resume an interrupted upload of THIS package, retry fresh only before any entry existed,
# stop if an entry exists without a pending upload; up to 6 attempts, never past this day's cap (__CAP__).
_pend536(){ python3 - "$PKG" <<'PY'
import json, os, sys
try:
    d = json.load(open(os.path.expanduser("~/.config/vcc/state.json")))
except Exception:
    sys.exit(0)
for eid, r in d.get("profiles", {}).get("default", {}).get("pending_uploads", {}).items():
    if os.path.realpath(str(r.get("local_path", ""))) == os.path.realpath(sys.argv[1]):
        print(eid)
PY
}
OUT="$(vcc submit "$PKG" -m "$NAME" --wait 2>&1)"; RC=$?; N536=1
while [ $RC -ne 0 ] && [ $N536 -lt 6 ] && [ "$(date -u +%Y%m%d%H%M)" -lt __CAP__ ]; do
  P536="$(_pend536 | head -1)"
  if [ -n "$P536" ]; then
    log "submit attempt $N536 failed rc=$RC; resuming the pending upload $P536 in 90 s: $(echo "$OUT" | tail -c 200)"; sleep 90
    OUT="$(vcc submit "$PKG" -m "$NAME" --wait --resume "$P536" 2>&1)"; RC=$?
  elif echo "$OUT" | grep -qi "entry\|scor"; then
    log "submit attempt $N536 failed rc=$RC after an entry existed and no upload is pending: not retrying (a fresh submit could double-submit): $(echo "$OUT" | tail -c 200)"; break
  else
    log "submit attempt $N536 failed rc=$RC before any entry; fresh retry in 90 s: $(echo "$OUT" | tail -c 200)"; sleep 90
    OUT="$(vcc submit "$PKG" -m "$NAME" --wait 2>&1)"; RC=$?
  fi
  N536=$((N536 + 1)); done
'''

for f in sys.argv[1:]:
    p = Path(f); t = p.read_text()
    if "_pend536" in t:
        print(f"{f}: already patched"); continue
    m = re.search(r'-ge (\d{12}) \] && stop "clock is past 12:00Z', t)
    if not m:
        raise SystemExit(f"{f}: no late-submit cap line")
    cap = m.group(1)
    lines = t.split("\n")
    i = [k for k, l in enumerate(lines) if l.startswith('OUT="$(vcc submit "$PKG" -m "$NAME" --wait 2>&1)"; RC=$?')]
    if len(i) != 1 or not lines[i[0] + 1].startswith('if [ $RC -ne 0 ] && ! echo "$OUT" | grep -qi "entry\\|scor"; then log "submit attempt 1 failed'):
        raise SystemExit(f"{f}: submit block not found exactly once")
    k = i[0]
    new = lines[:k] + BLOCK.replace("__CAP__", cap).rstrip("\n").split("\n") + lines[k + 2:]
    tmp = p.with_suffix(".sh.tmp536"); tmp.write_text("\n".join(new)); tmp.replace(p)
    print(f"{f}: patched (cap {cap})")
