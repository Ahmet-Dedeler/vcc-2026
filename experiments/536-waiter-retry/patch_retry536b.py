#!/usr/bin/env python3
"""536b (2026-10-04, [4d0207]): the slot waiters' retry loop keeps trying until the day's cap, with exponential backoff.

Why: 536's loop allows 6 attempts 90 s apart, about 10 minutes of retrying, while its cap is ~12 h away. On 10-04 a hotspot
drop outlasted that. 498p2g's waiter gave up at 01:20Z after 6 attempts and needed a manual `--resume` at 05:21Z, and slot 2
hit `[Errno 8] nodename nor servname` at 06:00Z. Final-phase uploads (Oct 22 - Nov 5) run on whatever network Ahmet has.

What changes in a 536-patched waiter (nothing else):
  - backoff 90 s, 180, 360, 720, then 1,440 s (24 min) per attempt (`VCC_RETRY_BASE_S` overrides the 90 s base, for tests);
  - attempts are bounded by the day's cap, not by a count (hard bound 60);
  - each sleep is clamped to end >= 60 s before the cap, and the loop stops ("no time left before the cap") instead of
    sleeping across it;
  - unchanged: resume a pending upload of THIS package; a fresh retry only before any entry existed; stop if an entry
    exists without a pending upload (a fresh submit could double-submit).
Refuses a file a live bash process is executing (bash reads scripts incrementally; AGENTS 6.2), and a file without
536's block. Idempotent.
usage: uv run python experiments/536-waiter-retry/patch_retry536b.py <arm*.sh> [...]     (--check: report only)
"""
import re
import subprocess
import sys
from pathlib import Path

OLD = r'''OUT="$(vcc submit "$PKG" -m "$NAME" --wait 2>&1)"; RC=$?; N536=1
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
  N536=$((N536 + 1)); done'''

NEW = r'''OUT="$(vcc submit "$PKG" -m "$NAME" --wait 2>&1)"; RC=$?; N536=1
# 536b (2026-10-04): back off 90 s -> 24 min and keep trying until this day's cap; never sleep across the cap
CAPS536=$(date -j -u -f %Y%m%d%H%M __CAP__ +%s 2>/dev/null); [ -n "$CAPS536" ] || CAPS536=$(( $(date +%s) + 21600 ))
B536=${VCC_RETRY_BASE_S:-90}
while [ $RC -ne 0 ] && [ $N536 -lt 60 ] && [ "$(date -u +%Y%m%d%H%M)" -lt __CAP__ ]; do
  S536=$(( B536 << (N536 < 5 ? N536 - 1 : 4) )); L536=$(( CAPS536 - $(date +%s) - 60 ))
  if [ $L536 -le 0 ]; then log "submit attempt $N536 failed rc=$RC; no time left before the cap: stopping"; break; fi
  [ $S536 -gt $L536 ] && S536=$L536
  P536="$(_pend536 | head -1)"
  if [ -n "$P536" ]; then
    log "submit attempt $N536 failed rc=$RC; resuming the pending upload $P536 in $S536 s: $(echo "$OUT" | tail -c 200)"; sleep $S536
    OUT="$(vcc submit "$PKG" -m "$NAME" --wait --resume "$P536" 2>&1)"; RC=$?
  elif echo "$OUT" | grep -qi "entry\|scor"; then
    log "submit attempt $N536 failed rc=$RC after an entry existed and no upload is pending: not retrying (a fresh submit could double-submit): $(echo "$OUT" | tail -c 200)"; break
  else
    log "submit attempt $N536 failed rc=$RC before any entry; fresh retry in $S536 s: $(echo "$OUT" | tail -c 200)"; sleep $S536
    OUT="$(vcc submit "$PKG" -m "$NAME" --wait 2>&1)"; RC=$?
  fi
  N536=$((N536 + 1)); done'''


def live(path: Path) -> list[str]:
    """pids of bash processes executing this file: interpreter + "<dir>/<name>" in the arguments (waiters are launched as
    `bash experiments/<exp>/arm100X.sh ...`), never a bare pgrep -f (AGENTS 6.2)."""
    key = f"{path.parent.name}/{path.name}"
    out = subprocess.run(["ps", "-axo", "pid=,comm=,args="], capture_output=True, text=True).stdout
    rows = [ln.split(None, 2) for ln in out.splitlines()]
    return [r[0] for r in rows if len(r) == 3 and r[1] in ("bash", "/bin/bash") and key in r[2]]


def main() -> int:
    check = "--check" in sys.argv
    files = [a for a in sys.argv[1:] if a != "--check"]
    rc = 0
    for f in files:
        p = Path(f).resolve(); t = p.read_text()
        if "536b (2026-10-04)" in t:
            print(f"{f}: already 536b"); continue
        # 536's loop, with ANY attempt count (the coordinator raised 6 to 120 in arm1004..1007) and the file's own cap
        pat = re.escape(OLD).replace(re.escape("[ $N536 -lt 6 ]"), r"\[ \$N536 -lt (\d+) \]").replace(re.escape("__CAP__"), r"(\d{12})")
        ms = list(re.finditer(pat, t))
        if len(ms) != 1:
            print(f"{f}: 536's block not found exactly once -- refused"); rc = 1; continue
        m = ms[0]; cnt, cap = m.group(1), m.group(2)
        pids = live(p)
        if pids:
            print(f"{f}: a live bash process is executing it (pid {', '.join(pids)}) -- refused (AGENTS 6.2)"); rc = 1; continue
        if check:
            print(f"{f}: would patch (cap {cap}, was {cnt} attempts x 90 s)"); continue
        tmp = p.with_suffix(".sh.tmp536b")
        tmp.write_text(t[:m.start()] + NEW.replace("__CAP__", cap) + t[m.end():])
        tmp.replace(p)   # atomic rename
        print(f"{f}: patched to 536b (cap {cap}; was {cnt} attempts x 90 s)")
    return rc


if __name__ == "__main__":
    sys.exit(main())
