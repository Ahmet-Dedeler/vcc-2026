#!/usr/bin/env python3
"""Print the board entry id a slot waiter captured for <tag>, or "none".

Slot logs (`submit<day>s<n>.log`, written by the arm*.sh waiters) accumulate across kills and re-arms, so the last
"CAPTURED own entry" line in a file need not belong to <tag>. This scans each log in order and returns the CAPTURED id that
follows the LAST "ARMED slot <n>: <tag> " line, stopping at the next ARMED line (a re-arm) in the same file. When several
logs hold the tag, the last capture found wins (later arguments are later days).

usage: slot_entry_id.py <tag> <log> [<log> ...]
"""
import re
import sys
from pathlib import Path

tag, logs = sys.argv[1], sys.argv[2:]
armed = re.compile(r"ARMED slot [12]: (\S+) ")
cap = re.compile(r"CAPTURED own entry ([A-Za-z0-9]{20})\b")
found = "none"
for log in logs:
    p = Path(log)
    if not p.is_file():
        continue
    current, got = None, None
    for line in p.read_text(errors="replace").splitlines():
        m = armed.search(line)
        if m:
            current = m.group(1)
            if current == tag:
                got = None            # a (re-)arm of this tag: only a capture after it counts
            continue
        c = cap.search(line)
        if c and current == tag:
            got = c.group(1)
    if got:
        found = got
print(found)
