#!/usr/bin/env python3
"""Extract wait deadlines from the armed chains and committer specs and flag inversions.

A successor that gives up at T while its predecessor is still allowed to be running at T
will exit on the deadline and never build. A committer that gives up at the chain's
*release* deadline dies while that chain is still doing the work the commit is for.
"""
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]

CHAINS = {
    "chain479": ("experiments/479-pertarget-sharpen/chain479.sh", "start-gate"),
    "commit479": ("experiments/479-pertarget-sharpen/commit479.spec", "committer"),
    "chain495": ("experiments/495-lodo-cap6/chain495.sh", "start-gate"),
    "commit495": ("experiments/495-lodo-cap6/commit495.spec", "committer"),
    "chain479b": ("experiments/479-pertarget-sharpen/chain479b.sh", "release"),
    "chain476": ("experiments/476-def-day-merge/chain476.sh", "release"),
    "commit476": ("experiments/476-def-day-merge/commit476.spec", "committer"),
    "chain497": ("experiments/497-october-switches/chain497.sh", "release"),
    "commit497": ("experiments/497-october-switches/commit497.spec", "committer"),
    "chain498": ("experiments/498-combined-bulk/chain498.sh", "release"),
    "chain496p": ("experiments/496-uncalled-sign-hook/chain496p.sh", "release"),
    "commit496": ("experiments/496-uncalled-sign-hook/commit496.spec", "committer"),
    "chain497e": ("experiments/497-october-switches/chain497e.sh", "release"),
    "relaunch479b": ("experiments/479-pertarget-sharpen/relaunch479b.sh", "release"),
}

def deadlines(path: Path):
    text = path.read_text()
    if path.suffix == ".spec":
        m = re.search(r"^DEADLINE=(\d{12})", text, re.M)
        return [m.group(1)] if m else []
    return re.findall(r"-ge (\d{12})", text)

print(f"{'name':<16} {'kind':<12} deadlines")
found = {}
for name, (rel, kind) in CHAINS.items():
    ds = deadlines(ROOT / rel)
    found[name] = ds
    print(f"{name:<16} {kind:<12} {', '.join(ds) or '?'}")

# Code order after the 06:40Z rewire (chain479b.sh header, chain495.sh header, relaunch479b.sh):
#   495 runs first; 479b is started once 495 logs 'patched:' and then waits for 495 to EXIT;
#   476 waits for 495's pid AND for no chain479b.sh (its awk); 496p waits for 498 to exit.
print()
print("inversions (successor deadline <= the moment the predecessor is still allowed to run):")
pairs = [
    ("commit479", "chain479", "committer dies at chain479's release cap, while the three KAs run after that cap"),
    ("chain476", "chain479b", "476's awk waits until chain479b.sh is gone; 479b may keep waiting until its own later cap, then still has to build"),
    ("commit495", "chain479b", "479b will not build until the 495 builder is committed; the committer can die while 495 is still building"),
    ("chain497", "commit476", "497 requires 476's files at HEAD; both give up at the same minute, and 476's rehearsals run after its release"),
    ("chain496p", "chain498", "496p waits for chain498 to exit, and both give up at the same minute — 498 is allowed to start at that minute"),
]
for succ, pred, why in pairs:
    ts, tp = found[succ][0], found[pred][-1]
    flag = "INVERSION" if ts <= tp else "ok"
    print(f"  {flag}: {succ} {ts}  vs  {pred} {tp}")
    print(f"           {why}")
