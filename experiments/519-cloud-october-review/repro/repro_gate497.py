#!/usr/bin/env python3
"""finish497dr.sh / gate497.py versus the checks chain497.sh actually ran, plus the stale-log ready().

Reads source. Does not invoke gate497.py (that reads 497's artifacts and --write rewrites its log).
"""
from __future__ import annotations

import re
import subprocess
import sys
import tempfile
from pathlib import Path

REPO = Path(__file__).resolve().parents[3]
bad = 0


def expect(label, cond, detail=""):
    global bad
    print(f"{'OK' if cond else 'UNEXPECTED'}  {label}" + (f": {detail}" if detail != "" else ""))
    if not cond:
        bad += 1


def main() -> int:
    chain = (REPO / "experiments/497-october-switches/chain497.sh").read_text()
    finish = (REPO / "experiments/497-october-switches/finish497dr.sh").read_text()
    gate = (REPO / "experiments/497-october-switches/gate497.py").read_text()
    expect("chain497 requires DONE. Next before relabel",
           'grep -q "DONE. Next"' in chain, "")
    expect("finish497dr does not look for DONE. Next", "DONE. Next" not in finish, "")
    expect("gate497 does not look for DONE. Next or postbuild.json",
           "DONE. Next" not in gate and "postbuild" not in gate, "")
    expect("gate497 docstring names gate460.py", "gate460.py" in gate, "")
    gate_line = next(l for l in finish.splitlines() if "gate497.py" in l)
    expect("finish logs the gate rc and does not exit on it",
           gate_line.strip().endswith('log "gate497 rc=$?"') and "||" not in gate_line, gate_line.strip())
    expect("chain497 also ignores the gate rc",
           'log "gate497 rc=$?"' in chain and "gate497.py --write" not in chain, "")

    day = (REPO / "experiments/497-october-switches/work/run/def_day.log").read_text()
    expect("497's def_day.log has no DONE. Next", "DONE. Next" not in day, "")
    expect("the log's postbuild PASS line is the manual continuation",
           "postbuild PASS (manual re-run" in day, "")
    # the zero-gene regex, taken from the source, matches a hand-typed line
    m = re.search(r'ZERO block total \\\(target excluded\\\): 0 in 0 blocks', gate)
    expect("gate source contains the zero-gene regex", m is not None, "")
    forged = "2099-01-01T00:00:00Z someone typed: ZERO block total (target excluded): 0 in 0 blocks; lfc < -5: 99"
    expect("a forged log line satisfies that regex",
           re.search(r"ZERO block total \(target excluded\): 0 in 0 blocks", forged) is not None, "")

    # ready() in chain1003/1004/1005: the gate's rc is not the function's rc.
    ready_src = {}
    for name in ("chain1003.sh", "chain1004.sh", "chain1005.sh"):
        path = next((REPO / "experiments").glob(f"*/{name}"))
        text = path.read_text()
        ready_src[name] = "gate497.py --write" in text and "|| return" not in text.split("gate497.py --write", 1)[1][:200]
    expect("1003/1004/1005 call gate497 --write and do not return on its rc",
           all(ready_src.values()), ready_src)

    with tempfile.TemporaryDirectory() as td:
        log = Path(td) / "gates497.log"
        log.write_text("[497dr] CLIFF OK\n[497dr] PASS: config diff vs 460dr\nnovelty rc=0\n")
        script = r"""
            set -uo pipefail
            LOG=%s
            ready() {
              false
              grep -q '^\[497dr\] PASS' "$LOG" && grep -q '^\[497dr\] CLIFF OK' "$LOG" && grep -q 'novelty rc=0' "$LOG"
            }
            if ready; then echo STALE_PASS; else echo REJECTED; fi
        """ % log
        r = subprocess.run(["bash", "-c", script], text=True, capture_output=True, check=True)
        expect("a failing gate then a stale PASS log still readies", r.stdout.strip() == "STALE_PASS", r.stdout.strip())
        script2 = r"""
            set -uo pipefail
            LOG=%s
            ready() {
              false || return 1
              grep -q '^\[497dr\] PASS' "$LOG"
            }
            if ready; then echo STALE_PASS; else echo REJECTED; fi
        """ % log
        r = subprocess.run(["bash", "-c", script2], text=True, capture_output=True, check=True)
        expect("returning on the gate rc rejects the stale log", r.stdout.strip() == "REJECTED", r.stdout.strip())

    print("FAILS" if bad else "ALL DEMONSTRATED", bad)
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
