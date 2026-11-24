#!/usr/bin/env python3
"""552 (R7 owner [d80a10], 10-16): run r7a551.access() up to the truth boundary and stop there, opening nothing.

panel() is replaced by a sentinel that raises before its first line, so a PASS means every check access() makes before
it touches a panel passed: the claim's sha256, the ledger's single 'R7a armed' line, the claim's F, the instrument record
(with plan()'s pinned admission diff), the packages, the external scripts (with the ledger's driver pin) and arm_gate()'s
window. Nothing is written: say() prints to this process's stdout, plan() runs with write=False, no bytecode is cached.
The 10-16 06:34Z start failed at the instrument-record check, which a plan() dry run does not exercise; this test does.
usage: uv run python experiments/552-r7-reserve/test_access552.py
"""
from __future__ import annotations

import sys
from pathlib import Path

sys.dont_write_bytecode = True
ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "experiments/551-r6-c10b"))
import r7a551 as M  # noqa: E402


class Boundary(Exception):
    pass


def boundary(p, F):
    raise Boundary(f"{p}, F {F}")


def main() -> int:
    M.panel = boundary
    try:
        M.access()
    except Boundary as e:
        ok, why = M.G.window_ok()
        print(f"window now: {'OPEN' if ok else 'CLOSED'} ({why}); driver pin {M.driver_pin()}")
        print(f"PASS: access() reached the truth boundary (first panel {e}) with every pre-access check passing; nothing opened")
        return 0 if ok else 2
    print("FAIL: access() returned without reaching a panel")
    return 1


if __name__ == "__main__":
    sys.exit(main())
