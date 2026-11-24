"""Does r7a551.arm_gate actually replace every gate551 lookup, or can C10b's window and the
10-06 12:00Z limit still fire?

gate551 functions call window_ok / thresholds / GATE_LOG / HEAVY_LOG as bare names. This loads the real
gate551 (it does not import guard551 and does not write unless wait_gate/heavy run) and applies the same
assignments arm_gate makes, then calls the patched callables with synthetic timestamps.
"""
from __future__ import annotations

import datetime as dt
import dis
import importlib
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "experiments/551-r6-c10b"))
G = importlib.import_module("gate551")

UTC = dt.timezone.utc


def arm(start="2026-10-08T06:00Z", end="2026-10-08T18:00Z"):
    ws = dt.datetime.strptime(start, "%Y-%m-%dT%H:%MZ").replace(tzinfo=UTC)
    we = dt.datetime.strptime(end, "%Y-%m-%dT%H:%MZ").replace(tzinfo=UTC)
    claim = {"window": {"start": start, "end": end}}

    def window_ok(now=None):
        now = now or dt.datetime.now(UTC)
        if not (ws <= now < we):
            return False, f"R7a's window is {start} - {end}"
        return True, "R7a window open"

    G.window_ok = window_ok
    G.thresholds = lambda now=None: (22.0, 55.0, 400.0, "R7a access floors")
    G.GATE_LOG = Path("/tmp/r7a-gate-should-not-be-written.log")
    G.HEAVY_LOG = Path("/tmp/r7a-heavy-should-not-be-written.log")
    return claim


def globals_of(fn):
    """Names the function loads from its globals (call-time), vs names bound as locals or defaults."""
    load, store = set(), set()
    for ins in dis.get_instructions(fn):
        if ins.opname in ("LOAD_GLOBAL", "LOAD_GLOBAL_MODULE", "LOAD_GLOBAL_BUILTIN"):
            load.add(ins.argval)
        elif ins.opname in ("STORE_FAST", "STORE_DEREF"):
            store.add(ins.argval)
    return sorted(load), sorted(store)


def main():
    names = ("window_ok", "thresholds", "floor_gib", "gate", "wait_gate", "heavy")
    print("=== bare-name lookups (LOAD_GLOBAL) in gate551 callables ===")
    interesting = {"window_ok", "thresholds", "GATE_LOG", "HEAVY_LOG", "floor_gib", "OVERRIDE",
                   "NO_START_AFTER", "BLACKOUT", "PHASE", "ACCESS"}
    for n in names:
        load, store = globals_of(getattr(G, n))
        hit = [x for x in load if x in interesting]
        print(f"{n}: reads globals {hit}; stores locals {[s for s in store if s in interesting]}")

    # floor_gib calls thresholds() by global lookup; confirm the bytecode name
    print("floor_gib LOAD_GLOBAL names:", globals_of(G.floor_gib)[0])

    arm()
    # also set PHASE the way drive551 does, to show the patch ignores it
    G.PHASE["name"] = "access"

    cases = [
        ("inside R7a window", dt.datetime(2026, 10, 8, 12, 0, tzinfo=UTC)),
        ("C10b window 10-04 08:00Z (would be open for drive551)", dt.datetime(2026, 10, 4, 8, 0, tzinfo=UTC)),
        ("after 10-06 12:00Z ledger limit", dt.datetime(2026, 10, 7, 0, 0, tzinfo=UTC)),
        ("coordinator blackout 10-04 13:00Z", dt.datetime(2026, 10, 4, 13, 0, tzinfo=UTC)),
        ("before the recorded window", dt.datetime(2026, 10, 8, 5, 59, tzinfo=UTC)),
        ("at window end (exclusive)", dt.datetime(2026, 10, 8, 18, 0, tzinfo=UTC)),
    ]
    print("=== patched window_ok (PHASE is 'access', the C10b setting) ===")
    for label, now in cases:
        ok, msg = G.window_ok(now)
        print(f"  {label}: ok={ok} msg={msg}")

    print("=== patched thresholds, PHASE='access' ===")
    print(" ", G.thresholds(dt.datetime(2026, 10, 4, 8, 0, tzinfo=UTC)))

    print("=== floor_gib with no override follows the patched thresholds ===")
    # point OVERRIDE at a path that does not exist
    G.OVERRIDE = ROOT / "experiments/555-cloud-r7-review/repro/no-such-override.json"
    print(" ", G.floor_gib())

    print("=== floor_gib WITH an override file ignores the patched 22 GiB floor ===")
    ov = ROOT / "experiments/555-cloud-r7-review/repro/tmp_override.json"
    ov.write_text(json.dumps({"min_free_gib": 10, "by": "stale", "reason": "left over from C10b"}))
    G.OVERRIDE = ov
    try:
        print(" ", G.floor_gib())
    finally:
        ov.unlink()

    print("=== original window_ok still rejects 10-07 if the patch is not applied ===")
    # reimport fresh
    import importlib
    importlib.reload(G)
    # reload restores originals only if we didn't keep the patched functions... reload re-executes the module. Good.
    ok, msg = G.window_ok(dt.datetime(2026, 10, 7, 0, 0, tzinfo=UTC))
    print(f"  unpatched 2026-10-07: ok={ok} msg={msg}")
    ok, msg = G.window_ok(dt.datetime(2026, 10, 4, 8, 0, tzinfo=UTC))
    print(f"  unpatched 2026-10-04 08:00 PHASE=None: ok={ok} msg={msg}")
    G.PHASE["name"] = "access"
    ok, msg = G.window_ok(dt.datetime(2026, 10, 4, 8, 0, tzinfo=UTC))
    print(f"  unpatched 2026-10-04 08:00 PHASE=access: ok={ok} msg={msg}")
    ok, msg = G.window_ok(dt.datetime(2026, 10, 8, 12, 0, tzinfo=UTC))
    print(f"  unpatched 2026-10-08 12:00 PHASE=access: ok={ok} msg={msg}")


if __name__ == "__main__":
    main()
