#!/usr/bin/env python3
"""Slot clocks for the 12 h-late 09-28 fire, and chain1001b's whole-file wait.

A tiny interpreter of the until / in-loop cap / post-loop check. It does not run arm*.sh
and it does not call vcc. slot_entry_id is the real script, on the real 09-28 logs.
"""
from __future__ import annotations

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


def arm0928_slot1(now: int) -> str:
    """arm0928.sh:20-23. The 12:00 cap is inside the until, which is skipped once now >= floor."""
    floor, cap = 202609280001, 202609281200
    if now >= floor:
        return "submit"  # until body never runs
    # body would have stopped at the cap; not reached for any now < floor that we care about
    if now >= cap:
        return "stop-cap"
    return "wait"


def arm0928_slot2(now: int, captured: bool) -> str:
    """arm0928.sh:25-29. Cap is inside the until. No post-loop date check."""
    cap = 202609282330
    if captured:
        return "submit"
    if now >= cap:
        return "stop-cap"
    return "wait"


def arm0929_slot1(now: int) -> str:
    """arm0929.sh:20-23 and the post-loop check at line 33."""
    floor, cap = 202609290001, 202609291200
    if now >= cap:
        return "stop-post-loop"
    if now >= floor:
        return "submit"
    return "wait"


def chain0929_releases(now: int, log_has_armed: bool, arm0928_alive: bool) -> str:
    """chain0929.sh:12-16. 18:00 break arms anyway."""
    if now >= 202609280005 and log_has_armed and not arm0928_alive:
        return "arm-0929"
    if now >= 202609281800:
        return "arm-0929-anyway"
    return "wait"


def chain_abs(now: int, floor: int) -> str:
    return "release" if now >= floor else "wait"


def main() -> int:
    # Observed fire: re-arm 12:11Z, slot 1 submitting 12:11:29Z, slot 2 12:38Z. Both are 09-28.
    expect("slot 1 at 12:11Z submits (cap is dead)", arm0928_slot1(202609281211) == "submit", "")
    expect("slot 1 at 12:11Z is still 09-28", 202609281211 < 202609290000, "")
    expect("slot 2 at 12:38Z with slot 1 captured submits", arm0928_slot2(202609281238, True) == "submit", "")
    expect("slot 2 still waiting at 18:00 does not submit", arm0928_slot2(202609281800, False) == "wait", "")
    expect("slot 2 at 23:30 with no capture stops", arm0928_slot2(202609282330, False) == "stop-cap", "")
    # no post-loop check: a capture that arrives at 23:29 still submits, possibly past 23:30
    expect("slot 2 submits whenever CAPTURED is already present, even after 23:30",
           arm0928_slot2(202609282331, True) == "submit", "")

    # chain0929 on the afternoon of 09-28, after the waiters have exited.
    expect("chain0929 releases the same afternoon",
           chain0929_releases(202609281310, True, False) == "arm-0929", "")
    expect("the arm0929 it starts does not submit until 09-29 00:01",
           arm0929_slot1(202609281310) == "wait" and arm0929_slot1(202609290001) == "submit", "")
    expect("arm0929 started after 12:00Z 09-29 refuses", arm0929_slot1(202609291300) == "stop-post-loop", "")
    expect("chain0929 does not release while arm0928 is alive, before 18:00",
           chain0929_releases(202609281200, True, True) == "wait", "")
    expect("18:00 break arms 09-29 even if arm0928 is alive",
           chain0929_releases(202609281800, True, True) == "arm-0929-anyway", "")

    for label, floor in (
        ("chain1002", 202610012000), ("chain1003", 202610022000),
        ("chain1004", 202610032000), ("chain1005", 202610042000),
    ):
        expect(f"{label} does not release on 09-28", chain_abs(202609281310, floor) == "wait", "")

    # Real entry ids. A re-arm must not keep the killed waiter's capture.
    r = subprocess.run(
        [sys.executable, str(REPO / "scripts/slot_entry_id.py"), "453dr",
         str(REPO / "experiments/456-zs-riders/work/submit0928s1.log"),
         str(REPO / "experiments/456-zs-riders/work/submit0928s2.log")],
        text=True, capture_output=True, check=True)
    expect("453dr id is the post-re-arm capture", r.stdout.strip() == "q3VuTImh4FRGKXH2waCp", r.stdout.strip())
    r = subprocess.run(
        [sys.executable, str(REPO / "scripts/slot_entry_id.py"), "446bj2z",
         str(REPO / "experiments/456-zs-riders/work/submit0928s1.log"),
         str(REPO / "experiments/456-zs-riders/work/submit0928s2.log")],
        text=True, capture_output=True, check=True)
    expect("446bj2z id is the post-re-arm capture", r.stdout.strip() == "be90smZUDcIV7j5MgQ03", r.stdout.strip())

    # Synthetic: an earlier ARMED+CAPTURED is dropped when the tag is re-armed with no later capture.
    with tempfile.TemporaryDirectory() as td:
        log = Path(td) / "s.log"
        log.write_text(
            "2026-09-28T00:02:00Z ARMED slot 1: 453dr pkg=old\n"
            "2026-09-28T00:20:00Z CAPTURED own entry AAAAAAAAAAAAAAAAAAAA\n"
            "2026-09-28T12:11:05Z ARMED slot 1: 453dr pkg=new\n"
        )
        r = subprocess.run([sys.executable, str(REPO / "scripts/slot_entry_id.py"), "453dr", str(log)],
                           text=True, capture_output=True, check=True)
        expect("re-arm without a new capture yields none", r.stdout.strip() == "none", r.stdout.strip())

    # chain1001b.sh:13. The pattern matches a stale line; there is no process check and no date floor.
    pattern = r"CAPTURED own entry|NOT SUBMITTING|exiting| END$"
    with tempfile.TemporaryDirectory() as td:
        stale = Path(td) / "submit0930s2.log"
        stale.write_text("2026-09-26T17:00:00Z quoted from somewhere else ... END\n")
        r = subprocess.run(["grep", "-qE", pattern, str(stale)])
        expect("a quoted line ending in END releases chain1001b's until", r.returncode == 0, r.returncode)
        stale.write_text("2026-09-30T00:00:00Z slot 2 is the PI's call; exiting\n")
        r = subprocess.run(["grep", "-qE", pattern, str(stale)])
        expect("the word exiting anywhere in the file also releases it", r.returncode == 0, r.returncode)
        missing = Path(td) / "nope.log"
        r = subprocess.run(["grep", "-qE", pattern, str(missing)], stderr=subprocess.DEVNULL)
        expect("a missing 09-30 log does not release it", r.returncode != 0, r.returncode)

    # BSD-grep trap, on this box's grep. The slot chains use grep -E, which anchors.
    line = "2026-09-28T12:37:36Z END\nnot the end of the line END more\n"
    e = subprocess.run(["grep", "-E", "-c", r" END$"], input=line, text=True, capture_output=True)
    expect("grep -E ' END$' matches one line", e.stdout.strip() == "1", e.stdout.strip())

    print("FAILS" if bad else "ALL DEMONSTRATED", bad)
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
