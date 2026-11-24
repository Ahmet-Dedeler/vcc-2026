#!/usr/bin/env python3
"""Restart window in def_pin_panel_artifacts.py.

The apply path replaces the builder, and only then writes the JSON record and the
.builder-before.py sibling. A kill in that gap leaves a pinned builder and nothing
the runbook's --restore --record can read. A later def_day in a fresh work dir records
the pinned bytes as builder_sha256_before, so its exit trap "restores" them and exits 0.

This does not run --apply. It replays the three writes in the source's order against
temp files. Exits 0 while the source order is still replace-then-record.
"""
from __future__ import annotations

import hashlib
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
SRC = (ROOT / "scripts/def_pin_panel_artifacts.py").read_text()
DEF_DAY = (ROOT / "scripts/def_day.sh").read_text()


def sha(b: bytes) -> str:
    return hashlib.sha256(b).hexdigest()


def order() -> tuple[int, int, int]:
    # The apply block only. Earlier write_text calls are the dry-run helpers.
    apply = SRC.index("if args.apply:")
    replace = SRC.index("tmp.replace(BUILDER)", apply)
    record = SRC.index("args.out.write_text", apply)
    copy = SRC.index("copy.write_bytes(before)", apply)
    return replace, record, copy


def replay(order_name: str, steps: list[str], kill_after: str) -> dict:
    """steps are 'replace', 'record', 'copy'. Kill means stop before the next step."""
    import tempfile
    td = Path(tempfile.mkdtemp())
    builder = td / "builder.py"
    original, pinned = b"HEAD-BYTES\n", b"PINNED-BYTES\n"
    builder.write_bytes(original)
    record = td / "pin-artifacts.json"
    copy = td / "pin-artifacts.builder-before.py"
    before = builder.read_bytes()
    done = []
    for step in steps:
        if step == kill_after:
            break
        if step == "replace":
            tmp = builder.with_suffix(".py.tmp")
            tmp.write_bytes(pinned)
            tmp.replace(builder)
        elif step == "copy":
            copy.write_bytes(before)
        elif step == "record":
            record.write_text(json.dumps({"builder_sha256_before": sha(before)}))
        done.append(step)
    restorable = record.is_file() and copy.is_file()
    restored = None
    if restorable:
        blob = copy.read_bytes()
        if sha(blob) == json.loads(record.read_text())["builder_sha256_before"]:
            builder.write_bytes(blob)
            restored = sha(builder.read_bytes()) == sha(original)
    return {"order": order_name, "completed": done, "builder_sha": sha(builder.read_bytes())[:12],
            "original": sha(original)[:12], "pinned": sha(pinned)[:12],
            "record": record.is_file(), "before_copy": copy.is_file(), "restored_to_head": restored}


def main() -> int:
    replace, record, copy = order()
    print(f"source offsets: replace={replace} record={record} copy={copy}")
    assert replace < record < copy, (replace, record, copy)
    assert "pin-lodo.json" not in DEF_DAY.split("restore()")[1].split("# 1. setup")[0]
    assert 'def_pin_panel_artifacts.py --restore --record "$W/pin-artifacts.json"' in DEF_DAY

    current = ["replace", "record", "copy"]
    gap = replay("current", current, kill_after="record")
    print("kill after replace, before the record:", gap)
    assert gap["builder_sha"] == gap["pinned"]
    assert gap["restored_to_head"] is None

    # Fresh work dir after that kill. The new pin's "before" is the already-pinned file.
    # Restoring it puts the pin back and reports success.
    import tempfile
    td = Path(tempfile.mkdtemp())
    pinned = b"PINNED-BYTES\n"
    fresh_before = td / "fresh.builder-before.py"
    fresh_before.write_bytes(pinned)
    fresh_rec = td / "fresh.json"
    fresh_rec.write_text(json.dumps({"builder_sha256_before": sha(pinned)}))
    builder = td / "builder.py"
    builder.write_bytes(pinned)
    bsha0 = sha(builder.read_bytes())
    # trap: restore from THIS work dir's record, then compare to BSHA0
    builder.write_bytes(fresh_before.read_bytes())
    now = sha(builder.read_bytes())
    print(f"fresh work dir: BSHA0 {bsha0[:12]} after trap {now[:12]} match={now == bsha0}")
    assert now == bsha0 and now == sha(pinned)

    patched = replay("record-and-copy first", ["copy", "record", "replace"], kill_after="replace")
    print("patched order, kill before replace:", patched)
    assert patched["restored_to_head"] is True
    print("PIN WINDOW CONFIRMED: runbook restore needs both files, and they are written after the replace")
    return 0


if __name__ == "__main__":
    sys.exit(main())
