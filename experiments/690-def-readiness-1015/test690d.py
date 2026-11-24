#!/usr/bin/env python3
"""690d tests: validate() on synthetic records, imported from a patched scratch copy of capture_board_entry.py.
usage: uv run python experiments/690-def-readiness-1015/test690d.py
"""
from __future__ import annotations

import copy
import importlib.util
import json
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
RESULTS: list[tuple[bool, str]] = []


def rec(ok: bool, name: str) -> None:
    RESULTS.append((bool(ok), name))
    print(f"{'PASS' if ok else 'FAIL'} {name}", flush=True)


def main() -> int:
    td = Path(tempfile.mkdtemp(prefix="t690d_"))
    try:
        (td / "scripts").mkdir()
        shutil.copy2(ROOT / "scripts/capture_board_entry.py", td / "scripts/capture_board_entry.py")
        r = subprocess.run([sys.executable, str(HERE / "patch690d.py"), "--root", str(td), "--check"], capture_output=True, text=True)
        rec(r.returncode == 0 and "final_warnings" in r.stdout, "D1 --check on the live script")
        r = subprocess.run([sys.executable, str(HERE / "patch690d.py"), "--root", str(td)], capture_output=True, text=True)
        rec(r.returncode == 0, "D1 apply")
        r = subprocess.run([sys.executable, str(HERE / "patch690d.py"), "--root", str(td)], capture_output=True, text=True)
        rec(r.returncode == 1, "D1 second apply refused")
        spec = importlib.util.spec_from_file_location("cbe", td / "scripts/capture_board_entry.py")
        m = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(m)
        val = json.loads((ROOT / "data/board-entries/dVHLws6KdMxudq0Yn5ap.json").read_text())["entry"]
        eid = val["entry_id"]
        rec(m.validate(val, eid) == [] and m.final_warnings(val) == [], "D2 a real validation entry stays valid, no warnings")
        v2 = copy.deepcopy(val); v2["rank"] = None
        rec("rank missing" in m.validate(v2, eid), "D2 a validation entry without rank still fails")
        v3 = copy.deepcopy(val); v3["anchor_version"] = None
        rec(any("anchor_version" in p for p in m.validate(v3, eid)), "D2 a validation entry without its stamp still fails")
        f1 = copy.deepcopy(val); f1.update(is_final=True, partition="test", rank=None)
        rec(m.validate(f1, eid) == [] and m.final_warnings(f1) == ["rank missing (dark final board)"],
            "D3 a final entry with members and no rank is VALID, warned")
        f2 = copy.deepcopy(f1); f2["anchor_version"] = None
        rec(m.validate(f2, eid) == [] and "anchor_version missing" in m.final_warnings(f2), "D3 a final entry without anchor: valid, warned")
        f3 = copy.deepcopy(f1); f3["score_fid"] = None
        rec(any("score_fid" in p for p in m.validate(f3, eid)), "D4 a final entry with a PARTIAL member set still fails")
        f4 = copy.deepcopy(f1)
        for k in m.RAW + m.SCALED + ("score_avg",):
            f4[k] = None
        rec(m.members_hidden(f4) and m.validate(f4, eid) == [], "D5 the all-null final entry is members_hidden and valid (531 F8)")
        f5 = copy.deepcopy(f1); f5["status"] = "scoring"
        rec(any("not terminal" in p for p in m.validate(f5, eid)), "D6 a non-terminal final entry still fails")
        src = (td / "scripts/capture_board_entry.py").read_text()
        rec('"final_warnings": final_warnings(entry),' in src, "D7 the wrapper records final_warnings")
    finally:
        shutil.rmtree(td, ignore_errors=True)
    n = sum(ok for ok, _ in RESULTS)
    print(f"\n{n}/{len(RESULTS)} PASS")
    return 0 if n == len(RESULTS) else 1


if __name__ == "__main__":
    sys.exit(main())
