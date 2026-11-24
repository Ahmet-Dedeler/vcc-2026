#!/usr/bin/env python3
"""708 AMENDMENT 4 test for decide708.py's robustness rule (science [d80a10], 10-16).

Six synthetic score708 records (tags tdec708-<k>, written to work/ and removed afterwards) and a synthetic R-KA proof that
quotes the LIVE def_day.sh sha256 with "690i APPLIED", so only the rule under test can block:
  1 both versions pass everything            -> weights MOVE, rider ON, swap ON
  2 the swap passes on all, fails floor-kept   -> swap OFF (floor-fragile)
  3 the rider's veto holds on all, not floor   -> rider ON (no veto)
  4 the rider's veto holds in both             -> rider VETOED
  5 no arms_floor at all                       -> FLOOR line says blocked; weights BLOCKED, rider ON (no veto), swap OFF
  6 the weights confirm on all, not floor-kept -> weights BLOCKED (floor-fragile)
usage: uv run python experiments/708-def-panel-repricing/test_decide708.py      (prints PASS or raises)
"""
from __future__ import annotations

import copy
import hashlib
import json
import subprocess
import sys
import tempfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
FOLDS = ("k562_gwps", "orion_hct116", "orion_hek293t")


def arm(mean, n_pos=3, lt3=None):
    return dict(per_fold={f: mean for f in FOLDS}, mean=mean, n_pos=n_pos, lt3_mean=mean if lt3 is None else lt3)


def record(**over):
    a = {"r": arm(0.0006), "s": arm(0.0020), "rs": arm(0.0026), "s|r": arm(0.0020), "A:h2.5c1.5": arm(0.0010)}
    rec = {"arms": copy.deepcopy(a), "arms_floor": copy.deepcopy(a),
           "floor": {"threshold_cpm": 13.4, "kept": {f: 250 for f in FOLDS}}}
    for k, v in over.items():
        ver, name = k.split("__")
        if v is None:
            rec.pop(ver, None)
        else:
            rec[ver][name.replace("_or_", "|").replace("_colon_", ":")] = v
    return rec


def run(tag, rec, r705="CLEARED", r707="NULL", des="h2.5c1.5"):
    p = HERE / "work" / f"score708-{tag}.json"
    p.write_text(json.dumps(rec))
    try:
        out = subprocess.run([sys.executable, str(HERE / "decide708.py"), tag, str(DAYOF), des, r705, r707, str(RLOG)],
                             capture_output=True, text=True, cwd=ROOT, check=True, env=ENV).stdout
    finally:
        p.unlink()
    return {ln.split()[0]: ln for ln in out.splitlines() if ln.strip()}


ENV = None


def main():
    global DAYOF, RLOG, ENV
    import os
    tmp = Path(tempfile.mkdtemp())
    ENV = dict(os.environ, R708_REHEARSAL="1")      # cases 1-6 test AMENDMENT 4's rule inside a rehearsal-exempt run
    DAYOF = tmp / "dayof.txt"
    DAYOF.write_text("RECOMMENDATION: H1 x2.5, CD4 x1.5 (+0.00217 over H1 x2 / CD4 x1)\n")
    live = hashlib.sha256((ROOT / "scripts/def_day.sh").read_bytes()).hexdigest()
    RLOG = tmp / "run690.log"
    RLOG.write_text(f"2026-10-16T00:00:00Z R-KA PASS: test\n2026-10-16T00:00:01Z R-KA tree: def_day.sh {live} 690i APPLIED\n")
    L = run("tdec708-1", record())
    assert L["WEIGHTS"].endswith("MOVE"), L["WEIGHTS"]
    assert "-> ON" in L["RIDER"] and "-> ON" in L["SWAP707"], (L["RIDER"], L["SWAP707"])
    assert "FINAL    weights h2.5c1.5 | rider ON | swap707 ON" in L["FINAL"], L["FINAL"]
    L = run("tdec708-2", record(arms_floor__s_or_r=arm(0.0003)))
    assert "-> OFF (floor-fragile" in L["SWAP707"], L["SWAP707"]
    L = run("tdec708-3", record(arms__r=arm(-0.0008, 0), arms_floor__r=arm(-0.0002, 1)))
    assert "-> ON (veto floor-fragile" in L["RIDER"], L["RIDER"]
    assert "judged on s|r" in L["SWAP707"], L["SWAP707"]
    L = run("tdec708-4", record(arms__r=arm(-0.0008, 0), arms_floor__r=arm(-0.0007, 0)))
    assert "VETOED" in L["RIDER"] and "judged on s:" in L["SWAP707"], (L["RIDER"], L["SWAP707"])
    rec = record(arms__r=arm(-0.0008, 0))
    rec.pop("arms_floor")
    L = run("tdec708-5", rec)
    assert "BLOCKED (AMENDMENT 4)" in L["FLOOR"], L["FLOOR"]
    assert "BLOCKED" in L["WEIGHTS"] and "-> OFF" in L["SWAP707"], (L["WEIGHTS"], L["SWAP707"])
    assert "-> ON (veto floor-fragile or floor missing" in L["RIDER"], L["RIDER"]
    L = run("tdec708-6", record(**{"arms__A_colon_h2.5c1.5": arm(0.0003), "arms_floor__A_colon_h2.5c1.5": arm(-0.0001, 1)}))
    assert "BLOCKED: floor-fragile" in L["WEIGHTS"], L["WEIGHTS"]
    # AMENDMENT 5: on a real bundle (no exemption) every action needs the rehearsal on file
    empty = tmp / "noreh"; empty.mkdir()
    ENV = {k: v for k, v in os.environ.items() if k != "R708_REHEARSAL"}
    ENV["R708_REHEARSAL_DIR"] = str(empty)
    L = run("tdec708-7", record())
    assert "MISSING -> 708 OFF" in L["REHEARSAL"], L["REHEARSAL"]
    assert "BLOCKED" in L["WEIGHTS"] and "-> OFF" in L["SWAP707"] and "-> ON" in L["RIDER"], (L["WEIGHTS"], L["SWAP707"])
    assert "FINAL    weights h2c1 | rider ON | swap707 OFF" in L["FINAL"], L["FINAL"]
    reh = tmp / "reh"; reh.mkdir()
    (reh / "kaA708-r708abcdef.json").write_text(json.dumps({"verdict": "PASS"}))
    (reh / "rehearse708.log").write_text("2026-10-16T21:40:00Z B PASS\n2026-10-16T21:45:00Z REHEARSAL PASS\n2026-10-16T21:45:00Z END\n")
    (reh / "rehearse708b.log").write_text("2026-10-16T22:50:00Z day708 rc 0; its last line: x\n2026-10-16T22:50:01Z DAYPATH same choice\n")
    ENV["R708_REHEARSAL_DIR"] = str(reh)
    L = run("tdec708-8", record())
    assert L["REHEARSAL"].startswith("REHEARSAL ON FILE"), L["REHEARSAL"]
    assert "FINAL    weights h2.5c1.5 | rider ON | swap707 ON" in L["FINAL"], L["FINAL"]
    (reh / "rehearse708.log").write_text("2026-10-16T12:24:33Z REHEARSAL FAIL A\n2026-10-16T12:24:33Z END\n")
    L = run("tdec708-9", record())
    assert "MISSING" in L["REHEARSAL"], L["REHEARSAL"]
    print("test_decide708 PASS: 9 cases (both pass; floor-fragile swap; floor-fragile veto; veto in both; floor missing; "
          "floor-fragile weights; rehearsal missing -> 708 OFF; rehearsal on file -> normal; only a pre-AMENDMENT-5 FAIL -> OFF)")


if __name__ == "__main__":
    main()
