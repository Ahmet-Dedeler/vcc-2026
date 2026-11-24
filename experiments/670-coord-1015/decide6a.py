#!/usr/bin/env python3
"""The 6a replacement-entry rule (SIXA-DECISION.md + Amendments 1-5), applied mechanically.

Written 2026-10-20 by the coordinator [3c7a7e], before any D/E/F data. It decides nothing that the file does not already fix; it
only removes hand-reading on the day. Inputs:
  --decision708 <file>   708's day decision (experiments/708-def-panel-repricing/work/decision708-d708def.txt): its FINAL line
                         "FINAL    weights <hXcY> | rider ON|OFF | swap707 ON|OFF" (+ the RIDER line for a veto).
  board reads            bins from experiments/670-coord-1015/read695.py on the captured entries (R-ABC2H, C1, C0, W1, W2).
                         W1 / W2 default to the entries captured in submit-1021s1/s2.log; a missing read blocks nothing.
  --agg <json> ...       the lab's aggtag811 outputs (keys "<cand>|<base>|<TAG>|pds"); DEFo and DEF, pairs C1:B1, C0:B0, C1:C0.
  --lz PASS|FAIL|NONE    the LZ sign check on the deciding increment (NONE = not read yet).
  --lab-clean YES|NO     the lab's K / DE-H1 clean item (runbook §5 item 4).
  --x1-gates / --x2-gates PASS|FAIL|NOTBUILT   6a gate results, once built (default NOTBUILT).
Output: the swap branch, X1 and X2 with their env changes over R0's recipe and eligibility reasons, the build order, and what
ships given the gate results. Rule summary:
  X1 = R0 + 708's decided changes (weights move, rider veto, swap), minus blocked parts (W1 BELOW: H1 <= x2; W2 BELOW: CD4 = x1).
  X2 = the priced construction only (Amendment 5): C1 = R0 + swap + 811 or C0 = R0 + 811.
     swap OFF: C0 if C0:B0 passes (DEFo ci95[0] > 0 and (all folds > 0 or LT3 > 0); DEF point > 0) and C0's board read is not HARM.
     swap ON:  C1:B1 must pass the same bars and C1 not HARM; then C1 if C1:C0 on DEFo has ci95[0] > 0 (Amendment 4b),
               else C0 if C0:B0 passes and C0 not HARM, else C1.
     and the lab's clean item YES and LZ PASS.
  Ship: X2 if eligible and its gates PASS, else X1 if eligible and its gates PASS, else R0 stays final (Amendment 5's order).
usage: uv run python experiments/670-coord-1015/decide6a.py --decision708 <f> --agg <json> [<json> ...] --lz NONE --lab-clean YES
"""
from __future__ import annotations

import argparse
import json
import re
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
W = ROOT / "experiments/670-coord-1015/work"
BOARD = ROOT / "data/board-entries"
READ = ROOT / "experiments/670-coord-1015/read695.py"
RABC2, RABC2H, C1E, C0E = "nnHQfcZHsT3as9EafA3i", "ndt8swt7ILYn3qMJYZ9a", "RMxJtslBpY2GkwIMHOiB", "4rq3qHIuv8bHAxpTkUeu"
HOOK707 = "experiments/707-strong-tier-direction/dumphook707.sh"
HOOK_C1 = "experiments/811-lolo-model/dumphook811.sh"
HOOK_C0 = "experiments/811-lolo-model/dumphook811c0.sh"


def read695(*args):
    """Run read695.py; return its verdict text (first word is the bin)."""
    out = subprocess.run([sys.executable, str(READ), *args], capture_output=True, text=True, cwd=ROOT)
    if out.returncode != 0:
        raise SystemExit(f"read695 {' '.join(args)} failed: {out.stderr.strip()[-300:]}")
    return json.loads(out.stdout)["verdict"]


def captured(sl):
    """The entry id arm1015 captured for slot <sl>, if any."""
    p = W / f"submit-{sl}.log"
    if not p.exists():
        return None
    ids = re.findall(r"CAPTURED own entry ([A-Za-z0-9]{20})", p.read_text())
    return ids[-1] if ids else None


def board_bins(w1=None, w2=None):
    b = {"rabc2h": read695("--rabc2h", str(BOARD / f"{RABC2H}.json"), str(BOARD / f"{RABC2}.json")),
         "c1": read695("--c1", str(BOARD / f"{C1E}.json"), str(BOARD / f"{RABC2H}.json")),
         "c0": read695("--c0", str(BOARD / f"{C0E}.json"), str(BOARD / f"{RABC2}.json"))}
    w1 = w1 or captured("1021s1")
    w2 = w2 or captured("1021s2")
    b["w1"] = read695("--w1", str(BOARD / f"{w1}.json"), str(BOARD / f"{RABC2}.json")) if w1 else "NOT READ"
    b["w2"] = (read695("--w2", str(BOARD / f"{w2}.json"), str(BOARD / f"{w1}.json")) if (w1 and w2) else "NOT READ")
    return b


def parse708(path):
    t = Path(path).read_text()
    m = re.search(r"^FINAL\s+weights (\S+) \| rider (ON|OFF) \| swap707 (ON|OFF)", t, re.M)
    if not m:
        raise SystemExit(f"no FINAL line in {path}")
    rider_line = next((l for l in t.splitlines() if l.startswith("RIDER")), "")
    return dict(weights=m.group(1), rider=m.group(2) == "ON", vetoed="VETOED" in rider_line, swap=m.group(3) == "ON")


def load_agg(paths):
    res = {}
    for p in paths or []:
        res.update(json.loads(Path(p).read_text())["results"])
    return res


def bars(res, cand, base):
    """runbook §5 items 1-2 (Amendment 1): DEFo literal-pds ci95[0] > 0 and (all folds > 0 or LT3 > 0); DEF point > 0."""
    o, d = res.get(f"{cand}|{base}|DEFo|pds"), res.get(f"{cand}|{base}|DEF|pds")
    if o is None or d is None:
        return False, f"{cand}:{base} missing on {'DEFo' if o is None else 'DEF'}"
    ok_o = o["ci95"][0] > 0 and (o["folds_pos"] == o["n_folds"] or o["leave_top3_out"] > 0)
    ok_d = d["point"] > 0
    return ok_o and ok_d, (f"{cand}:{base} DEFo point {o['point']:+.5f} ci95[0] {o['ci95'][0]:+.5f} folds+ {o['folds_pos']}/{o['n_folds']} "
                           f"LT3 {o['leave_top3_out']:+.5f} -> {'ok' if ok_o else 'FAIL'}; DEF point {d['point']:+.5f} -> {'ok' if ok_d else 'FAIL'}")


def c1_over_c0(res):
    """Amendment 4b: the swap stays on top of 811 only if C1 - C0 on DEFo has ci95[0] > 0 (literal pds)."""
    o = res.get("C1|C0|DEFo|pds")
    if o is None:
        return None, "C1:C0 missing on DEFo (Amendment 4 needs it)"
    pos = sum(1 for f in o["per_fold"] if f["overall"] > 0)
    n = len(o["per_fold"])
    ok = o["ci95"][0] > 0
    return ok, (f"C1:C0 DEFo point {o['point']:+.5f} ci95[0] {o['ci95'][0]:+.5f}, folds+ {pos}/{n} -> "
                f"{'swap stays (clears its bound)' if ok else 'swap does not clear its bound on top of 811'}")


def decide(d708, bins, res, lz, clean, g1, g2):
    notes = []
    harm = lambda v: v.startswith("HARM") or v.startswith("BELOW")
    swap = d708["swap"] and not bins["rabc2h"].startswith("BELOW")
    if d708["swap"] and not swap:
        notes.append("708 turned the swap ON but R-ABC2H read BELOW: swap OFF")
    # X1
    env1 = {}
    m = re.match(r"h([0-9.]+)c([0-9.]+)$", d708["weights"])
    h, c = (float(m.group(1)), float(m.group(2))) if m else (2.0, 1.0)
    if h > 2.0 and bins["w1"].startswith("BELOW"):
        notes.append(f"W1 BELOW: H1 x{h:g} -> x2 (SIXA Amendment 2 backstop)"); h = 2.0
    if c > 1.0 and bins["w2"].startswith("BELOW"):
        notes.append(f"W2 BELOW: CD4 x{c:g} -> x1 (SIXA Amendment 3)"); c = 1.0
    if h != 2.0:
        env1["DEF_DAY_H1_STACK"] = f"x{h:g}"
    if c != 1.0:
        env1["DEF_DAY_CD4_STACK"] = f"x{c:g}"
    if d708["vetoed"]:
        env1["DEF_DAY_LODO_STRENGTH"] = "UNSET (708 vetoed 705's rider)"
    if swap:
        env1["DEF_DAY_DUMP_HOOK"] = HOOK707
    x1 = dict(eligible=bool(env1), env=env1, why="708 moved nothing" if not env1 else "708's decided changes")
    nested = not any(k in env1 for k in ("DEF_DAY_H1_STACK", "DEF_DAY_CD4_STACK", "DEF_DAY_LODO_STRENGTH"))
    # X2 (the priced construction only, Amendment 5)
    why2, cons = [], None
    c0ok, c0why = bars(res, "C0", "B0")
    c0board = not harm(bins["c0"])
    if not swap:
        why2.append(c0why)
        if not c0board:
            why2.append(f"C0's board read {bins['c0'][:20]}: swap-OFF X2 closed (Amendment 1)")
        cons = "C0" if (c0ok and c0board) else None
    else:
        c1ok, c1why = bars(res, "C1", "B1")
        why2.append(c1why)
        c1board = not harm(bins["c1"])
        if not c1board:
            why2.append(f"C1's board read {bins['c1'][:20]}: 811's swap-ON route closed")
        if c1ok and c1board:
            k, kwhy = c1_over_c0(res)
            why2.append(kwhy)
            if k is None:
                cons = None; why2.append("cannot apply Amendment 4 without C1:C0: X2 waits for it")
            elif k:
                cons = "C1"
            elif c0ok and c0board:
                cons = "C0"; why2.append(c0why + " -> X2 = C0 (Amendment 4)")
            else:
                cons = "C1"; why2.append(c0why + " -> C0 fails its own bars: X2 = C1 (Amendment 4 item 4)")
    if cons and clean != "YES":
        why2.append("the lab's K / DE-H1 clean item is not YES"); cons_ok = False
    else:
        cons_ok = cons is not None
    if cons_ok and lz != "PASS":
        why2.append(f"LZ sign check {lz}: " + ("X2 closed" if lz == "FAIL" else "not read yet (open LZ on the increment "
                                                                                + ("C1 - B1" if cons == "C1" else "C0 - B0") + ")"))
    x2_ok = cons_ok and lz == "PASS"
    x2 = dict(eligible=x2_ok, construction=cons, env={"DEF_DAY_DUMP_HOOK": HOOK_C1 if cons == "C1" else HOOK_C0} if cons else {},
              lz_increment=("C1 - B1" if cons == "C1" else "C0 - B0") if cons else None, why=why2)
    order = [n for n, x in (("X1", x1), ("X2", x2)) if x["eligible"]]
    if x2_ok and g2 == "PASS":
        ship = "X2"
    elif x2_ok and g2 == "NOTBUILT":
        ship = "PENDING: build " + " then ".join(n for n, g in (("X1", g1), ("X2", g2)) if n in order and g == "NOTBUILT")
    elif x1["eligible"] and g1 == "PASS":
        ship = "X1"
    elif x1["eligible"] and g1 == "NOTBUILT":
        ship = "PENDING: build X1"
    else:
        ship = "R0 (stays final)"
    return dict(swap_branch="ON" if swap else "OFF", weights=d708["weights"], rider_vetoed=d708["vetoed"], nested=nested,
                board=bins, X1=x1, X2=x2, build_order=order, gates=dict(X1=g1, X2=g2), ship=ship, notes=notes)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--decision708", required=True)
    ap.add_argument("--agg", nargs="*", default=[])
    ap.add_argument("--lz", choices=("PASS", "FAIL", "NONE"), default="NONE")
    ap.add_argument("--lab-clean", choices=("YES", "NO"), default="NO")
    ap.add_argument("--x1-gates", choices=("PASS", "FAIL", "NOTBUILT"), default="NOTBUILT")
    ap.add_argument("--x2-gates", choices=("PASS", "FAIL", "NOTBUILT"), default="NOTBUILT")
    ap.add_argument("--w1-entry")
    ap.add_argument("--w2-entry")
    a = ap.parse_args()
    out = decide(parse708(a.decision708), board_bins(a.w1_entry, a.w2_entry), load_agg(a.agg), a.lz, a.lab_clean,
                 a.x1_gates, a.x2_gates)
    print(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
