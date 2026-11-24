#!/usr/bin/env python3
"""198 — has each hard-flagged closure already been revisited by a LATER experiment?

A closure that a later experiment has already corrected is not reopenable; it is
closed correctly and the map is self-healing. One that nothing later cites is the
real candidate. Counts forward citations of `NNN` in every later RESULTS.md plus
the two live docs.
"""
from __future__ import annotations
import json
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
W = Path(__file__).resolve().parent / "work"
hard = json.loads((W / "hardflags.json").read_text())
rows = {r["exp"]: r for r in json.loads((W / "enumeration.json").read_text())["rows"]}
TERMINAL = {"NO-GO", "VETO", "NEGATIVE", "CLOSED"}

docs = {p.name: p.read_text(errors="replace")
        for p in [ROOT / "AGENTS.md", ROOT / "docs/CURRENT_STATE.md",
                  ROOT / "docs/LOG.md", ROOT / "docs/BOARD-LEDGER.md"]}
res = {p.parent.name: p.read_text(errors="replace")
       for p in sorted((ROOT / "experiments").glob("*/RESULTS.md"))}
mandates = {p.parent.name: p.read_text(errors="replace")
            for p in sorted((ROOT / "experiments").glob("*/MANDATE.md"))}

OVERTURN = re.compile(
    r"overturn|refut|contradict|retract|falsif|superseded|no longer|"
    r"was wrong|is wrong|corrected|invalidat", re.I)

out = []
for o in hard:
    exp = o["exp"]
    if rows[exp]["verdict"] not in TERMINAL:
        continue
    num = exp.split("-")[0]
    n = int(num)
    pat = re.compile(rf"(?<!\d){num}(?!\d)")

    later, later_overturn = [], []
    for name, txt in list(res.items()) + list(mandates.items()):
        try:
            m = int(name.split("-")[0])
        except ValueError:
            continue
        if m <= n:
            continue
        for ln in txt.split("\n"):
            if pat.search(ln):
                later.append(name)
                if OVERTURN.search(ln):
                    later_overturn.append((name, ln.strip()[:170]))
                break
    doc_overturn = []
    for dn, txt in docs.items():
        for ln in txt.split("\n"):
            if pat.search(ln) and OVERTURN.search(ln):
                doc_overturn.append((dn, ln.strip()[:170]))

    out.append(dict(exp=exp, verdict=rows[exp]["verdict"],
                    endpoint=rows[exp]["primary_endpoint"],
                    surface=rows[exp]["primary_surface"],
                    flags=sorted(o["hard"]),
                    n_later_cites=len(set(later)),
                    later_overturn=later_overturn[:3],
                    doc_overturn=doc_overturn[:3]))

out.sort(key=lambda r: (r["n_later_cites"], len(r["later_overturn"])))
print(f"{'exp':32s} {'verdict':8s} {'endpt':7s} {'surface':13s} {'cites':>5s} "
      f"{'ovr':>3s}  flags")
for r in out:
    print(f"{r['exp']:32s} {r['verdict']:8s} {r['endpoint']:7s} {r['surface']:13s} "
          f"{r['n_later_cites']:5d} {len(r['later_overturn']) + len(r['doc_overturn']):3d}  "
          f"{','.join(f.replace('H_', '') for f in r['flags'])}")

print("\n=== flagged closures with the FEWEST later citations (candidates) ===")
for r in out[:8]:
    print(f"\n--- {r['exp']} ({r['verdict']}, {r['n_later_cites']} later cites)")
    for nm, ln in r["later_overturn"] + r["doc_overturn"]:
        print(f"    [{nm}] {ln}")

json.dump(out, open(W / "forward198.json", "w"), indent=1)
