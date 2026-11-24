#!/usr/bin/env python3
"""617 (read-only): for every bullet of FINDINGS §2 "Closed", the experiments it cites and, per experiment, how often its
RESULTS.md names each evidence surface. A first pass only: the flagged closures are then read by hand (AUDIT.md).
Surfaces: thinned-Orion harness truths (the artefact, 616), PLOLO / donor-table folds (clean), H1, board, Broad, K562."""
import json
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
F = (ROOT / "docs/FINDINGS.md").read_text().splitlines()
a = next(i for i, l in enumerate(F) if l.startswith("## 2. Closed"))
b = next(i for i, l in enumerate(F) if i > a and l.startswith("# History"))
sec = F[a:b]
bullets, cur = [], None
for l in sec:
    if l.startswith("- "):
        cur = [l]
        bullets.append(cur)
    elif cur is not None and (l.startswith("  ") or l.strip() == ""):
        cur.append(l)
PAT = {
    "orion_thinned": r"orion[-_ ]?(hct116|hek293t)?[-_ ]?thinned|thinned (orion|harness)|orion harness|orion-hct116\.h5ad|depthmatched|gate-orion|cancer[- ]fold|hct116|hek293t",
    "plolo_donor": r"plolo|kolf2|donor[- ]table|october-like|arc-like|efp|efr\b",
    "h1": r"h1-test|h1_test|\bh1\b|h1 real|h1's",
    "board": r"board|BOARD-RESULT|entry ",
    "broad": r"broad",
    "k562_truth": r"k562[_ ]gwps fold|k562 fold|k562 harness|k562_c250|k562-gwps-c250|held-out[- ]k562",
}
exps = {p.name[:3]: p for p in (ROOT / "experiments").iterdir() if p.is_dir() and re.match(r"\d{3}-", p.name)}
out = []
for bl in bullets:
    text = " ".join(bl)
    nums = sorted(set(re.findall(r"\b(0\d\d|[1-6]\d\d)\b", text)))
    rng = re.findall(r"\b(\d{3})\s*[–-]\s*(\d{3})\b", text)
    for x, y in rng:
        if int(y) - int(x) < 60:
            nums += [f"{k:03d}" for k in range(int(x), int(y) + 1)]
    nums = sorted(set(n for n in nums if n in exps))
    rows = {}
    for n in nums:
        r = exps[n] / "RESULTS.md"
        t = r.read_text(errors="ignore").lower() if r.exists() else ""
        rows[n] = {k: len(re.findall(v, t)) for k, v in PAT.items()} | {"dir": exps[n].name, "has_results": r.exists()}
    out.append({"bullet": bl[0][:160], "experiments": rows})
(Path(__file__).parent / "work/scan617.json").write_text(json.dumps(out, indent=1))
for o in out:
    print("\n" + o["bullet"])
    for n, r in o["experiments"].items():
        flag = "ORION" if r["orion_thinned"] else ""
        print(f"   {n} {r['dir'][:38]:38s} orion {r['orion_thinned']:3d} plolo {r['plolo_donor']:3d} h1 {r['h1']:3d} board {r['board']:3d} broad {r['broad']:2d} k562 {r['k562_truth']:2d} {flag}")
