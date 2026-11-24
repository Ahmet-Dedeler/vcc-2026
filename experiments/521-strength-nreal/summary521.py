"""521 — compact summary for RESULTS.md: base DE members per regime, and each headline policy's value on the three
emission seeds (seed-0 target-bootstrap CI) from policy521_s{0,1,2}.json.
    PYTHONDONTWRITEBYTECODE=1 uv run python experiments/521-strength-nreal/summary521.py
"""
import sys
sys.dont_write_bytecode = True
import json
import numpy as np
import lib521 as L
from policy521 import Table, means
T = Table(0)
out = {"base": {}, "rows": {}}
for reg in ("1.0", "0.5", "0.35", "cmp"):
    m = means(T.rows[f"{reg}|base"])
    out["base"][reg] = dict(zip(("fid", "jac", "reach", "nmae"), m.tolist()))
    print(f"base {reg:4s}: fid {m[0]:.5f} jac {m[1]:.5f} reach {m[2]:.5f} nmae {m[3]:.5f}")
P = [json.load(open(L.WORK / f"policy521_s{s}.json")) for s in (0, 1, 2)]
labels = ["ORACLE per-target best K", "NEG-ORACLE n<20 -> 25", "NULL-ORACLE n=0 -> 0", "S4 RANK q10 K25", "S4 RANK q20 K400",
          "S4 RANK q30 K400", "S4 TWO 10/25+20/400", "S4 BAYES (LOTO)", "S3 BAYES (LOTO)", "K562 BAYES (LOTO)"]
for var in ("b", "a"):
    print(f"\nvariant {var}: value s0 / s1 / s2  [seed-0 95 % target CI]  top3")
    for lab in labels:
        cells = []
        for reg in ("1.0", "0.5", "0.35", "cmp"):
            k = f"{reg}|headcpm|{var}|{lab}"
            v = [p[k]["value"] for p in P if k in p]
            ci = P[0][k]["ci"]; t3 = P[0][k]["top3"]
            out["rows"][f"{var}|{reg}|{lab}"] = dict(seeds=v, ci=ci, top3=t3, d=P[0][k]["d"])
            cells.append(f"{reg}: " + "/".join(f"{x:+.4f}" for x in v) + f" [{ci[0]:+.4f},{ci[1]:+.4f}] t3 {t3 if t3 == t3 else float('nan'):.2f}")
        print(f"  {lab:26s} " + " | ".join(cells))
(L.WORK / "summary521.json").write_text(json.dumps(out, indent=1, default=float))
