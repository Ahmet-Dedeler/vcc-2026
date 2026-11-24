#!/usr/bin/env python3
"""618 stacks decision (PRED-618-stacks.md; coordinator's rule): overall-eq vs plain (price618 = price630's value_mean,
mean of the seeds) per arm and regime; MIN over {f 0.35, f 0.5}; the higher MIN wins, |diff| <= 0.0005 goes to (b).
-> work/stacks618.json"""
import json
from pathlib import Path
W = Path(__file__).resolve().parent / "work" / "stacks"
ARMS = {"(a) x1.25 + 630": ("t618a", "b_k100_q20r3p10"), "(b) 630 + DN": ("t618b", "b_k100_q20r3p10_D200"),
        "630 alone (x1)": ("t618u", "b_k100_q20r3p10"), "x1.25 alone": ("t618x", "x125")}
REG = ("1.0", "0.5", "0.35", "cmp")
out = {}
for arm, (name, v) in ARMS.items():
    p = W / f"price630_{name}.json"
    if not p.exists() or v not in json.loads(p.read_text()):
        out[arm] = None
        print(f"{arm}: not priced")
        continue
    rec = json.loads(p.read_text())[v]
    vals = {r: rec[r]["value_mean"] for r in REG if r in rec}
    out[arm] = dict(values=vals, seeds=rec.get("seeds"), min_05_035=min(vals["0.5"], vals["0.35"]),
                    ci={r: rec[r]["ci"] for r in REG if r in rec})
    print(f"{arm:18s} " + "  ".join(f"f {r} {vals[r]:+.5f}" for r in REG if r in vals)
          + f"  | MIN(0.5, 0.35) {out[arm]['min_05_035']:+.5f}  seeds {rec.get('seeds')}")
a, b = out.get("(a) x1.25 + 630"), out.get("(b) 630 + DN")
if a and b:
    d = a["min_05_035"] - b["min_05_035"]
    pick = "(a) x1.25 + 630" if d > 0.0005 else "(b) 630 + DN"
    h1 = "TRANSFERS" if 0.0004 <= d <= 0.0024 else "does NOT transfer"
    out["decision"] = dict(min_a_minus_min_b=d, pick=pick, rule="higher MIN over {f 0.35, f 0.5}; |diff| <= 0.0005 -> (b)",
                           registered_h1_prediction=+0.0014, h1_ranking=h1)
    print(f"MIN(a) - MIN(b) = {d:+.5f} -> PICK {pick} (registered H1 prediction +0.0014: H1's stack ranking {h1})")
(W.parent / "stacks618.json").write_text(json.dumps(out, indent=1))
