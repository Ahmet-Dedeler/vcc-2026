#!/usr/bin/env python3
"""553 — the KA-H1 table for RESULTS.md from work/ka553.json (every row written by ka553.py).
usage: uv run python experiments/553-r7b-harness/katable553.py
"""
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
rows = json.loads((HERE / "work/ka553.json").read_text())
print("| surface | mode | seed | form | recorded run | arrays bit-equal (p, q, lfc, sums, jk, targets) | max \\|d\\| member (which) | pass |")
print("|---|---|---|---|---|---|---|---|")
for r in rows:
    eq = all(r["arrays_equal"].values())
    d = r["abs_d"]
    mk = max(d, key=d.get)
    rec = r["recorded"].split("experiments/")[-1]
    dm = "n/a (sealed stub scores nothing)" if r.get("surface") == "h1cb" else f"{d[mk]:.1e} ({mk})"
    print(f"| {r.get('surface', 'h1')} | `{r['mode']}` | {r['seed']} | {r['form']} | `{rec}` | {eq} | {dm} | "
          f"{'**PASS**' if r['pass_'] else 'FAIL'} |")
