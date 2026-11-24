#!/usr/bin/env python3
"""198 — emit the enumeration as TSV + the hard-flagged subset as markdown."""
from __future__ import annotations
import json
from pathlib import Path
from collections import Counter

HERE = Path(__file__).resolve().parent
W = HERE / "work"
rows = json.loads((W / "enumeration.json").read_text())["rows"]
hard = {o["exp"]: o for o in json.loads((W / "hardflags.json").read_text())}

TERMINAL = {"NO-GO", "VETO", "NEGATIVE", "CLOSED"}

# full TSV
with open(W / "enumeration.tsv", "w") as f:
    f.write("exp\tverdict\tprimary_endpoint\tscored_members\tprimary_surface\t"
            "hard_flags\ttitle\n")
    for r in sorted(rows, key=lambda r: r["exp"]):
        hf = ",".join(sorted(hard.get(r["exp"], {}).get("hard", {}))) or "-"
        f.write(f"{r['exp']}\t{r['verdict']}\t{r['primary_endpoint']}\t"
                f"{','.join(r['scored_members']) or '-'}\t{r['primary_surface']}\t"
                f"{hf}\t{r['title'][:90]}\n")

term = [r for r in rows if r["verdict"] in TERMINAL]
tflag = [r for r in term if r["exp"] in hard]

print(f"total RESULTS.md            {len(rows)}")
print(f"terminal verdict            {len(term)}")
print(f"terminal AND hard-flagged   {len(tflag)}")
print(f"terminal, NO hard flag      {len(term) - len(tflag)}  <- survive the sweep untouched")
print()
print("verdicts:", dict(Counter(r["verdict"] for r in rows).most_common()))
print()
print("primary surface of TERMINAL verdicts:",
      dict(Counter(r["primary_surface"] for r in term).most_common()))
print()
print("primary endpoint of TERMINAL verdicts:",
      dict(Counter(r["primary_endpoint"] for r in term).most_common()))

lines = ["| exp | verdict | endpoint | surface | hard flags |",
         "|---|---|---|---|---|"]
for r in sorted(tflag, key=lambda r: r["exp"]):
    hf = ", ".join(k.replace("H_", "") for k in sorted(hard[r["exp"]]["hard"]))
    lines.append(f"| `{r['exp']}` | {r['verdict']} | `{r['primary_endpoint']}` | "
                 f"`{r['primary_surface']}` | {hf} |")
(W / "flagged.md").write_text("\n".join(lines) + "\n")
print(f"\nwrote work/enumeration.tsv ({len(rows)} rows) and work/flagged.md "
      f"({len(tflag)} rows)")
