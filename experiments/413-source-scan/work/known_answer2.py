"""Known answer through overlap.py (literal and alias-resolved)."""
import json, pathlib
from overlap import report, ROOT
out = [report(d, [l.strip() for l in open(ROOT / f"data/donors/{d}/perts.txt")]) for d in ["k562_gwps", "orion_hct116", "orion_hek293t", "nadig_jurkat"]]
for o in out: print(o["source"], o["n_unique"], o["panel_exact"], o["panel_resolved"])
json.dump(out, open("known_answer_overlap.json", "w"), indent=1)
