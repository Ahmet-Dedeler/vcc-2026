"""413 known-answer row: panel overlap of held donor perturbation lists.

Expected (mandate): k562_gwps 269, orion_hct116 290, orion_hek293t 294, nadig_jurkat 0.
"""
import csv, json, pathlib
ROOT = pathlib.Path(__file__).resolve().parents[3]
panel = [r["target_gene"] for r in csv.DictReader(open(ROOT / "data/challenge/pert_counts.csv"))]
P = set(panel)
assert len(panel) == 300 and len(P) == 300
out = {}
for d in ["k562_gwps", "orion_hct116", "orion_hek293t", "nadig_jurkat"]:
    perts = [l.strip() for l in open(ROOT / f"data/donors/{d}/perts.txt") if l.strip()]
    out[d] = {"n_perts": len(perts), "n_unique": len(set(perts)), "panel_overlap": len(P & set(perts))}
print(json.dumps(out, indent=1))
json.dump(out, open(pathlib.Path(__file__).parent / "known_answer.json", "w"), indent=1)
